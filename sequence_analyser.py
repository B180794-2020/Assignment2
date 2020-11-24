#!/usr/bin/python3
import os
import re
import string
import subprocess

import pandas as pd

# important to import readline despite none of its functions being spesifically called:
import readline
# enables the use of backspace in the terminal without them being part of user inputs
# Finding path to file location and directory, changing working directory to location
abspath = os.path.abspath(__file__)
mydir = os.path.dirname(abspath)
os.chdir(mydir)


# Defining main body that will run the other functions in correct order/manner
def main():
    # defining variables to be used by functions
    new_search = True
    threads = 16
    number_of_top = 250
    max_sequences = 1000
    max_species = 300

    # fields to limit search results
    field_one = "[Protein Family]"
    field_two = "[Organism]"
    pullseq = "/localdisk/data/BPSM/Assignment2/pullseq"

    while new_search:
        # Get user input of Taxonomy and Protein for search, choose to continue with search or not
        mysearch, progress = user_search(field_one, field_two)

        if progress:
            # Conduct search and fetch number of results, evaluate is there the required minimum of 3 seq
            progress, max_seq, number_of_results = fetch_data(mysearch, max_sequences, max_species)

            # if there are enough sequences and the user wishes to continue
            if progress:
                # fetch the fasta files
                filename = fetch_fasta(mysearch, number_of_results)
                # align the sequences, find consensus sequence and conduct blast search using them
                # generates blast results, the consensus sequence and aligned sequence fastas
                aligned, consensus, blast_results = conserved_sequence_analysis(filename, max_seq, threads)
                # plots the conservation of the top desired number of results, by default 250. Uses bit scores to sort top
                accnumbers, save = plot_top_250(filename, blast_results, aligned, number_of_top, pullseq)
                # finds known motifs from the PROSITE database based on the aligned sequences
                find_motifs(aligned, accnumbers, pullseq)
                wildcard(filename, aligned, accnumbers, pullseq)
            else:
                print("Search has been cancelled")
        else:
            print("Search has been cancelled")

        new_search = yes_no("\033[1;32;40m Would you like to do another search? Y/N ", "Exiting")
        print("\033[1;37;40m ")


# function for determining parameters for user search. Takes two arguments: two [fields] used on NCBI in searches
# Returns a string mysearch containing the search query for esearch
def user_search(field_one, field_two):
    # declaring variables in case somehow exit loop without declaring them
    input_one = input_two = ""
    # loop for giving an option to change search after input in case a mistake was made
    progress = False
    while not progress:
        # user input for Protein family and Taxonomy
        # if keyboard is interrupted loop is broke, one or the other == "", progress is therefore False
        # by default getting a input for the protein family
        try:
            input_one = input(f"\033[1;32;40m Enter your desired {field_one}: ")
            print("\033[1;37;40m ")
            if not check_input(input_one):
                return "", False

        except KeyboardInterrupt:
            print("Error: Keyboardinterrupt")
            return "", False

        # by default getting an input for the taxonomy
        try:
            input_two = input(f"\033[1;32;40m Enter {field_two}: ")
            print("\033[1;37;40m ")
            if not check_input(input_two):
                return "", False

        except KeyboardInterrupt:
            print("Error: Keyboardinterrupt")
            return "", False

        print(f"Protein family: {input_one}, Taxonomic group: {input_two}")
        # choice to continue or not
        progress = yes_no("Are these correct? Y/N: ", "Please re-enter protein family and group.")

        # test to see if blank parameters
        if input_one.strip(" ") == "" or input_two.strip(" ") == "":
            print("Warning one of the fields was left blank. Please re-enter the Protein family and Taxonomic group.")
            progress = False
    # prompting the user to see if they want to exclude partial and or predicted sequences from their analysis

    pred = part = ""
    ex_predict = yes_no("Do you wish to exclude predicted sequences? Y/N: ", "")
    ex_partial = yes_no("Do you wish to exclude partial sequences? Y/N: ", "")

    if ex_predict:
        pred = "NOT predicted NOT hypothetical"
    if ex_partial:
        part = "NOT partial"
    # Term that will be searched for on NCBI
    mysearch = f"{input_two + ' ' + field_two} AND {input_one + ' ' + field_one} {pred} {part}"

    # return the search term
    return mysearch, progress


# function for fetching the data of how many search results there are for the search.
# takes a string containing the query as an argument. Returns whether the search should progress (boolean)
# maximum number of sequences (int) and number of total results (int)
def fetch_data(mysearch, max_sequences, max_species):
    # calling shell command of esearch with specified parameters
    # piping results into grep that selects titles of each result
    try:
        print("Conducting search...")
        res = subprocess.check_output(
            f"esearch -db protein -query \"{mysearch}\" | efetch -format docsum | grep -E \"<AccessionVersion>|<Title>\" ",
            shell=True)
    except subprocess.CalledProcessError:
        print("Error: There were no results for search")
        progress = False
        return progress, 0, 0

    # find [species names] that are in square brackets and accession numbers
    species = re.finditer(r"\[.*?\]", str(res))
    accession = re.finditer(r'<AccessionVersion>.*?</AccessionVersion>', str(res))
    # put species into a list without the brackets
    species_list = []
    acession_list = []

    for i in species:
        species_list.append(i.group(0).strip("[]"))
    for i in accession:
        acession_list.append(i.group(0).strip("<AccessionVersion></AccessionVersion>"))
    # total number of results in the list and number of unique species names
    total_results = len(species_list)
    species_number = len(set(species_list))
    # conditions for continuing process

    progress = True
    # prompt user to continue based on n of seq and species
    print(f"Number of Sequences: {total_results}\nNumber of Species: {species_number}")

    if total_results > max_sequences:
        progress = yes_no(
            f"Warning! Search resulted in more than {max_sequences} sequences. \n do you wish to continue? Y/N: ",
            "Exiting")
        max_sequences = total_results
    if species_number > max_species:
        progress = yes_no(
            f"Warning! Search resulted in more than {max_species} species. \n do you wish to continue? Y/N: ",
            "Exiting")
    if total_results < 3:
        print("Not enough sequences in search result to conduct analysis")
        progress = False

    return progress, max_sequences, total_results


# Full esearch that fetches fasta files for each result, defines a filename that is used throughout the analysis
# Takes a search query (str) and number of results for the search as arguments.
# Option to exclude duplicate sequences by calling EMBOSS con. Returns a filename of the fasta files
def fetch_fasta(mysearch, number_of_results):
    # let user choose file name where fasta will be saved
    while True:
        filename = input("\033[1;32;40m Enter filename: ").lower()
        print("\033[1;37;40m ")
        if check_input(filename):
            break

    # get fasta files
    print("Fetching Fasta files from NCBI protein database...")
    subprocess.call(f"esearch -db protein -query \"{mysearch}\" | efetch -format fasta > {filename}", shell=True)
    remove = False
    # Only gives option to remove seq if there are more than 3 of them
    if number_of_results > 3:
        remove = yes_no(
            "Do you wish to remove duplicate sequences? This will also remove isoforms of the same protein. Y/N: ", "")
    if remove:
        out = filename + ".keep"
        # EMBOSS skip redundant to identify duplicate sequences, keeps longer of two if identical
        print("Removing redundant...")
        subprocess.call(
            f"skipredundant -maxthreshold 100.0 -minthreshold 100.0 -mode 2 -gapopen 0.0 -gapextend 0.0 -outseq {out} -datafile EBLOSUM62 -redundant \"\" {filename}",
            shell=True)
        # Counting the number of sequences in the .keep file as well as creating the file
        print("Checking remaining files...")
        with open(out) as f:
            sequence_count = 0
            #list containing each line
            lines = f.readlines()
            for line in lines:
                # each line that starts with a > should indicate a sequence
                if ">" in line:
                    sequence_count += 1
        # checking that there are at least 3
        print(f"Remaining files: {sequence_count}")
        if sequence_count < 3:
            print("Not enough sequences after removing redundant sequences.")
            print("Reverting to using original full list of sequences.")
            subprocess.call(f"rm {out}", shell=True)
            return filename

        return out
    else:
        return filename


# Function for conducting sequence alignment, blast database and blast search constructing a consensus sequence
# Takes a filename (str) max number of sequences and number of threads to be used as arguments
# Outputs 3 filenames for aligned, consensus and blast results respectively
def conserved_sequence_analysis(filename, max_seq, threads):
    # Names of files that will be created
    aligned = filename + ".aligned"
    consensus = filename + ".con"
    blast_results = filename + ".blast"

    # sequence alignment and finding a consensus sequence
    print("Aligning sequences...")
    subprocess.call(f"clustalo --force --threads {threads} --maxnumseq {max_seq} -i {filename} -o {aligned}",
                    shell=True)
    print("Finding consensus sequence...")
    subprocess.call(f"cons -datafile EBLOSUM62 -sequence {aligned} -outseq {consensus}", shell=True)

    # BLAST
    print("Constructing Blast database...")
    subprocess.call(f"makeblastdb -in {filename} -dbtype prot -out {filename}", shell=True)
    print("Conducting blast search with consensus sequence...")
    subprocess.call(f"blastp -db {filename} -query {consensus} -outfmt 7 > {blast_results}", shell=True)

    return aligned, consensus, blast_results


# Function for plotting the top results of based on the bit scores of the blast search
# Takes the filename of the blast results and aligned sequences and the number of sequences required as arguments
# Returns the filename of the saved image as save.ps and a list: accession numbers of top bit scores in order
def plot_top_250(filename, blast_results, aligned, n, pullseq):
    # setting headings for dataframe, assumes blast with n rows and 12 columns (-outfmt 7)
    headings = ["queryacc.", "subjectacc.", "% identity", "alignment_length",
                "mismatches", "gap_opens", "q.start", "q.end", "s.start",
                "s.end", "e-value", "bit_score"]
    # setting up dataframe using pandas
    df = pd.read_csv(f"{blast_results}", skiprows=5, names=headings, sep="\t")
    # sorting according to bitscores
    df.sort_values('bit_score', ascending=False, inplace=True)

    # taking top n number of sequences
    max_seq = n
    dfsubset = df[0:max_seq]
    print(dfsubset)

    print(f"Finding top {n} for plotting sequence conservation...")
    # collecting accession numbers of the top 250
    acc_numbers = dfsubset["subjectacc."].tolist()

    # removing nan in case there are less than n sequences in the initial test
    if len(acc_numbers) < n:
        acc_numbers = acc_numbers[:-1]

    # removing weird characters from acc numbers
    for i in range(len(acc_numbers)):
        if "|" in acc_numbers[i]:
            res = re.search('\|(.*)\|', acc_numbers[i])
            acc = str(res.group(1))
            acc_numbers[i] = acc

    acc_numbers = set(acc_numbers)

    # filename
    top250 = filename + ".250"
    top_fasta = top250 + ".fasta"

    # creating a file containing the accnum of top 250
    with open(top250, "w") as f:
        for num in acc_numbers:
            f.write(f"{num}\n")

    # preparing for a new seach of just the top 250
    print(f"Fetching top {n} for plotting sequence conservation...")
    subprocess.call(f"{pullseq} -i {aligned} -n {top250}> {top_fasta}", shell=True)

    while True:
        save = input("\033[1;32;40m Choose filename for saving Conservation plot: ") + ".plot"
        print("\033[1;37;40m ")
        if check_input(save):
            break

    # searching for the top 250, aligning them and plotting the conservation using EMBOSS plotcon
    subprocess.call(f"plotcon -goutfile {save} -scorefile EBLOSUM62 -winsize 4 -graph ps {top_fasta}", shell=True)

    # removing temporary files that were generated for each to be taken as an input
    subprocess.call(f"rm {top250}", shell=True)
    subprocess.call(f"rm {top_fasta}", shell=True)

    # returns a list containing the top accnumbers
    return acc_numbers, (save + ".ps")


# Function to search the PROSITE database for motifs for each of the top sequences
# Requires the filename containing the fasta sequences and a list containing the accession numbers as arguments
# Creates a table of the results that is saved
def find_motifs(aligned, acc_numbers, pullseq):
    print(f"Seaching for protein motifs...")
    motif_list = []

    # Blast may have more than one alignment for each sequence multiple times
    # Only need to search for them once -> set contains unique elements
    acc_numbers = set(acc_numbers)

    # creating temp files for input into pullseq containing the sequence accnumber
    for number in acc_numbers:
        with open(number, "w") as f:
            f.write(f"{number}")

        # file names for patmatmotif output and temp
        motifs = number + ".motif"
        temp = "temp"

        # finding the fasta of a protein based on accesion number and storing it in a temporary file
        subprocess.call(f"{pullseq} -i {aligned} -n {number} > {temp}", shell=True)
        subprocess.call(f"patmatmotifs {temp} -outfile {motifs}", shell=True)

        # removing temporary files that were generated for each to be taken as an imput
        subprocess.call(f"rm {number}", shell=True)
        subprocess.call(f"rm {temp}", shell=True)
        motif_list.append(motifs)

    my_dic = {}
    # Constructing a dictionary containing the information from patmatmotifs
    for motif in motif_list:
        # open report file of patmat, scan each line for key words, add only value
        try:
            with open(motif) as f:
                lines = f.readlines()
                # If no motifs then return as empty strings
                mot = length = start = end = ""
                for line in lines:
                    if "Length" in line:
                        length = int(line.split(" ")[2].strip("\n"))
                    if "Start" in line:
                        start = int(line.split(" ")[3])
                    if "End" in line:
                        end = int(line.split(" ")[3])
                    if "Motif" in line:
                        mot = line.split(" ")[2].strip("\n")
                my_dic[motif] = [mot, length, start, end]

            # remove the report file after information has been extracted
            subprocess.call(f"rm {motif}", shell=True)
        except FileNotFoundError:
            print(f"Error: Pullseq was unable to find Fasta for accession number {motif}")
            print(f"Returned empty file. Accession number is excluded from search.")

    headings = ["Motif", "Length", "Start", "End"]
    df = pd.DataFrame.from_dict(my_dic, orient='index', columns=headings)
    print(df)

    while True:
        save = input("\033[1;32;40m Input filename for motifs: ")
        print("\033[1;37;40m ")
        if check_input(save):
            break

    with open(save, "w") as f:
        df.to_string(f)
    return True


# protein function analysis
def wildcard(filename, aligned, acc_numbers, pullseq):

    # filename
    top250 = filename + ".250"
    top_fasta = top250 + ".fasta"

    # creating a file containing the accnum of top 250
    with open(top250, "w") as f:
        for num in acc_numbers:
            f.write(f"{num}\n")

    # preparing for a new seach of just the top 250
    subprocess.call(f"{pullseq} -i {aligned} -n {top250}> {top_fasta}", shell=True)

    while True:
        save = input("\033[1;32;40m Choose filename for Predicted transmembrane report: ")

        print("\033[1;37;40m ")
        if check_input(save):
            break

    subprocess.call(f"tmap {top_fasta} -out {save} -goutfile {save + '.plot'} -graph ps", shell=True)
    subprocess.call(f"rm {top250}", shell=True)
    subprocess.call(f"rm {top_fasta}", shell=True)


# function for determining True/False based on user input
# takes two strings as arguments. The first is displayed as the prompt, the other is False is chosen
def yes_no(question, reprompt):
    # list of possible approved answers for each option
    # For automated convenienve a function returning true or false could also be used as str(True)
    yes = ["y", "Y", "Yes", "YES", "yes", "True"]
    no = ["n", "N", "No", "NO", "no", "False"]

    invalid = False

    # infinite loop if an answer other than one in the list is given
    while not invalid:
        while True:
            try:
                answer = input(f"\033[1;32;40m {question}")
                print("\033[1;37;40m ")
                break
            except KeyboardInterrupt:
                print("Error: KeyboardInterrupr")
                print("Try again: ")

        # returns true or false depending on answer
        if answer in yes:
            return True
        elif answer in no:
            print(reprompt)
            return False

        print("Invalid input. Please answer Yes or No")


def check_input(user_input):
    accepted_specials = [" ", ".", "-"]

    for char in user_input:
        if char not in string.digits and char not in string.ascii_letters and char not in accepted_specials:
            print('Special characters are not allowed')
            return False

    return True


if __name__ == '__main__':
    main()

