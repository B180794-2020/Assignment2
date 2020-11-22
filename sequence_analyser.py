#!/usr/bin/python3
import re
import subprocess
import readline
import pandas as pd
import os

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# Defining main body that will run the other functions in correct order/manner
def main():
    # Get user input of Taxonomy and Protein for search, choose to continue with search or not
    mysearch = user_search()
    progress, max_seq = fetch_data(mysearch)

    if progress == True:
        filename = fetch_fasta(mysearch)
        aligned, consensus, blastResults = conserved_sequence_analysis(filename, max_seq)
        accNumbers, top250 = Version_2_lot_top_250(filename, blastResults, aligned, 250)
        findMotifs(aligned, accNumbers)

# function for determining paramaters for user search
def user_search():
    # loop for giving an option to change search after input in case a mistake was made
    progress = False
    while progress == False:
        # user input for Protein family and Taxonomy
        family = input("Enter protein family: ")
        tax = input("Enter Taxanomic group: ")

        print(f"Protein family: {family}, Taxanomic group: {tax}")
        # choice to continue or not
        progress = yesNo("Are these correct? Y/N: ", "Please re-enter protein family and group.")
    # promting the user to see if they want to exclude partial and or predicted sequences from their analysis
    pred = part = ""
    ex_predict = yesNo("Do you wish to exclude predicted sequences? Y/N: ", "")
    ex_partial = yesNo("Do you wish to exclude partial sequences? Y/N: ", "")
    if ex_predict == True:
        pred = "NOT predicted"
    if ex_partial == True:
        part = "NOT partial"
    # Term that will be searched for on NCBI
    mysearch = f"{tax}[Organism] AND {family}[Protein] {pred} {part}"

    return mysearch


def fetch_data(mysearch):
    # calling shell command of esearch with specified paramaters and piping results into grep that selects titles of each result
    res = subprocess.check_output(f"esearch -db protein -query \"{mysearch}\" | efetch -format docsum | grep -E \"<AccessionVersion>|<Title>\" ", shell=True)
    # find [species names] that are in square brackets and accession numbers
    species = re.finditer(r"\[.*?\]", str(res))
    accession = re.finditer(r'<AccessionVersion>.*?</AccessionVersion>', str(res))
    # put species into a list wihout the brackets
    speciesList = []
    acessionList = []

    for i in species:
        speciesList.append(i.group(0).strip("[]"))
    for i in accession:
        acessionList.append(i.group(0).strip("<AccessionVersion></AccessionVersion>"))
    # total number of results in the list and number of unique species names
    totalResults = len(speciesList)
    speciesNumber = len(set(speciesList))
    # conditions for continuing process
    max_sequences = 1000
    max_species = 500
    progress = True
    # prompt user to continue based on n of seq and species
    print(f"Number of Sequences: {totalResults}\nNumber of Species: {speciesNumber}")

    if totalResults > max_sequences:
        progress = yesNo("Warning! Search resulted in more than 1000 sequences. \n do you wish to continue? Y/N: ",
                         "Exiting")
    if speciesNumber > max_species:
        progress = yesNo("Warning! Search resulted in more than 1000 sequences. \n do you wish to continue? Y/N: ",
                         "Exiting")

    return progress, max_sequences


def fetch_fasta(mysearch):
    # let user choose file name where fasta will be saved
    filename = input("Enter filename: ")
    # get fastafile
    subprocess.call(f"esearch -db protein -query \"{mysearch}\" | efetch -format fasta > {filename}", shell=True)
    keep = yesNo("Do you wish to remove duplicate sequences? This will also remove isoforms of the same protein. Y/N: ",
                 "")
    if keep == True:
        out = filename + ".keep"
        # EMBOSS skip redundant to identify duplicate sequences, keeps longer of two if identical
        subprocess.call(
            f"skipredundant -maxthreshold 100.0 -minthreshold 100.0 -mode 2 -gapopen 0.0 -gapextend 0.0 -outseq {out} -datafile EBLOSUM62 -redundant \"\" {filename}",
            shell=True)
        return out
    else:
        return filename


def conserved_sequence_analysis(filename, max_seq):
    # Names of files that will be created
    aligned = filename + ".out"
    consensus = filename + ".con"
    blastResults = filename + ".blast"

    # sequence alignment and finding a consensus sequence
    subprocess.call(f"clustalo --force --threads 8 --maxnumseq {max_seq} -i {filename} -o {aligned}", shell=True)
    subprocess.call(f"cons -datafile EBLOSUM62 -sequence {aligned} -outseq {consensus}", shell=True)

    # BLAST
    subprocess.call(f"makeblastdb -in {filename} -dbtype prot -out {filename}", shell=True)
    subprocess.call(f"blastp -db {filename} -query {consensus} -outfmt 7 > {blastResults}", shell=True)

    return aligned, consensus, blastResults


def plot_top_250(blastResults, n):
    # setting headings for dataframe, assumes blast with n rows and 12 columns (-outfmt 7)
    headings = ["queryacc.", "subjectacc.", "% identity", "alignment_length",
                "mismatches", "gap_opens", "q.start", "q.end", "s.start",
                "s.end", "e-value", "bit_score"]
    # setting up dataframe using pandas
    df = pd.read_csv(f"{blastResults}", skiprows=5, names=headings, sep="\t")
    # sorting according to bitscores
    df.sort_values('bit_score', ascending=False, inplace=True)

    # taking top n number of sequences
    max_seq = n
    dfsubset = df[0:max_seq]

    # collecting accession numbers of the top 250
    accNumbers = dfsubset["subjectacc."].tolist()
    if len(accNumbers) < n:
        accNumbers = accNumbers[:-1]

    # preparing for a new seach of just the top 250
    mysearch = ' '.join(accNumbers)
    filename = "top250"
    aligned = filename + ".aligned"

    # searching for the top 250, aligning them and plotting the conservation using EMBOSS plotcon
    subprocess.call(f"esearch -db protein -query \"{mysearch}\" | efetch -format fasta > {filename}", shell=True)
    subprocess.call(f"clustalo --force --threads 8 --maxnumseq {max_seq} -i {filename} -o {aligned}", shell=True)
    subprocess.call(f" plotcon -winsize 4 -graph x11  {aligned}", shell=True)

    return aligned


def Version_2_lot_top_250(filename, blastResults, aligned, n):
    # setting headings for dataframe, assumes blast with n rows and 12 columns (-outfmt 7)
    headings = ["queryacc.", "subjectacc.", "% identity", "alignment_length",
                "mismatches", "gap_opens", "q.start", "q.end", "s.start",
                "s.end", "e-value", "bit_score"]
    # setting up dataframe using pandas
    df = pd.read_csv(f"{blastResults}", skiprows=5, names=headings, sep="\t")
    # sorting according to bitscores
    df.sort_values('bit_score', ascending=False, inplace=True)

    # taking top n number of sequences
    max_seq = n
    dfsubset = df[0:max_seq]

    # collecting accession numbers of the top 250
    accNumbers = dfsubset["subjectacc."].tolist()
    if len(accNumbers) < n:
        accNumbers = accNumbers[:-1]

    top250 = filename + ".250"
    topFasta = top250 + ".fasta"

    # automatically closes after loop
    with open(top250, 'w') as f:
        for num in accNumbers:
            f.write(f"{num}\n")

    # preparing for a new seach of just the top 250

    subprocess.call(f"/localdisk/data/BPSM/Assignment2/pullseq -i {aligned}  -n > {topFasta}", shell = True )

    # searching for the top 250, aligning them and plotting the conservation using EMBOSS plotcon
    subprocess.call(f"plotcon -winsize 4 -graph x11  {topFasta}", shell=True)

    return accNumbers, top250

def findMotifs(aligned, accnumbers):
    nameList= []

    for number in accnumbers:
        motifs = number + ".motif"
        subprocess.call(f"/localdisk/data/BPSM/Assignment2/pullseq -i {aligned}  -n > {number}", shell = True)
        subprocess.call(f"patmatmotifs {number} -outfile {motifs}", shell = True)       
        subprocess.call(f"rm {number}", shell = True)
        nameList.append(motifs)

    myDic = {}

    for name in nameList:
        with open(name, "r") as f:
            lines = f.readlines()
            for line in lines:
                if "length" in line:
                    myDic[name]["length"] = line.strip("length = ")
    print(myDic)
    return None

def yesNo(question, reprompt):
    yes = ["y", "Y", "Yes", "YES", "yes"]
    no = ["n", "N", "No", "NO", "no"]

    invalid = False
    while invalid == False:
        answer = input(question)
        if answer in yes:
            return True
        elif answer in no:
            print(reprompt)
            return False

        print("Invalid input. Please answer Yes or No")

if __name__ == '__main__':
    main()

