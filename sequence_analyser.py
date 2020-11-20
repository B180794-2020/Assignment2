#!/usr/bin/python3
import subprocess
import re

#Defining main body that will run the other functions in correct order/manner
def main():
    #Get user input of Taxonomy and Protein for search, choose to continue with search or not
    mySearch, progress = user_search()
    if progress == True:
        print(mySearch)
#function for determining paramaters for user search
def user_search():
    #loop for giving an option to change search after input in case a mistake was made

    progress = False
    while progress == False:
        #user input for Protein family and Taxonomy
        family = input("Enter protein family: ")
        tax = input("Enter Taxanomic group: ")

        print(f"Protein family: {family}, Taxanomic group: {tax}")
        #choice to continue or not
        progress = yesNo("Are these correct? Y/N:", "Please re-enter protein family and group.")

    #promting the user to see if they want to exclude partial and or predicted sequences from their analysis
    pred = part = ""
    ex_predict = yesNo("Do you wish to exclude predicted sequences? Y/N: ", "")
    ex_partial = yesNo("Do you wish to exclude partial sequences? Y/N: ", "")

    if ex_predict == True:
        pred = "NOT predicted"

    if ex_partial == True:
        part = "NOT partial"

    #Term that will be searched for on NCBI
    mySearch = "Aves[Organism] AND glucosE-6-phosPhatase[Protein Family] NOT predicted NOT partial"
    #calling shell command of esearch with specified paramaters and piping results into grep that selects titles of each result
    res = subprocess.check_output(f"esearch -db protein -query \"{mySearch}\" | efetch -format docsum | grep \"<Title>\" ", shell=True)
    
    #find species names that are in square brackets
    species = re.finditer(r'\[.*?\]', str(res))
    
    #put species into a list wihout the brackets
    speciesList = []
    for i in species:
        speciesList.append(i.group(0).strip("[]"))
    
    #total number of results in the list and number of unique species names
    totalResults = len(speciesList)
    speciesNumber = len(set(speciesList))

    #conditions for continuing process
    max_sequences = 1000
    max_species = 500
    progress = True
    
    #prompt user to continue based on n of seq and species
    print(f"Number of Sequences: {totalResults}\nNumber of Species: {speciesNumber}" )
    if totalResults > max_sequences:
        progress = yesNo("Warning! Search resulted in more than 1000 sequences. \n do you wish to continue? Y/N: ", "Exiting")
    if speciesNumber > max_species:
        progress = yesNo("Warning! Search resulted in more than 1000 sequences. \n do you wish to continue? Y/N: ", "Exiting")

    return mySearch, progress

def fetch_data(mySearch):

    return None





def yesNo(question,reprompt):
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
