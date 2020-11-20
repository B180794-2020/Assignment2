#!/usr/bin/python3
import subprocess
import re

def main():
    mySearch, progress = user_search()
    if progress == True:
        print(mySearch)
        
def user_search():
    happy = False
    while happy == False:
        family = input("Enter protein family: ")
        tax = input("Enter Taxanomic group: ")

        print(f"Protein family: {family}, Taxanomic group: {tax}")
        happy = yesNo("Are these correct? Y/N:", "Please re-enter protein family and group.")

    pred = part = ""
    ex_predict = yesNo("Do you wish to exclude predictive sequences? Y/N: ", "")
    ex_partial = yesNo("Do you wish to exclude partial sequences? Y/N: ", "")

    if ex_predict == True:
        pred = "NOT predicted"

    if ex_partial == True:
        part = "NOT partial"

    mySearch = "Aves[Organism] AND glucosE-6-phosPhatase[Protein Family] NOT predicted NOT partial"
    res = subprocess.check_output(f"esearch -db protein -query \"{mySearch}\" | efetch -format docsum | grep \"<Title>\" ", shell=True)

    species = re.finditer(r'\[.*?\]', str(res))

    speciesList = []
    for i in species:
        speciesList.append(i.group(0).strip("[]"))

    totalResults = len(speciesList)
    speciesNumber = len(set(speciesList))
    
    max_sequences = 1000
    max_species = 500
    progress = True

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
