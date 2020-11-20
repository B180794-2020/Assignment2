#!/usr/bin/python3
import subprocess
import re

#Main body that will run the functions
def main():
    mySearch, count, con = user_search()
    print( mySearch, count, con)

#Get initial user search parameters for NCBI protein database
def user_search():
    happy = False
    while happy == False:
        family = input("Enter protein family: ")
        tax = input("Enter Taxanomic group: ")

        print(f"Protein family: {family}, Taxanomic group: {tax}")
        happy = yesNo("Are these correct? Y/N:", "Please re-enter protein family and group.")

    pred = part = ""
    ex_predict = yesNo("Do you wish to exclude predictive sequences? Y/N: ","")
    ex_partial = yesNo("Do you wish to exclude partial sequences? Y/N: ","")

    if ex_predict == True:
        pred = "NOT predicted"

    if ex_partial == True:
        part = "NOT partial"

    mySearch = f"{tax}[Organism] AND {family}[Protein Family] {pred} {part}"

    res = str(subprocess.check_output(f"esearch -db protein -query \"{mySearch}\" | efetch -format docsum | grep \"<Title>\" | less", shell = True))
    count = get_species(res)
    print (count)

    return mySearch, count, con
#get total number of results and number of different species
def get_species(res):
    new = re.sub("\[(.*?)\]", res)
    return count

#fetch the actual data?
def fetch_data(mySearch):



    return None




#function to get user choice returns True or False
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
#run
if __name__ == '__main__':
    main()
