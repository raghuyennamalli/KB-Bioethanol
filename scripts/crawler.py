from Bio import Entrez as en
import json
import re
import mongoConn as mc

#Change this value to our email ID
en.email = "122013041@sastra.ac.in"

#This variable is a list of your queries
queries = []

user_input = input("Enter query:")
queries.append(user_input)

#Outer loop, to obtain matching hits from GEO for each query
metadata = []
for query in queries:
    handle = en.esearch(db="gds",term=query,retmax=10000)
    result = en.read(handle)
    handle.close()

    print (len(result['IdList'])," hits \n")

    #Inner loop, to read metadata for each hit from GEO
    for ID in result['IdList']:
        handle2 = en.esummary(db="gds",id=ID,rettype='xml',retmode='text')
        details = en.read(handle2)[0]
        metadata.append(details)
        handle2.close()

#Write metadata to a JSON file
with open("metadata.json","w") as file:
    json.dump(list(metadata),file)

#Insert metadata into database
mc.write_metadata('metadata.json')

print ("\n\nSuccessfully completed !!")
from datetime import datetime
print (datetime.now().strftime('%d-%m-%Y %H:%M:%S'))
