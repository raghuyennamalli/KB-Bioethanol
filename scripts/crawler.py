from Bio import Entrez as en
import json
from datetime import datetime
import re

en.email = "122013041@sastra.ac.in"

#This variable's value should be your query
queries = ['''("Neurospora crassa"[Organism] OR neurospora crassa[All Fields]) AND ("glucose"[MeSH Terms] OR glucose[All Fields]) AND "Expression profiling by high throughput sequencing"[Filter]''',
         '''("Neurospora crassa"[Organism] OR neurospora crassa[All Fields]) AND ("cellulose"[MeSH Terms] OR cellulose[All Fields]) AND "Expression profiling by high throughput sequencing"[Filter]''',
         '''("Neurospora crassa"[Organism] OR neurospora crassa[All Fields]) AND ("cellulose"[MeSH Terms] OR avicel[All Fields]) AND "Expression profiling by high throughput sequencing"[Filter]''',
         '''("Neurospora crassa"[Organism] OR neurospora crassa[All Fields]) AND ("xylose"[MeSH Terms] OR xylose[All Fields]) AND "Expression profiling by high throughput sequencing"[Filter]''',
         '''("Neurospora crassa"[Organism] OR neurospora crassa[All Fields]) AND ("carbon"[MeSH Terms] OR carbon[All Fields]) AND "Expression profiling by high throughput sequencing"[Filter]''',
         '''("aspergillus fumigatus"[MeSH Terms] OR "Aspergillus fumigatus"[Organism] OR Aspergillus fumigatus[All Fields]) AND ("glucose"[MeSH Terms] OR glucose[All Fields]) AND "Expression profiling by high throughput sequencing"[Filter]''',
         '''("aspergillus fumigatus"[MeSH Terms] OR "Aspergillus fumigatus"[Organism] OR Aspergillus fumigatus[All Fields]) AND ("cellulose"[MeSH Terms] OR cellulose[All Fields]) AND "Expression profiling by high throughput sequencing"[Filter]''',
         '''("aspergillus fumigatus"[MeSH Terms] OR "Aspergillus fumigatus"[Organism] OR Aspergillus fumigatus[All Fields]) AND ("cellulose"[MeSH Terms] OR avicel[All Fields]) AND "Expression profiling by high throughput sequencing"[Filter]''',
         '''("aspergillus fumigatus"[MeSH Terms] OR "Aspergillus fumigatus"[Organism] OR Aspergillus fumigatus[All Fields]) AND ("xylose"[MeSH Terms] OR xylose[All Fields]) AND "Expression profiling by high throughput sequencing"[Filter]''',
         '''("aspergillus fumigatus"[MeSH Terms] OR "Aspergillus fumigatus"[Organism] OR Aspergillus fumigatus[All Fields]) AND ("carbon"[MeSH Terms] OR carbon[All Fields]) AND "Expression profiling by high throughput sequencing"[Filter]''',
         '''("Chaetomium thermophilum"[All Fields] OR "Chaetomium thermophilum"[Organism] OR Chaetomium thermophilum[All Fields]) AND ("glucose"[MeSH Terms] OR glucose[All Fields]) AND "Expression profiling by high throughput sequencing"[Filter]'''
        ]

metadata = []
for query in queries:
    handle = en.esearch(db="gds",term=query,retmax=10000)
    result = en.read(handle)
    handle.close()

    print(query)
    print (len(result['IdList'])," hits \n")

    for ID in result['IdList']:
        handle2 = en.esummary(db="gds",id=ID,rettype='xml',retmode='text')
        details = en.read(handle2)[0]
        handle2.close()

        row = dict(organism = re.findall(r"\"\w+\s*\w*\"",query)[0].replace("\"",""), source = re.findall(r"\"\w+\s*\w*\"",query)[1].replace("\"",""),
        title = details['title'], Accession = details['Accession'], summary = details['summary'])
        metadata.append(row)

with open("metadata.json","w") as file:
    json.dump(list(metadata),file)

print ("\n\nSuccessfully completed !!")
from datetime import datetime
print (datetime.now().strftime('%d-%m-%Y %H:%M:%S'))
