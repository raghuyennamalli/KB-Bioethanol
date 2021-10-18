#Database creation and data insertion in MongoDB 
import pymongo
import json

client = pymongo.MongoClient('mongodb://localhost:27017/')
db = client['test1']
tab = db['meta_data']


#To read JSON files
def write_metadata(filename):
    with open(filename) as file:
        metadata = json.load(file)

    tab.insert_many(metadata)
