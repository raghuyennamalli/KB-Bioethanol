import urllib.request, time, configparser, os
from bs4 import BeautifulSoup
from tqdm import tqdm
# reading config file
config = configparser.ConfigParser()
config.read('../config.ini')
# 
def get_time(function_name): 
    def inner1(*args, **kwargs): 
        begin = time.time() 
        function_output = function_name(*args, **kwargs) 
        end = time.time() 
        print("Total time taken in : ", function_name.__name__, end - begin,flush=True) 
        return function_output
    return inner1
#
# --------------------------------------------------------------------
# Function to fetch data from NCBI for the provided query and database
str_ncbi_baseurl = config["NCBI"]["base_url"]
str_ncbi_apikey = config["NCBI"]["api_key"]
def get_data_from_ncbi(str_term=None,str_db_name=None,int_retmax=9900,str_output_path=None):
    """
    Description: To search in NCBI portal using esearch API for provided term in the provided database
    Input:
        str_term: query to be searched in NCBI portal. For two terms follow below pattern:
                  EGFR AND Breast Cancer
        str_db_nam: name of the NCBI database where query has to be searched
        int_retmax: total number of documents to be retrieved in one batch. Should be less than 10,000.
        str_output_path: path to a directory where output will be stored in xml file format
    Output:
    Function return True if no issues are encountered else False.
    In case of True, XML files containing data for the searched query will be stored at provided output path.
    """
    if str_term is None or str_db_name is None or str_output_path is None:
        print("Input parameters were not passed as expected by the function. Function will exit !")
        return False
    # format term query: replace spaces with "+"; convert string to lowercase
    str_term = str_term.replace(" ","+").lower().replace("+and+","+AND+").replace("+or+","+OR+").replace("+not+","+NOT+")
    # form esearch url using base url and query
    str_esearch_url = str_ncbi_baseurl+"esearch.fcgi?db={}&term={}&retmax={}&usehistory=y&api_key={}".format(str_db_name,str_term,int_retmax,str_ncbi_apikey)
    # get response of esearch url
    try:
        with urllib.request.urlopen(str_esearch_url) as response:
            bs_esearch_response = response.read()
    except Exception as e:
        print("Error encountered while searching in NCBI: {}".format(str(e)))
        return False
    # parse response using beautifulsoup and get values
    bs_esearch_response = BeautifulSoup(bs_esearch_response,"lxml")
    str_webenv = bs_esearch_response.find("webenv").get_text()
    str_querykey = bs_esearch_response.find("querykey").get_text()
    int_count = int(bs_esearch_response.find("count").get_text())
    if int_count > 0:
        os.makedirs(str_output_path, exist_ok=True)
    int_retstart = int(bs_esearch_response.find("retstart").get_text())
    print("Total number of hits for the term {} is: {}".format(str_term,int_count))
    # iterate till int_count and fetch data in batches using retmax
    for counter,i in (enumerate(range(int_retstart,int_count,int_retmax))):
        str_fetch_url = str_ncbi_baseurl+"efetch.fcgi?db={}&WebEnv={}&query_key={}&retstart={}&retmax={}&retmode={}&api_key={}".format(str_db_name,str_webenv,str_querykey,i,int_retmax,"xml",str_ncbi_apikey)
        try:
            # get data
            with urllib.request.urlopen(str_fetch_url) as response:
                bs_efetch_response = response.read()
            bs_efetch_response = BeautifulSoup(bs_efetch_response,"lxml")
            # write output data to xml file batch-wise
            str_filepath = os.path.join(str_output_path,"batch_{}.xml".format(counter))
            with open(str_filepath,"w") as f:
                f.write(bs_efetch_response.prettify())
        except Exception as e:
            print("Error encountered while fetching data from NCBI for {} to {}: {}".format(i,i+int_retmax,str(e)))
    return True
#
