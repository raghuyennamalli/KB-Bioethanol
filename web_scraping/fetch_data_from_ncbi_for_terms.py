import utils_webscraping as utils
import os
if __name__ == '__main__':
    list_term = ["Casey S Greene","Jayashree Ramana"]
    for term in list_term:
        str_dir_name = os.path.join("../data/",term.replace(" ","_"))
        #os.makedirs(str_dir_name, exist_ok=True)
        bool_response = utils.get_data_from_ncbi(str_term=term,
                                                 str_db_name="pubmed",
                                                 str_output_path=str_dir_name)
        print("Response for term {} is: {}".format(term,bool_response))