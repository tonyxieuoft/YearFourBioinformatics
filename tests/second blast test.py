from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML
import os.path


def unique(path):
    final_path = path
    counter = 1
    while os.path.exists(final_path):
        final_path = path + str(counter)
        counter += 1
    return final_path


def remove_create_view(file):
    """
    Takes in either the path to an XML file or an open XML file, removes any
    instances of the phrase "CREATE_VIEW" then returns none if a path was
    inputted, or the file contents without "CREATE_VIEW" if an open XML was
    inputted

    :param file: either an open XML or the path to an XML file

    :return: None if a path was entered, or the open XML without
    "CREATE_VIEW" if an open XML was entered
    """

    try:  # stay within the try block if a path was entered
        file_string = open(file, 'r').read()
        file_string = file_string.replace('CREATE_VIEW', '')
        with open(file, 'w') as out:
            out.write(file_string)

    except TypeError:  # if a TypeError is thrown, file is an open XML

        file_string = file.read()
        file_string = file_string.replace('CREATE_VIEW', '')
        file.close()

        # write the contents of the file to a temporary file, open that and
        # delete the temporary file, then return the open file
        path = unique(r'C:\Users\tonyx\Downloads\temp')
        with open(path, 'w') as out:
            out.write(file_string)
        out.close()
        newfile = open(path, 'r')
        return newfile


NCBIWWW.email = "xiaohan.xie@mail.utoronto.ca"
query_sequences = r"C:\Users\tonyx\Downloads\blasttest.txt"
print("querying")
result_handle = NCBIWWW.qblast("blastn", "nt",
                               query_sequences,
                               entrez_query="txid9733[orgn]")
print("done")
result_handle = remove_create_view(result_handle)
blast_records = NCBIXML.parse(result_handle)


