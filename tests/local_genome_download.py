import os
from Basic_Tools.basic_dictionaries import json_to_dict, print_dict
import subprocess

if __name__ == "__main__":
    # tempjson.txt location will instead be the folder the user specifies
    temp_file_name = r'tempjson.txt'
    os.system("datasets summary genome taxon 'orcinus orca' > " + temp_file_name)
    dct = json_to_dict(temp_file_name)

    accession = ""
    species = ""
    for genome_record in dct["reports"]:
        if "refseq_category" in genome_record["assembly_info"]:
            print(genome_record["organism"]["organism_name"])
            species = genome_record["organism"]["organism_name"]

            print(genome_record["accession"])
            accession = genome_record["accession"]

    os.system(r"datasets download genome accession " + accession)
    os.system("unzip ncbi_dataset.zip '" + species + "'")
    os.system("rm ncbi_dataset.zip")
    os.system("cd '" + species + "'/data/" + accession)
    genome_file = str(subprocess.check_output(["ls"])).strip()
    print("This is the genome file: " + genome_file)


