import os
from Basic_Tools.basic_dictionaries import json_to_dict, print_dict


if __name__ == "__main__":
    # tempjson.txt location will instead be the folder the user specifies
    temp_file_name = r'tempjson.txt'
    os.system("datasets summary genome taxon 'orcinus orca' > " + temp_file_name)
    dct = json_to_dict(temp_file_name)

    for genome_record in dct["reports"]:
        if "refseq_category" in genome_record["assembly_info"]:
            print(genome_record["organism"]["organism_name"])
            print(genome_record["accession"])


