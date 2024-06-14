import os
import subprocess

from Basic_Tools.basic_dictionaries import json_to_dict


if __name__ == "__main__":

    working_path = subprocess.check_output(["pwd"], shell=True). \
        decode("utf-8").strip()
    prev_blastdb_path = subprocess.check_output(["echo", "$BLASTDB"], shell=True). \
        decode("utf-8").strip()
    os.environ["BLASTDB"] = prev_blastdb_path + working_path + "/local_blast_db"

    os.system("datasets summary genome taxon 7797 > 'temp.txt'")
    genome_summary_dict = json_to_dict("temp.txt")
    os.remove("temp.txt")

    accession = ""
    for genome_record in genome_summary_dict["reports"]:

        if "refseq_category" in genome_record["assembly_info"]:

            print(genome_record["organism"]["organism_name"])
            print(genome_record["accession"])
            accession = genome_record["accession"]

    os.system(r"datasets download genome accession " +
              accession + " --dehydrated")

    # unzip it into a generic, temporary directory called ncbi_dataset
    os.system("unzip ncbi_dataset.zip")

    # remove zip file
    os.system("rm ncbi_dataset.zip")
    os.system("rm README.md")

    os.system("datasets rehydrate --directory " + working_path)


    # get genome fasta file path nested within the temp directory
    genome_file_directory = os.path. \
        join("ncbi_dataset", "data", accession)
    genome_filename = subprocess. \
        check_output(["ls", genome_file_directory]). \
        decode("utf-8").strip()
    genome_filepath = os.path. \
        join(genome_file_directory, genome_filename)

    # create local blast database using the genome file
    os.system("mkdir local_blast_db")
    os.system("makeblastdb -dbtype nucl -in " + genome_filepath +
              " -out local_blast_db/local")

    # remove original genome file
    os.system("rm -r ncbi_dataset")

    # run blastn on the local database
    reference_filepath = "'/mnt/c/Users/tonyx/Downloads/of_interest/NCBI_exon_pull_results (11)/CNGA3/elasmobranchii/Amblyraja radiata'"
    os.system("blastn -db local -outfmt 5 -query " +
              reference_filepath + " > out.xml")

    # delete the local blast database
    os.system("rm -r local_blast_db")

