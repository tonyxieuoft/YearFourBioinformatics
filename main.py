import os
from Bio import Entrez
from Bio.Blast import NCBIWWW

from After_BLAST.concatenate_gene_results import concatenate_gene_results
from Basic_Tools.lists_and_files import file_to_list, make_unique_directory, \
    unique_filepath
from Basic_Tools.numeric_user_input import numeric_user_input

# Press the green button in the gutter to run the script.
from NCBI_Exon_Puller.handle_ncbi_exon_puller import handle_ncbi_exon_puller
from NCBI_Genome_Blaster.driver_genome_blaster import driver_genome_blaster
from Prepare_For_BLAST.prepare_query_files import prepare_query_files
from User_Interaction.user_exon_pulling import enter_gene_filepath, \
    get_generic_directory, get_generic_filepath

if __name__ == '__main__':

    print("Welcome to the pipeline...")

    email_confirm = 2
    email = ""
    while email_confirm == 2:

        print("NCBI Entrez and BLAST require an email to contact in case any "
              "issues arise. Please enter your email. ")
        email = input()

        print("Inputted email: \"" + email + "\"")
        email_confirm_prompt = "Enter 1 to confirm, or 2 to re-enter your email"
        email_confirm = numeric_user_input(1, 2, email_confirm_prompt)

    Entrez.email = email

    valid_directory = False
    download_dir = ""
    while not valid_directory:

        print("Please specify a valid directory to which files created by the "
              "pipeline can be downloaded to.")
        download_dir = input()

        valid_directory = os.path.isdir(download_dir)
        if not valid_directory:
            print("Invalid directory: " + download_dir)

    print("Great! Now, we will begin pulling exons from NCBI. If you would "
          "like to skip this step and start blasting immediately, enter 1. "
          "Otherwise, to begin pulling exons from NCBI, enter 2.")

    skip_pull_step = numeric_user_input(1, 2, "")
    if skip_pull_step == 1:
        print("skipping exon pulling...")
        print("Enter a path to a directory containing reference sequences "
              "previously pulled out using the NCBI exon puller.")
        exon_pull_path = get_generic_directory()
    else:
        gene_query_filepath = enter_gene_filepath()

        print("Now, enter a valid file path containing taxa of interest.")
        taxa_filepath = get_generic_filepath()

        print("Ok, we are all set to pull reference sequences from "
              "the specified taxa for the specified genes.")

        exon_pull_path = make_unique_directory(download_dir,
                                               "NCBI_exon_pull_results")
        print("The results will be at the following path: " +
              exon_pull_path)

        print("Starting NCBI Exon Puller...")
        handle_ncbi_exon_puller(exon_pull_path, gene_query_filepath,
                                taxa_filepath)
        print("Finished pulling exons!")

        termination_message = "Enter 1 to terminate the program here, or enter " \
                              "2 to prepare for blast"
        if numeric_user_input(1, 2, termination_message) == 1:
            print("exiting...")
            exit(0)
    print("============================================")
    print("Now, we'll prepare query files for blasting.")
    print("")
    assign_message = "If you would to assign reference species to subject taxa " \
                     "manually, enter 2. Otherwise, for automatic assignment, " \
                     "enter 1."
    auto_assign = numeric_user_input(1, 2, assign_message)

    print("")
    fill_in_message = "Sometimes, a reference species will not have sequences " \
                      "for certain genes. In cases where this occurs, the " \
                      "sequences from the closest available species will be " \
                      "used instead. Enter 1 for this fill-in process to be " \
                      "automatic. Enter 2 if you would like to manually select " \
                      "the closest species during the fill-in process."
    auto_fill = numeric_user_input(1, 2, fill_in_message)

    queries_path = make_unique_directory(download_dir, "query_files")
    print("Great. Query files will be at the path: " + queries_path)
    print("=============================================")
    print("Making query files...")

    taxa_blast_order, reference_species = \
        prepare_query_files(auto_assign, auto_fill, exon_pull_path, queries_path)

    print(taxa_blast_order)

    blast_results_path = make_unique_directory(download_dir, "blast_results")
    print("Blast results will be at the path: " + blast_results_path)
    print("Enter any key to begin the blasting process")
    input()

    NCBIWWW.email = email
    driver_genome_blaster(blast_results_path, queries_path, taxa_blast_order,
                          reference_species)

    print("BLAST complete.")
    print("=============================================")
    alignments_path = make_unique_directory(download_dir, "alignments")
    print("Organizing BLAST results into alignments, will be at the path: " +
          alignments_path)
    print("...")

    concatenate_gene_results(exon_pull_path, blast_results_path, alignments_path)

    print("=============================================")
    print("Program finished.")










