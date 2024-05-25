import os

from prepare_for_blast.get_longest_transcript import get_longest_transcript
from basic_tools.numeric_user_input import numeric_user_input


def fill_in_missing_genes(gene_folder: str, gene_path: str, species_name: str,
                          taxa: str, file) -> None:
    """
    Called when building a query sequence file for taxa, and the reference
    species being collected does not have a particular gene.

    :param gene_folder: name of the missing gene
    :param gene_path:
    :param species_name: name of the species missing the gene
    :param taxa: the taxon that the query sequences currently being built will
    be BLASTed against
    :param file: contains the query
    :return: nothing. via user input, appends sequences from available species
    to the query file when the specified species isn't present
    """

    print("Gene: " + gene_folder + " not found for " + species_name +
          ", which will be used to BLAST against " + taxa)

    choose_alt_statement = "Would you like to use an alternative species for " + gene_folder + "? Press 1 for yes, and 2 for no."
    choose_alt = numeric_user_input(1, 2, choose_alt_statement)

    if choose_alt == 1:

        # lists taxa to choose from
        print("Available taxa:")
        taxa_present = os.listdir(gene_path)
        for a in range(len(taxa_present)):
            print("(" + str(a) + ")" + " " + taxa_present[a])
        taxa_choice = numeric_user_input(0, len(taxa_present),
                                         "Enter a number to choose")
        taxa_path = gene_path + "\\" + taxa_present[taxa_choice]

        # lists species that have a reference sequence for the gene
        print("Available reference species for this taxon:")
        species_present = os.listdir(taxa_path)
        # menu listing
        for a in range(len(species_present)):
            print("(" + str(a) + ")" + " " + species_present[a])
        # option for selecting none of them
        print("(" + str(len(species_present)) + ") None of the "
                                                "above, and move on")

        species_choice = numeric_user_input(0, len(species_present),
                                            "Enter a number to choose")
        # grabs the exons of the gene from the specified species
        if species_choice != len(species_present):
            print("obtaining from " + species_present[species_choice] +
                  "...")
            species_path = taxa_path + "\\" + species_present[species_choice]
            transcript_file = get_longest_transcript(species_path)
            if transcript_file != "":
                transcript_path = species_path + "\\" + transcript_file
                file.write(open(transcript_path, "r").read() + "\n")
            else:
                print("Something went wrong")

    else:
        print("moving on...")


def combine_species_sequences(species_name: str, reference_seq_path: str,
                              save_path: str, taxa: str) -> None:
    """
    Assumes that exons for reference sequences have been pulled out and are in
    the folder format: general folder -> gene -> taxon -> species -> transcript.
    Concatenate all results for one species into a single file, and name the
    file after the taxon that the query will be blasted against. If no reference
    sequence of a species exists for a given gene, asks the user to input an
    alternative species to draw the sequence from.

    :param species_name: the species to concatenate results for
    :param reference_seq_path: the general folder containing reference sequences
    :param save_path: str
    :param taxa: the taxa that the query will be blasted against
    :return: nothing. a file will be created in the specified save_path
    """

    # query file to contain reference sequences to blast against a taxon
    file = open(save_path + "\\" + taxa + ".fas", "a")

    # follows through the folder structure
    for gene_folder in os.listdir(reference_seq_path):
        gene_path = reference_seq_path + "\\" + gene_folder
        species_found = False
        for taxa_folder in os.listdir(gene_path):
            taxa_path = gene_path + "\\" + taxa_folder

            for species_folder in os.listdir(taxa_path):

                if species_folder.upper() == species_name.upper():
                    species_path = taxa_path + "\\" + species_folder

                    transcript_file = get_longest_transcript(species_path)
                    if transcript_file != "":
                        transcript_path = species_path + "\\" + transcript_file
                        file.write(open(transcript_path, "r").read() + "\n")
                        species_found = True

        if not species_found:
            fill_in_missing_genes(gene_folder, gene_path, species_name, taxa,
                                  file)


if __name__ == "__main__":

    save_path = r"C:\Users\tonyx\Downloads\test_taxa4"
    path = r"C:\Users\tonyx\Downloads\api_pull_complete_results8"
    if not os.path.isdir(save_path):
        os.mkdir(save_path)

    command_file = open(r"C:\Users\tonyx\Downloads\concatenate_species.txt","r")
    commands = command_file.readlines()
    for i in range(len(commands)):
        if commands[i].strip() != "":
            commands[i] = commands[i].strip().split(",")
            combine_species_sequences(commands[i][0], path, save_path, commands[i][1])

    combine_species_sequences("Amblyraja radiata", path, save_path, "elasmobranchii")

