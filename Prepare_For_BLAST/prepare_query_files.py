import os
from typing import Dict, List, Tuple

from Prepare_For_BLAST.get_longest_transcript import get_longest_transcript
from Basic_Tools.numeric_user_input import numeric_user_input
from Prepare_For_BLAST.taxonomy_browser import get_taxonomy_lineage
from Basic_Tools.lists_and_files import file_to_list, list_to_string
from Basic_Tools.basic_dictionaries import dict_get_values


def select_fill_in_manual(species_name: str, gene_name: str, gene_path: str,
                          subject_taxa: str) -> Dict or None:
    """
    Called when building a query sequence file for taxa, and the reference
    species being collected does not have a particular gene.

    :param gene_name: name of the missing gene
    :param gene_path: the filepath to the folder containing sequences for that
    gene
    :param species_name: name of the species missing the gene
    :param subject_taxa: the taxon that the query sequences currently being built will
    be BLASTed against
    :return: a dictionary storing the taxon and species of the user-selected
    species, or None if no species is selected
    """

    print("Gene: " + gene_name + " not found for " + species_name +
          ", which will be used to BLAST against " + subject_taxa)

    # given that the user has requested to fill in missing sequences manually,
    # asks the user whether the would like to fill in something for this
    # particular instance
    choose_alt_statement = "Would you like to use an alternative species for " + gene_name + "? Press 1 for yes, and 2 for no."
    choose_alt = numeric_user_input(1, 2, choose_alt_statement)

    if choose_alt == 1:

        taxa_present = os.listdir(gene_path)

        # displays subject_taxa and asks for user input
        print("Available subject_taxa:")
        for a in range(len(taxa_present)):
            print("(" + str(a) + ")" + " " + taxa_present[a])
        taxa_choice = numeric_user_input(0, len(taxa_present), "Enter a number to choose")

        # path to check out species
        taxa_path = os.path.join(gene_path, taxa_present[taxa_choice])

        # lists species that have a reference sequence for the gene
        species_present = os.listdir(taxa_path)

        # displays species and asks for user input
        print("Available reference species for this taxon:")
        for a in range(len(species_present)):
            print("(" + str(a) + ")" + " " + species_present[a])
        print("(" + str(len(species_present)) + ") None of the above, and move on")
        species_choice = numeric_user_input(0, len(species_present), "Enter a number to choose")

        # grabs the exons of the gene from the specified species
        if species_choice != len(species_present):
            print("obtaining from " + species_present[species_choice] + "...")
            return {"taxon": taxa_present[taxa_choice],
                    "species": species_present[species_choice]}

    else:
        print("moving on...")

    return None


def select_fill_in_auto(species_name, available_species: List[Dict],
                        lineage_dict: Dict) -> Dict:
    """
    Called when building a query sequence file for species_name, but is missing
    sequences for some genes. Automatically pulls from available species the
    species with the closest lineage (as stored in NCBI).

    :param species_name: name of species to build a query sequence file for
    :param available_species: a list of dictionaries, each storing a species
    that has the missing gene and the taxa the species is a part of, for species
    that have the missing gene.
    :param lineage_dict: a dictionary where the keys are species and the values
    are lists corresponding to the species' lineages.
    :return: a dictionary storing the species with the closest lineage to
    species_name
    """
    # Stores the lineage of species_name in hashmap format for quick retrieval.
    # Each key represents a taxon in the lineage, and its value is the
    # corresponding depth (0 for species level)
    host_lineage_hash = {}

    host_depth = 0
    for taxon in lineage_dict[species_name]:
        host_lineage_hash[taxon] = host_depth
        host_depth += 1

    # trace up the lineages of available species to determine how similar they
    # are to species_name. store the closest one to return later
    closest_depth = float('inf')
    closest_species = ""
    for organism in available_species:

        # trace up the species lineage, stopping when encountering the first
        # taxa shared with the host
        avail_depth = 0
        curr_taxon = lineage_dict[organism["species"]][avail_depth]
        while curr_taxon not in host_lineage_hash:
            avail_depth += 1
            curr_taxon = lineage_dict[organism["species"]][avail_depth]

        # if the depth is smaller than the closest depth so far, this species
        # is more simialr
        if host_lineage_hash[curr_taxon] < closest_depth:
            closest_depth = host_lineage_hash[curr_taxon]
            closest_species = organism

    print("Species missing: " + species_name + " closest alternative: " + closest_species["species"])

    return closest_species


def concatenate_sequences_one_query(autofill, species_name: str, taxa: str,
                                    reference_seq_path: str, save_path: str,
                                    lineage_dict) -> bool:
    """
    Assumes that exons for reference sequences have been pulled out and are in
    the folder format: general folder -> gene -> taxon -> species -> transcript.
    Concatenate all results for one species into a single file, and name the
    file after the taxon that the query will be blasted against. If no reference
    sequence of a species exists for a given gene, asks the user to input an
    alternative species to draw the sequence from (or finds one automatically).

    :param autofill: '1' for automatically pulling from the closest available
    species when a sequence is missing, '0' for manual user input.
    :param species_name: the species to concatenate results for
    :param reference_seq_path: the general folder containing reference sequences
    :param save_path: str
    :param taxa: the subject_taxa that the query will be blasted against
    :return: True iff the species has reference sequences for every
    specified gene. Also, a file will be created in the specified save_path
    """

    # query file to contain reference sequences to blast against a taxon
    file = open(os.path.join(save_path, taxa + ".fas"), "a")

    complete_species = True

    # follows through the folder structure
    for gene_folder in os.listdir(reference_seq_path):
        gene_path = os.path.join(reference_seq_path, gene_folder)

        # just in case no reference sequence exists for the species for a given
        # gene
        species_found = False
        available_species = []

        for taxa_folder in os.listdir(gene_path):
            taxa_path = os.path.join(gene_path, taxa_folder)

            for species_folder in os.listdir(taxa_path):
                species_path = os.path.join(taxa_path, species_folder)

                available_species.append({"taxon": taxa_folder,
                                          "species": species_folder})
                # if the species is the one we are looking for
                if species_folder.upper() == species_name.upper():

                    if len(os.listdir(species_path)) != 0:
                    # get the "best" transcript for the species, then append it
                    # to the query file we are building up
                        transcript_file = get_longest_transcript(species_path)
                        transcript_path = os.path.join(species_path, transcript_file)
                        file.write(open(transcript_path, "r").read() + "\n")
                        species_found = True

        if not species_found:
            complete_species = False
            # determine the alternative species to pull from
            if autofill == 1:
                # automatic
                to_select = select_fill_in_auto(species_name, available_species, lineage_dict)
            else:
                # manual
                to_select = select_fill_in_manual(species_name, gene_folder, gene_path, taxa)

            species_path = os.path.join(gene_path, to_select["taxon"], to_select["species"])
            transcript_file = get_longest_transcript(species_path)

            if transcript_file != "":
                transcript_path = os.path.join(species_path, transcript_file)
                file.write(open(transcript_path, "r").read() + "\n")
            else:
                print("Something went wrong")

    return complete_species


def get_reference_species(ref_seq_path: str) -> Dict:
    """
    Given a ref_seq_path pointing to a folder of reference sequences, identifies
    the reference species present and organizes them by taxa

    :param ref_seq_path: a string path pointing to a folder of reference
    sequences generated by the NCBI_Exon_Puller package.
    :return: a dictionary where keys are taxa and values are reference species
    belonging to a taxon.
    """
    # species encountered so far. hashmap data structure for quick verification
    # that a species has already been encountered.
    encountered_species = {}
    # output dictionary
    taxon_to_species = {}
    # iterate through folder structure
    # gene level
    for gene_folder in os.listdir(ref_seq_path):
        gene_path =os.path.join(ref_seq_path, gene_folder)
        # taxon level
        for taxon_folder in os.listdir(gene_path):
            taxon_path = os.path.join(gene_path, taxon_folder)
            # species level
            for species_folder in os.listdir(taxon_path):
                # if we haven't already added the species to the output
                if species_folder not in encountered_species:
                    encountered_species[species_folder] = True
                    # if the taxon of the species is already in the output,
                    # just append
                    if taxon_folder in taxon_to_species:
                        taxon_to_species[taxon_folder].append(species_folder)
                    else:
                        # otherwise, create a new key
                        taxon_to_species[taxon_folder] = [species_folder]

    return taxon_to_species


def auto_assign_taxa_to_ref(lineage_dct: Dict, ref_species: List[str],
                            overhead_taxon: str) -> List:
    """
    Given an overarching taxon (ex. cetacea, elasmobranchii), delegate
    sub-branches of the taxon for reference species to blast against. This is to
    ensure that query sequences are as similar to the subject as possible.

    :param lineage_dct: a dictionary where keys are species and values are lists
    corresponding to a species' lineage.
    :param ref_species: a list of reference species to delegate sub-taxa to
    :param overhead_taxon: the overarching taxon
    :return: a list of tuples each containing a reference species and the
    corresponding sub-taxon in the overhead_taxon it was assigned to. Order
    matters, in that species-taxon combos at the top of the list will be blasted
    first, and NOT be reblasted later.
    """
    lineage_assignments = {}
    result = []
    for ref in ref_species:
        ref_lineage = lineage_dct[ref]
        depth = 0
        # traces up the lineage and assigns them to itself until bumping into
        # a subject_taxa that has already been assigned
        while ref_lineage[depth] != overhead_taxon and \
                ref_lineage[depth] not in lineage_assignments:
            lineage_assignments[ref_lineage[depth]] = ref
            depth += 1

        if ref_lineage[depth] == overhead_taxon and \
                overhead_taxon not in lineage_assignments:
            # if the overhead taxon is reached for the first time
            lineage_assignments[overhead_taxon] = ref
            result.insert(0, [ref, overhead_taxon])
        else:
            # this happens after bumping into an already assigned subject_taxa
            # therefore, its subject_taxa is the one just below it
            result.insert(0, [ref, ref_lineage[depth-1]])
            # insert at the front

    return result


def get_lineage_dict_and_taxid_codes(taxa_to_species_dict: Dict) -> [Dict, Dict]:
    """
    Given the reference species and associated taxa pulled out via the
    NCBI_Exon_Puller package, find lineages/taxids for each of them.

    :param taxa_to_species_dict: a dictionary, where keys are taxa and values
    are reference species within that taxa
    :return:
    """
    # get a list of all taxa present among reference species
    taxa_list = list(taxa_to_species_dict.keys())
    # get a list of all reference_species present
    species_list = dict_get_values(taxa_to_species_dict)

    # a list of taxa/species to lookup lineages for
    # ideally, only call get_taxonomy_lineage once, because it takes time to
    # contact the website
    species_string = list_to_string(taxa_list + species_list, "\n")
    print("contacting the NCBI taxonomy API...")
    lineage_dict = get_taxonomy_lineage(species_string)

    # a dictionary where keys are taxa and values are taxids corresponding to
    # the taxa
    taxid_codes = {}
    lineage_keys = list(lineage_dict.keys())
    # remove taxa from the results, we want lineage_dict to just have species
    for key in lineage_keys:
        if key in taxa_list: # slightly inefficient, O(T), but T is pretty short
            taxid_codes[key] = lineage_dict[key][0]
            lineage_dict.pop(key, None)

    return lineage_dict, taxid_codes


def get_assignments(auto: int, lineage_dict: Dict[str, List[str]],
                    taxid_codes: Dict[str, str],
                    taxa_to_species_dict: Dict[str, List[str]]) -> List:
    """
    Get the sub-branch assignments for each reference species.

    :param auto: '1' for automatic determination of assignments, '0' for
    manually entering assignments.
    :param lineage_dict: a dictionary where keys are species names and values
    are lists of taxa corresponding to the lineage of a species
    :param taxid_codes: a dictionary where keys are taxa and values are taxids
    :param taxa_to_species_dict: a dictionary where keys are taxa and values are
    reference species in the taxa
    :return: a list of reference -> taxa assignments
    """
    # automatic assignment of reference species to sub-branches of overhead taxa
    if auto == 1:
        assignments = []
        # for each overhead taxon, get automatic assignments. concatenate them
        # at the end as a list of entries to BLAST.
        for taxon in taxid_codes.keys():
            taxid = taxid_codes[taxon]
            assignments += auto_assign_taxa_to_ref(lineage_dict, taxa_to_species_dict[taxon], taxid)
        return assignments

    # user themselves inputs a file, where it's reference sequence + taxon for
    # name
    else:
        # ex. r"C:\Users\tonyx\Downloads\concatenate_species.txt"
        assignments_file = input("Please enter a valid directory for files of "
                             "assignments")
        assignments = file_to_list(assignments_file)
        return assignments


def prepare_query_files(auto_assign, auto_fill_in, ref_seq_path, save_path) \
        -> Tuple[List[str], List[str]]:
    """
    Prepare query files based on reference sequences pulled using the
    NCBI_Exon_Puller module.

    :param auto_assign: '1' for automatic delegation of sub-branches of an
    overhead taxon to reference species in the taxon. '0' for manual assignment.
    :param auto_fill_in: '1' for automatic pulling from the closest available
    species when a species is missing sequences for a gene. '0' for automatic
    filling.
    :param ref_seq_path: The directory containing reference sequences from the
    NCBI_Exon_Puller in the original folder hierarchy.
    :param save_path: A directory to save query files to.
    """
    # get a dictionary where the keys are taxa and values are reference species
    # corresponding to the taxa
    taxa_to_species_dict = get_reference_species(ref_seq_path)

    # get the lineage dictionary and taxids for the reference species
    lineage_dct, taxid_codes = get_lineage_dict_and_taxid_codes(taxa_to_species_dict)

    # get assignments for the reference species. this is based on the lineage
    # dictionary if automatic
    assignments = get_assignments(auto_assign, lineage_dct, taxid_codes,
                                  taxa_to_species_dict)

    # for each reference species -> sub-taxon assignment, make a query file for
    # it to use in BLAST. Sometimes, the reference species are missing sequences
    # for specific genes, and automatic/manual pulling from other available
    # species is offered.
    ordered_taxa_to_be_blasted = []
    complete_species = []
    for assignment in assignments:
        if concatenate_sequences_one_query(auto_fill_in, assignment[0],
                                           assignment[1], ref_seq_path,
                                           save_path, lineage_dct):
            complete_species.append(assignment[0])
        ordered_taxa_to_be_blasted.append(assignment[1])

    return ordered_taxa_to_be_blasted, complete_species


if __name__ == "__main__":

    path = r"C:\Users\tonyx\Downloads\api_pull_complete_results8"
    save_path = r"C:\Users\tonyx\Downloads\test_taxa6"

    #if not os.path.isdir(save_path):
    #    os.mkdir(save_path)
    #print(prepare_query_files(1, 1, path, save_path))
    print(get_reference_species(r'C:\Users\tonyx\Downloads\NCBI_exon_pull_results (2)'))

