import os

from ncbi_exon_puller.ncbi_exon_puller import ncbi_get_gene_sequence
from Bio import Entrez
from basic_tools.xml_extraction import xml_to_dictionary

SEQUENCE_INDICES_FROM_MRNA_TAG = 2


def get_directory(parent_directory: str, directory_name: str) -> str:

    directory_path = parent_directory + "\\" + directory_name
    if not os.path.isdir(directory_path):
        os.mkdir(directory_path)

    return directory_path


def parse_blast_xml(file: str, save_dir: str, taxon_name: str, curr_species):

    Entrez.email = "xiaohan.xie@mail.utoronto.ca"

    results_dict = xml_to_dictionary(file)
    iterations = results_dict['BlastOutput']['BlastOutput_iterations']
    for iteration in iterations:
        print(iteration['Iteration_iter-num'])
        query_title = iteration['Iteration_query-def'].split(" ")
        query_seq_length = int(iteration['Iteration_query-len'])
        gene_name = query_title[0]

        mrna_section_no = 1
        while len(query_title[mrna_section_no]) < 4 or \
                query_title[mrna_section_no][:4] != "mRNA":
            mrna_section_no += 1

        ref_transcript_var = query_title[mrna_section_no].split(":")[1]
        ref_sequence_range = query_title[mrna_section_no +
                                         SEQUENCE_INDICES_FROM_MRNA_TAG]

        if isinstance(iteration['Iteration_hits'], list) or \
                isinstance(iteration['Iteration_hits'], dict): # if not empty

            if isinstance(iteration['Iteration_hits'], list):
                top_hit = iteration['Iteration_hits'][0]
            else:
                top_hit = iteration['Iteration_hits']['Hit']

            if isinstance(top_hit["Hit_hsps"], list):
                seq_data = top_hit["Hit_hsps"][0]
            else:
                seq_data = top_hit["Hit_hsps"]["Hsp"]

            result_sequence = seq_data["Hsp_hseq"]
            query_bound1 = int(seq_data["Hsp_query-from"])
            query_bound2 = int(seq_data["Hsp_query-to"])


            missing_left = query_bound1 - 1
            missing_right = query_seq_length - query_bound2
            # given this, can easily fill up with "N"s or "-"s

            accession = top_hit['Hit_accession']

            hit_max = int(top_hit['Hit_len'])

            if missing_left > 0 or missing_right > 0:

                hit_bound1 = int(seq_data['Hsp_hit-from'])
                hit_bound2 = int(seq_data['Hsp_hit-to'])

                if hit_bound1 < hit_bound2:
                    strand = "1"
                else:
                    strand = "2"

                if missing_left > 0:

                    to_salvage = missing_left

                    if strand == "1":
                        lower_bound = hit_bound1 - missing_left
                        upper_bound = hit_bound1 - 1
                        if upper_bound < 1:
                            to_salvage = 0
                        elif lower_bound < 1:
                            lower_bound = 1
                            to_salvage = upper_bound - lower_bound + 1
                    else:
                        lower_bound = hit_bound1 + missing_left
                        upper_bound = hit_bound1 + 1
                        if upper_bound > hit_max:
                            to_salvage = 0
                        elif lower_bound > hit_max:
                            lower_bound = hit_max
                            to_salvage = lower_bound - upper_bound + 1

                    string = ""
                    if to_salvage > 0:
                        arr = ncbi_get_gene_sequence(accession, lower_bound,
                                                 upper_bound, strand)
                        for ch in arr:
                            string += ch

                    string = "-"*(missing_left - to_salvage) + string

                    print("missing left recovered: " + string + " iteration: " + iteration['Iteration_iter-num'])

                    result_sequence = string + result_sequence

                if missing_right > 0:

                    to_salvage = missing_right

                    if strand == "1":
                        lower_bound = hit_bound2 + 1
                        upper_bound = hit_bound2 + missing_right

                        if lower_bound > hit_max:
                            to_salvage = 0
                        elif upper_bound > hit_max:
                            upper_bound = hit_max
                            to_salvage = upper_bound - lower_bound + 1

                    else:
                        lower_bound = hit_bound2 - 1
                        upper_bound = hit_bound2 - missing_right

                        if lower_bound < 1:
                            to_salvage = 0
                        elif upper_bound < 1:
                            upper_bound = 1
                            to_salvage = lower_bound - upper_bound + 1

                    string = ""
                    if to_salvage > 0:
                        arr = ncbi_get_gene_sequence(accession, lower_bound, upper_bound, strand)

                        for ch in arr:
                            string += ch

                    string = string + "-"*(missing_right - to_salvage)

                    print("missing right recovered: " + string + " iteration: " + iteration['Iteration_iter-num'])

                    result_sequence = result_sequence + string

            gene_folder = get_directory(save_dir, gene_name)
            taxa_folder = get_directory(gene_folder, taxon_name)
            species_folder = get_directory(taxa_folder, curr_species)

            transcript_file = open(species_folder + "\\" + "Reference_" +
                                   ref_transcript_var + ".fas", "a")

            fasta_heading = ">" + gene_name + " " + curr_species + " reference_mrna:" + ref_transcript_var + " genome:" + accession + " " + ref_sequence_range
            transcript_file.write(fasta_heading + "\n")
            transcript_file.write(result_sequence + "\n")
            transcript_file.close()


if __name__ == "__main__":
    Entrez.email = "xiaohan.xie@mail.utoronto.ca"
    parse_blast_xml(r"C:\Users\tonyx\Downloads\51JMZR2N016-Alignment.xml", r"C:\Users\tonyx\Downloads\xml_readtest", "altered2", "some weird one")





