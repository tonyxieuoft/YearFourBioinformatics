from Bio.Blast import NCBIWWW, NCBIXML
from typing import List
import os


def local_genome_blaster(save_path: str, queries_path: str,
                         taxa_blast_order: List[str],
                         complete_reference_species: List[str],
                         expect_value):

    species_so_far = {}
    for species in complete_reference_species:
        species_so_far[species] = True

    for taxa in taxa_blast_order:

        reference_filepath = os.path.join(queries_path, taxa + ".fas")



if __name__ == "__main__":
    NCBIWWW.email = "xiaohan.xie@mail.utoronto.ca"
    query_sequences = open(r'C:\Users\tonyx\Downloads\query_api_test.txt').read()
    print(query_sequences)
    test(query_sequences)
