import os
from after_blast.concatenate_exons import concatenate_exons


def concatenate_gene_results(blast_path, save_path):

    save_dir = save_path + "\\" + "organized_blast_results2"
    os.mkdir(save_dir)
    for gene in os.listdir(blast_path):
        gene_file = open(save_dir + "\\" + gene + ".fas", "a")
        gene_path = blast_path + "\\" + gene
        for taxa in os.listdir(gene_path):
            taxa_path = gene_path + "\\" + taxa
            for species in os.listdir(taxa_path):
                species_path = taxa_path + "\\" + species
                for file in os.listdir(species_path):
                    to_write = concatenate_exons(species_path + "\\" + file)
                    gene_file.write(to_write)

if __name__ == "__main__":

    blast_path = r"C:\Users\tonyx\Downloads\api_pull_complete_results7"
    save_path = r"C:\Users\tonyx\Downloads"
    concatenate_gene_results(blast_path, save_path)
