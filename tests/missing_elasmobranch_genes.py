import os


def missing_elasmobranch_genes(general_folder: str):

    for gene in os.listdir(general_folder):

        gene_path = os.path.join(general_folder, gene)
        elasmobranch_path = os.path.join(gene_path, "elasmobranchii")
        if len(os.listdir(elasmobranch_path)) == 0:
            print(gene)


if __name__ == "__main__":
    pass
    # missing_elasmobranch_genes(r'C:\Users\tonyx\Downloads\NCBI_exon_pull_results (10)')
