from Bio.Blast import NCBIWWW, NCBIXML

def test(query_sequences):

    entrez_query = "txid9733[orgn]"
    result_handle = NCBIWWW.qblast("blastn", "refseq_rna",
                                   entrez_query=entrez_query)
    # ref_euk_rep_genomes
    contents = open(result_handle).read()
    new_file = open(r"C:\Users\tonyx\Downloads\blast_output_test.xml", "w")
    new_file.write(contents)


if __name__ == "__main__":
    NCBIWWW.email = "xiaohan.xie@mail.utoronto.ca"
    query_sequences = open(r'C:\Users\tonyx\Downloads\query_api_test.txt').read()
    print(query_sequences)
    test(query_sequences)
