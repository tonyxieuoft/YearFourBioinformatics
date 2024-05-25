from NCBI_Exon_Puller.api_based_strategy import file_to_list


def concatenate_exons(fasta_path):

    fasta_list = file_to_list(fasta_path)
    sections = fasta_list[0].split(" ")

    section_no = 0
    while len(sections[section_no]) < 6 or sections[section_no][:6] != "genome":
        section_no += 1
    exon_section = section_no + 2

    global_beginning = sections[exon_section].split("-")[0]

    sequence = ""
    last_range = None
    line_no = 0
    while line_no < len(fasta_list):
        coverage = fasta_list[line_no].split(" ")[exon_section].split("-")

        line_no += 1

        if last_range is None:
            if global_beginning != 1:
                sequence += (int(global_beginning) - 1)*"-"
            sequence += fasta_list[line_no]
        else:
            in_between = int(coverage[0]) - int(last_range[1]) - 1
            filler = "-"*in_between
            sequence += filler + fasta_list[line_no]

        last_range = coverage
        line_no += 1

    fasta_heading = ""
    for i in range(exon_section):
        if i == 0:
            fasta_heading = sections[i]
        else:
            fasta_heading += " " + sections[i]

    fasta_heading += " " + global_beginning + "-" + last_range[1] + "\n"

    return fasta_heading + sequence + "\n"


if __name__ == "__main__":
    fasta_path = r"C:\Users\tonyx\Downloads\blast_test_results1234567891011121314151617\rdh8\lamniformes\Carcharodon carcharias\Reference_XM_041177394.1.fas"
    concatenate_exons(fasta_path)
