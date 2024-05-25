import shutil
import time
import zipfile
import os
from selenium import webdriver
from selenium.webdriver import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

SEQID = 0
SOURCE = 1
REGION_TYPE = 2
SEQ_START = 3
SEQ_END = 4
SCORE = 5
STRAND = 6
PHASE = 7
ATTRIBUTES = 8


def exon_from_gene(line_features, ascending, endpoints, dna_array, directory,
                   file, line):
    sequence = ""

    if ascending:
        startpoint = int(line_features[SEQ_START]) - \
                     int(endpoints[0])
        endpoint = int(line_features[SEQ_END]) - \
                   int(endpoints[0])
    else:
        startpoint = int(endpoints[0]) - \
                     int(line_features[SEQ_END])
        endpoint = int(endpoints[0]) - \
                   int(line_features[SEQ_START])

    if startpoint >= 0 and endpoint <= len(dna_array) - 1:
        for i in range(startpoint, endpoint):
            sequence = sequence + dna_array[i]
    else:
        return False

    f = open(directory + "\\" + file, "w")
    fasta_heading = ">" + line_features[SEQID] + "_" + line_features[SEQ_START] \
                    + "-" + line_features[SEQ_END] + "_" + \
                    line_features[REGION_TYPE] +"\n"
    f.write(fasta_heading)
    f.write(sequence)
    f.close()

    return True


def extract_genome_from_file(directory):

    file = open(directory + r'\ncbi_dataset\data\gene.fna', 'r')
    line = file.readline()
    endpoints = line.split(":")[1]
    if endpoints[0] == "c":
        endpoints = endpoints[1:]

    endpoints = endpoints.split(" ")[0].split("-")

    # gets the gene
    dna_array = []
    line = file.readline()
    while line != "":
        for ch in line.strip():
            dna_array.append(ch)
        line = file.readline()

    file.close()

    return dna_array, endpoints


def genome_segment_extract(gene, save_location, directory): #outside makes the folder for the organism
    # we are provided the gene folder here, its just the spec. organism

    time.sleep(0.5)

    with zipfile.ZipFile(directory + r'\genome_segment.zip', 'r') as zip_ref:
        zip_ref.extractall(directory)

    os.remove(directory + r'\genome_segment.zip')

    dna_array, endpoints = extract_genome_from_file(directory)

    ascending = True
    if endpoints[0] > endpoints[1]:
        ascending = False

    # given that everything is sent to one directory, very consistent
    features_path = ""

    while features_path == "":
        for name in os.listdir(directory):
            if os.path.splitext(name)[1] == ".GFF3":
                features_path = directory + "\\" + name
                break

        time.sleep(0.05)

    features_file = open(features_path, "r")
    line = features_file.readline()
    while line.strip()[0] == "#":
        line = features_file.readline()

    current_seq_folder = ""
    file_counter = 0

    while line != "" and line[0] != "#":

        line_features = line.strip().split("\t")
        if line_features[REGION_TYPE] == "exon":

            attribute_features = line_features[ATTRIBUTES].split(";")

            # sometimes other genes are also in the track
            correct_gene = False
            reference = ""
            for af in attribute_features:
                key, value = af.split("=")[0], af.split("=")[1]
                if key == "gene" and value.upper() == gene.upper():
                    correct_gene = True
                elif key == "Parent":
                    reference = value

            if correct_gene:

                if reference == "":
                    reference = "unnamed" + str(file_counter)
                    file_counter += 1

                transcript_filepath = save_location + "\\" + str(reference)
                os.mkdir(transcript_filepath)

                exon_from_gene(line_features, ascending, endpoints,
                               dna_array, transcript_filepath, "exon1", line)

                line = features_file.readline()
                line_features = line.strip().split("\t")
                counter = 2
                while line != "" and line[0] != "#" and line_features[REGION_TYPE] == "exon":

                    exon_from_gene(line_features,ascending, endpoints,
                                   dna_array, transcript_filepath, "exon" +
                                   str(counter), line)

                    line = features_file.readline()
                    line_features = line.split("\t")
                    counter += 1

        line = features_file.readline()

    features_file.close()
    os.remove(features_path)


def file_to_list(filepath): # gene file is a txt file

    f = open(filepath)
    raw = f.readlines()
    output = []
    for i in range(0, len(raw)):
        if raw[i].strip() != "":
            output.append(raw[i].strip())

    return output


def click_downloads(driver, path):
    WebDriverWait(driver, 10).until(
        EC.presence_of_element_located((By.ID, "download-asm")))

    driver.find_element(By.ID, "download-asm").click()
    WebDriverWait(driver, 10).until(
        EC.element_to_be_clickable((By.ID, "datasets-download-submit")))

    WebDriverWait(driver, 10).until(
        EC.element_to_be_clickable((By.ID, "datasets-download-file-name")))
    filename_box = driver.find_element(By.ID, "datasets-download-file-name")
    # filename_box.clear()
    filename_box.click()
    filename_box.send_keys(Keys.CONTROL + "a")
    filename_box.send_keys(Keys.DELETE)
    filename_box.send_keys("genome_segment.zip")
    driver.find_element(By.ID, "datasets-download-submit").click()

    while not os.path.exists(path + "\\genome_segment.zip"):
        time.sleep(0.2)

    WebDriverWait(driver, 10).until(
        EC.element_to_be_clickable((By.ID, "button-1084")))
    driver.find_element(By.ID, "button-1084").click()

    WebDriverWait(driver, 10).until(
        EC.presence_of_element_located((By.ID, "menuitem-1056-itemEl")))
    driver.find_element(By.ID, "menuitem-1056-itemEl").click()

    WebDriverWait(driver, 10).until(
        EC.element_to_be_clickable((By.CSS_SELECTOR, ".x-toolbar.x-docked.x-toolbar-footer.x-box-layout-ct")))
    pop_up_box = driver.find_element(By.CSS_SELECTOR, ".x-toolbar.x-docked.x-toolbar-footer.x-box-layout-ct")

    buttons = pop_up_box.find_elements(By.XPATH, ".//a")
    buttons[1].click()

    downloaded = False
    iterations = 0
    while not downloaded and iterations < 200:
        textboxes = driver.find_elements(By.CSS_SELECTOR, ".x-toolbar-text.x-box-item.x-toolbar-item.x-toolbar-text-default")
        for textbox in textboxes:
            if textbox.text == "File created":
                downloaded = True

        iterations += 1
        time.sleep(0.1)

    buttons[2].click()


def get_exons():

    # all of these wil be input statements (so there are only 3)
    user_specified_directory = r'C:\Users\tonyx\Downloads'
    genes_filepath = r'C:\Users\tonyx\Downloads\genestest.txt'
    taxon_filepath = r'C:\Users\tonyx\Downloads\taxontest.txt'

    path = user_specified_directory + "\\exon_pull_complete_results"
    # make a folder beforehand, error check later to see if this folder exists
    os.mkdir(path)

    genes = file_to_list(genes_filepath) # try-except clauses here?
    taxa = file_to_list(taxon_filepath)

    # set the webdriver's download path
    chrome_options = webdriver.ChromeOptions()
    prefs = {'download.default_directory': path, 'profile.default_content_setting_values.automatic_downloads': 1}
    chrome_options.add_experimental_option('prefs', prefs)
    driver = webdriver.Chrome(chrome_options)

    # access the gene search page, which will have the refseq data
    driver.get("https://www.ncbi.nlm.nih.gov/gene/")

    for gene in genes:

        gene_folder = path + "\\" + gene
        os.mkdir(gene_folder)

        for taxon in taxa:

            taxon_folder = gene_folder + "\\" + taxon
            os.mkdir(taxon_folder)

            search_query = gene + "[gene] \"" + taxon + "\"[organism]"

            search_box = driver.find_element(By.NAME, "term")

            # clear previous search from box
            search_box.click()
            search_box.send_keys(Keys.CONTROL + "a")
            search_box.send_keys(Keys.DELETE)

            # enter new search
            search_box.send_keys(search_query)

            # click enter
            driver.find_element(By.ID, "search").click()

            # wait for page to load, then determine whether a single result,
            # a grid of results, or no results (error) showed up.
            page_type = ""
            iterations = 0
            while page_type == "" and iterations < 100:
                if len(driver.find_elements(By.CSS_SELECTOR, ".jig-ncbigrid.gene-tabular-rprt.ui-ncbigrid")) > 0:
                    page_type = "grid"
                elif len(driver.find_elements(By.ID, "gene-name")) > 0:
                    page_type = "single"
                elif len(driver.find_elements(By.CSS_SELECTOR, ".warn.icon")) > 0:
                    page_type = "error"

                time.sleep(0.1)
                iterations += 1

            print(page_type)

            if page_type == "grid":
                table = driver.find_element(By.CSS_SELECTOR,
                                            ".jig-ncbigrid.gene-tabular-rprt.ui-ncbigrid")
                elements = table.find_elements(By.CLASS_NAME, "rprt")
                for i in range(len(elements)):
                    gene_link = elements[i].find_element(By.XPATH, ".//a")
                    row_species = elements[i].find_element(By.XPATH, ".//em")

                    species_folder = taxon_folder + "\\" + row_species.text
                    os.mkdir(species_folder)

                    gene_link.click()
                    click_downloads(driver, path)
                    genome_segment_extract(gene, species_folder, path)
                    driver.execute_script("window.history.go(-1)")

                    #reset the page
                    WebDriverWait(driver, 10).until(
                        EC.presence_of_element_located(
                            (By.CSS_SELECTOR, ".jig-ncbigrid.gene-tabular-rprt.ui-ncbigrid")
                        )
                    )
                    table = driver.find_element(By.CSS_SELECTOR,
                                                ".jig-ncbigrid.gene-tabular-rprt.ui-ncbigrid")
                    elements = table.find_elements(By.CLASS_NAME, "rprt")

            elif page_type == "single":

                species = driver.find_element(By.CLASS_NAME, "tax")
                species_folder = taxon_folder + "\\" + species.text
                os.mkdir(species_folder)

                click_downloads(driver, path)
                genome_segment_extract(gene, species_folder, path)

    shutil.rmtree(path + "\\ncbi_dataset")
    os.remove(path + "\\README.MD")

"""

        


        time.sleep(6)
"""

    # ActionChains(driver).move_to_element(element).perform()


get_exons()
#genome_segment_extract("RHO")
