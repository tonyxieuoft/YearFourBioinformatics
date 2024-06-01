from selenium import webdriver
from selenium.webdriver.common.by import By
from typing import Dict

PHYLOGENY_COLUMN = 3
SPECIES_NAME_COLUMN = 1
CODE_COLUMN = 0


def get_taxonomy_lineage(species_string: str) -> Dict:
    """
    Given a string of species, returns complete lineages for each species as
    values in a dictionary, where the names of the species are the keys.

    :param species_string: a string of species names, separated by newlines
    :return: a dictionary, where keys are species names and values are lists of
    lineages.
    """
    # assemble a new driver with the "headless" option turned on
    chrome_options = webdriver.ChromeOptions()
    chrome_options.add_argument('--headless')
    chrome_options.add_argument('--disable-gpu')
    driver = webdriver.Chrome(chrome_options)

    # access the gene search page, which will have the refseq data
    # rate limiting factor
    driver.get("https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi")

    # text box for entering species' names
    driver.find_element(By.NAME, "tax").send_keys(species_string)

    # option for complete phylogeny
    driver.find_element(By.NAME, "lng").click()

    driver.find_element(By.XPATH, "//input[@value='Show on screen']").click()
    # table of complete phylogenies
    table = driver.find_element(By.XPATH, "//h2")
    # each row is for a separate species
    rows = table.find_elements(By.XPATH, ".//tr")
    phylo_dict = {}
    for row in rows[1:]:
        columns = row.find_elements(By.XPATH, ".//td")
        # if columns[CODE_COLUMN] == "1":
        phylo_dict[columns[SPECIES_NAME_COLUMN].text] = \
            columns[PHYLOGENY_COLUMN].text.split(" ")

    return phylo_dict


