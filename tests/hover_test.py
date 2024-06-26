from Bio import Entrez
from selenium import webdriver
from selenium.common import JavascriptException, StaleElementReferenceException, \
    WebDriverException
from selenium.webdriver import ActionChains
from selenium.webdriver.common.by import By
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC



def get_taxid(species: str):
    """
    Given a taxon name, return the list of orders within the taxon

    Input: taxon name
    Output: array of orders
    """
    # Open selenium and go to the NCBI taxonomy browser
    driver = webdriver.Chrome()
    driver.get("https://www.musicalmindsco.ca")

    WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.ID, "comp-jk0h0qvh2")))
    total = driver.find_element(By.ID, "comp-jk0h0qvh2")
    elements = total.find_elements(By.XPATH, "//*")

    print(len(elements))
    # Search the taxonomy name
    for element in elements:
        try:
            ActionChains(driver).move_to_element(element).perform()
        except JavascriptException:
            pass
        except StaleElementReferenceException:
            pass
    #WebDriverWait(driver, 1000).until(EC.presence_of_element_located((By.ID, "hello")))

if __name__ == "__main__":

    Entrez.email = "xiaohan.xie@mail.utoronto.ca"
    gene_table_file = Entrez.efetch(db='gene', id="122563057", rettype='gene_table',
                                    retmode="text")
    table_in_text = str(gene_table_file.read())
    print(table_in_text)
