from selenium import webdriver
from selenium.common import WebDriverException
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
    driver.get("https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi")

    # Search the taxonomy name
    driver.find_element(By.ID, "searchtxt").send_keys(species)
    driver.find_element(By.NAME, "a").click()

    try:
        WebDriverWait(driver, 10).until(EC.presence_of_element_located(
            (By.XPATH, "//strong")))
        all_clickable = driver.find_elements(By.XPATH, "//strong")
        # line order
        taxon = all_clickable[1]
        taxon.click()
    except WebDriverException:
        # a species name was entered so the taxonomy browser went directly to
        # the taxonomy page
        taxon = None

    WebDriverWait(driver, 10).until(
        EC.presence_of_element_located((By.XPATH, "//small[contains(text("
                                                  "), 'for references in "
                                                  "articles please "
                                                  "use')]")))
    taxid_text = driver.find_element(By.XPATH,
                                     "//small[contains(text("
                                     "), 'for references in "
                                     "articles please "
                                     "use')]").text
    print(taxid_text)
    taxid = ""
    for i in range(0, len(taxid_text)):
        if taxid_text[i].isdigit():
            taxid = taxid + taxid_text[i]

    return taxid

get_taxid("Cetacea")
