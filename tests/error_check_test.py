import time

from selenium import webdriver
from selenium.common import WebDriverException
from selenium.webdriver import ActionChains, Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

driver = webdriver.Chrome()
# access the gene search page, which will have the refseq data
driver.get("https://www.ncbi.nlm.nih.gov/datasets/genome/")

while True:
    try:
        driver.find_element(By.ID, "fake_id")
    except WebDriverException:
        print("oops")
    finally:
        time.sleep(2)


# it's just a button called "No Thanks"

