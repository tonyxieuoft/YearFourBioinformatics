import os.path
import time
import os
from selenium import webdriver
from selenium.common import WebDriverException
from selenium.webdriver import ActionChains, Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from test_fileread import parse_blast_xml

WAIT_CONSTANT = 0.1
MAX_NUM_TABS = 2


class DriverTimeoutException(Exception):
    pass


def clicking_wrapper(overall_driver, specific_driver, by_what, description, timer_limit):

    counter = 0
    while counter < timer_limit:
        try:
            specific_driver.find_element(by_what, description).click()
            print("successful click")
            break
        except WebDriverException:
            try:
                overall_driver.find_element(By.XPATH, "//button[text()[. = 'No Thanks']]").click()
                print("survey denied")
            except WebDriverException:
                pass

            time.sleep(WAIT_CONSTANT)
            counter += WAIT_CONSTANT

    if counter >= timer_limit:
        raise DriverTimeoutException("Timeout error for clicking wrapper, description:" +
                        description)


def get_element(specific_driver, by_what, description, timer_limit):

    counter = 0
    while counter < timer_limit:
        try:
            element = specific_driver.find_element(by_what, description)
            return element
        except WebDriverException:
            time.sleep(WAIT_CONSTANT)
            counter += WAIT_CONSTANT

    raise DriverTimeoutException("Timeout error for getting element, description: " +
                    description)


def try_click(specific_driver, by_what, description):

    try:
        specific_driver.find_element(by_what, description).click()
        return True
    except WebDriverException:
        return False


def blast_menu_button_clicker(driver, timer_limit):

    counter = 0
    while counter < timer_limit:
        try:
            clicking_wrapper(driver, driver, By.XPATH, "//li[text()[. = 'BLAST against this genome']]", 4)
            return True
        except DriverTimeoutException:
            if EC.element_to_be_clickable((By.XPATH, "//li[text()[. = 'View details']]")):
                return False

        time.sleep(WAIT_CONSTANT)
        counter += WAIT_CONSTANT

    raise DriverTimeoutException("blast menu button clicking timeout")


def xml_download_clicker(driver):

    time.sleep(1)
    while not try_click(driver, By.XPATH, "//a[text()[. = 'XML']]"):
        time.sleep(1)
        try_click(driver, By.XPATH,
                  "//a[text()[. = 'Download All']]")
        time.sleep(1)


def get_downloaded_xml_file(path):

    downloaded = False
    while not downloaded:
        for filepath in os.listdir(path):
            if not os.path.isdir(filepath) and os.path.splitext(filepath)[1] == ".xml":
                return path + "\\" + filepath


def juggle_blast_tabs(driver, by_what, description, num_tabs, sleep_time):

    switch_counter = 1
    driver.switch_to.window(driver.window_handles[1])

    while not try_click(driver, by_what, description):
        switch_counter = (switch_counter % num_tabs) + 1
        driver.switch_to.window(driver.window_handles[switch_counter])

        time.sleep(sleep_time)

    return switch_counter


if __name__ == "__main__":

    user_specified_directory = r'C:\Users\tonyx\Downloads'
    rt_directory_path = r'C:\Users\tonyx\Downloads\test_taxa' # CSV file

    path = user_specified_directory + "\\blast_test_results"
    # make a folder beforehand, error check later to see if this folder exists

    name_numbering = 1
    while os.path.exists(path):
        path = path + str(name_numbering)
        name_numbering += 1

    os.mkdir(path)

    # set the webdriver's download path
    chrome_options = webdriver.ChromeOptions()
    prefs = {'download.default_directory': path,
             'profile.default_content_setting_values.automatic_downloads': 1,
             "safebrowsing.enabled": True}
    chrome_options.add_experimental_option('prefs', prefs)
    driver = webdriver.Chrome(chrome_options)

    # access the gene search page, which will have the refseq data
    driver.get("https://www.ncbi.nlm.nih.gov/datasets/genome/")

    species_so_far = {}

    ordered_taxa_path = r"C:\Users\tonyx\Downloads\ordered_taxa.txt"
    taxa_in_order = open(ordered_taxa_path, "r").readlines()
    refined_taxa_in_order = []
    for line in taxa_in_order:
        if line.strip() != "":
            refined_taxa_in_order.append(line.strip())
    print(refined_taxa_in_order)

    for taxa in refined_taxa_in_order:

        reference_filepath = rt_directory_path + "\\" + taxa + ".fas"

        search_bar = get_element(driver, By.ID, "taxonomy_autocomplete", 40)
        clicking_wrapper(driver, driver, By.ID, "taxonomy_autocomplete", 40)
        search_bar.send_keys(Keys.BACKSPACE)
        # sometimes the page changes
        search_bar = driver.find_element(By.ID, "taxonomy_autocomplete")
        search_bar.send_keys(taxa)
        search_bar.send_keys(Keys.RETURN)

        # clicks the filter menu
        clicking_wrapper(driver, driver, By.ID, "panel-header", 40)
        # selects for only reference genomes to be displayed
        clicking_wrapper(driver, driver, By.XPATH, "//*[@value='reference_only']", 40)
        # waits for the reference genome to be displayed
        time.sleep(2)

        # clicks the "number of results displayed" menu button
        clicking_wrapper(driver, driver, By.CSS_SELECTOR,
                         ".MuiSelect-select.MuiSelect-outlined.MuiInputBase-"
                         "input.MuiOutlinedInput-input.css-zcubqt", 40)
        # selects the maximum number of results that can be displayed (100)
        clicking_wrapper(driver, driver, By.XPATH, "//li[@data-value='100']", 40)

        time.sleep(2)
        # gets the total number of rows to iterate over
        WebDriverWait(driver, 40).until(
            EC.presence_of_element_located((By.CSS_SELECTOR, ".MuiTableRow-root.css-1mc2v0b")))
        table_rows = driver.find_elements(By.CSS_SELECTOR, ".MuiTableRow-root.css-1mc2v0b")

        tabs_open = 0
        species_open = []

        for i in range(1, len(table_rows)):

            # waits for table results to show up, then finds all of them and -
            # puts them in an array
            WebDriverWait(driver, 40).until(
                EC.presence_of_element_located((By.CSS_SELECTOR, ".MuiTableRow-root.css-1mc2v0b")))
            table_rows = driver.find_elements(By.CSS_SELECTOR, ".MuiTableRow-root.css-1mc2v0b")
            current_row = table_rows[i]

            taxon_name = get_element(current_row, By.XPATH,
                                     ".//a[@data-ga-label='taxon_name']",
                                     40).text

            if taxon_name in species_so_far.keys():
                print("already searched: " + taxon_name)
            else:
                # open the menu
                print("on: " + taxon_name)
                species_so_far[taxon_name] = True

                clicking_wrapper(driver, current_row, By.XPATH,
                                 ".//button[@data-ga-action='click_open_menu']", 40)
                # wait until the menu item shows up

                if blast_menu_button_clicker(driver, 40):

                    tabs_open += 1
                    species_open.append(taxon_name)
                    print(species_open)
                    driver.switch_to.window(driver.window_handles[tabs_open])

                    clicking_wrapper(driver, driver, By.XPATH, "//*[text()[. = 'Somewhat similar sequences (blastn)']]", 40)
                    file_input = get_element(driver, By.ID, "upl", 40)
                    file_input.send_keys(reference_filepath)
                    clicking_wrapper(driver, driver, By.CLASS_NAME, "blastbutton", 40)

                if tabs_open == MAX_NUM_TABS:

                    tab_no = juggle_blast_tabs(driver, By.XPATH,
                                               "//a[text()[. = 'Download All']]",
                                               MAX_NUM_TABS, 2)
                    print("tab no: " + str(tab_no))

                    xml_download_clicker(driver)
                    file_to_analyze = get_downloaded_xml_file(path)

                    parse_blast_xml(file_to_analyze, path, taxa,
                                    species_open[tab_no - 1])
                    os.remove(file_to_analyze)

                    driver.close()
                    tabs_open -= 1
                    species_open.pop(tab_no - 1)
                    print(species_open)

                driver.switch_to.window(driver.window_handles[0])
                ActionChains(driver).move_by_offset(10, 10).click().perform()

        for i in range(tabs_open, 0, -1):

            tab_no = juggle_blast_tabs(driver, By.XPATH,
                              "//a[text()[. = 'Download All']]",
                              i, 2)

            print("tab no: " + str(tab_no))

            xml_download_clicker(driver)
            file_to_analyze = get_downloaded_xml_file(path)
            parse_blast_xml(file_to_analyze, path, taxa, species_open[tab_no - 1])
            os.remove(file_to_analyze)
            species_open.pop(tab_no - 1)

            print(species_open)
            driver.close()

        driver.switch_to.window(driver.window_handles[0])
        time.sleep(0.5)





