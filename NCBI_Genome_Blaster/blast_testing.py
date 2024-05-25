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


def clicking_wrapper(overall_driver, specific_driver, by_what, description,
                     timer_limit) -> None:
    """
    Clicks an element on the NCBI datasets page. Keeps attempting until
    timer_limit is reached or success occurs. On guard for a survey pop-up,
    which the program immediately closes upon detection. Note that this
    functions very similarly to selenium's built in WebDriverWait class, but
    has additional functionality tailored to NCBI

    :param overall_driver: general driver used
    :param specific_driver: a more specific driver that limits the scope of
    elements considered
    :param by_what: the type of the description (ex. By.ID, By.XPATH)
    :param description: descriptor of the element to click, as a string
    :param timer_limit: stop attempting once this much time has passed.
    :return:
    """

    counter = 0
    while counter < timer_limit:
        try:
            # finds the element and attempts to click
            specific_driver.find_element(by_what, description).click()
            print("successful click")
            break
        except WebDriverException:
            try:
                # checks if survey has popped up, not allowing clicking
                overall_driver.find_element(By.XPATH, "//button[text()[. = 'No Thanks']]").click()
                print("survey denied")
            except WebDriverException:
                pass

            # pauses for a period of time
            time.sleep(WAIT_CONSTANT)
            counter += WAIT_CONSTANT

    if counter >= timer_limit:
        # throw an exception if time limit is reached
        raise DriverTimeoutException("Timeout error for clicking wrapper, description:" +
                        description)


def get_element(specific_driver, by_what, description, timer_limit) -> object:
    """
    Wait until an element is present, then get the element as an object.

    :param specific_driver: driver, which allows for narrowing of the scope of
    elements searched for
    :param by_what: By.ID, By.XPATH, etc. the type of the descriptor used to
    search for the element
    :param description: the description of the element, ex. its xpath
    :param timer_limit: amount of time to attempt for
    :return: the element as an object
    """

    counter = 0
    while counter < timer_limit:
        try:
            element = specific_driver.find_element(by_what, description)
            return element
        except WebDriverException:
            time.sleep(WAIT_CONSTANT)
            counter += WAIT_CONSTANT

    # throws an exception if nothing is got within timer limit
    raise DriverTimeoutException("Timeout error for getting element, description: " +
                    description)


def try_click(specific_driver, by_what, description) -> bool:
    """
    Attempts a click, returning True if successful and False otherwise.

    :param specific_driver: driver
    :param by_what: type of description, ex. By.ID, By.XPATH
    :param description: description of element to click, in string form
    :return: True iff the click was successful
    """

    try:
        specific_driver.find_element(by_what, description).click()
        return True
    except WebDriverException:
        return False


def blast_menu_button_clicker(driver, timer_limit) -> bool:
    """
    Clicker function specifically for the menu button "BLAST against this
    genome. The reason why it requires its own function, rather than just the
    general clicker, is because sometimes the button will not appear as an
    option. TO verify that the menu has loaded but the button is not present,
    checks to see if the "view Details" button has loaded

    :param driver: driver
    :param timer_limit: time to keep attempting until an exception is thrown
    :return: True iff "Blast against this genome" is successfully clicked
    """

    counter = 0
    while counter < timer_limit:
        try:
            # tries clicking the "BLAST against this genome" button
            clicking_wrapper(driver, driver, By.XPATH, "//li[text()[. = 'BLAST against this genome']]", 4)
            return True
        except DriverTimeoutException:
            # if the View Details button is present, means that the "BLAST
            # against this genome" button is not available
            if EC.element_to_be_clickable((By.XPATH, "//li[text()[. = 'View details']]")):
                return False

        time.sleep(WAIT_CONSTANT)
        counter += WAIT_CONSTANT

    # throws exception if timer runs out
    raise DriverTimeoutException("blast menu button clicking timeout")


def xml_download_clicker(driver) -> None:
    """
    Clicker specifically to download the results of a BLAST file. After clicking
    "Download All", the "XML" button can be finicky to click, so the menu is
    looked over and over again until success is reached.

    :param driver: driver
    :return: None
    """

    time.sleep(1)
    while not try_click(driver, By.XPATH, "//a[text()[. = 'XML']]"):
        time.sleep(1)
        try_click(driver, By.XPATH,
                  "//a[text()[. = 'Download All']]")
        time.sleep(1)


def get_downloaded_xml_file(path: str) -> str:
    """
    Waits until the XML file is downloaded then returns its file path
    :param path: the path of the directory that will store the blast results
    (the XML file will be temporary and discorded soon after)
    :return: the file path of the XML file
    """

    downloaded = False
    while not downloaded:
        for filepath in os.listdir(path):
            if not os.path.isdir(filepath) and os.path.splitext(filepath)[1] == ".xml":
                return path + "\\" + filepath


def juggle_blast_tabs(driver, by_what, description, num_tabs,
                      sleep_time) -> int:
    """
    Switches between num_tabs tabs until the feature specified by the parameters
    is clickable. For this program, the only tabs that will be "juggled" are
    BLAST tabs while the results haven't finished loading. It is assumed the
    tabs to be juggled are in positions 1 - num_tabs (the NCBI datasets home
    page is position 0.

    :param driver: driver
    :param by_what: type of the description
    :param description: description of the element to click
    :param num_tabs: number of tabs to juggle between
    :param sleep_time:
    :return: the position of the tab that has loaded
    """

    switch_counter = 1
    driver.switch_to.window(driver.window_handles[1])

    while not try_click(driver, by_what, description):
        switch_counter = (switch_counter % num_tabs) + 1
        driver.switch_to.window(driver.window_handles[switch_counter])

        time.sleep(sleep_time)

    return switch_counter


if __name__ == "__main__":

    # this part will be in the user-inputted main program
    user_specified_path = r'C:\Users\tonyx\Downloads'
    # directory of files containing query sequences to blast. The taxa to blast
    # against are in the title
    rt_directory_path = r'C:\Users\tonyx\Downloads\test_taxa4' # CSV file

    path = user_specified_path + "\\" + "blast_test_results"

    # if the blast downloads folder already exists, use a number incrementing
    # system to name a new one
    name_numbering = 1
    while os.path.exists(path):
        path = path + str(name_numbering)
        name_numbering += 1

    os.mkdir(path)

    # set the webdriver's download path, allow automatic downloads and safe
    # browsing
    chrome_options = webdriver.ChromeOptions()
    prefs = {'download.default_directory': path,
             'profile.default_content_setting_values.automatic_downloads': 1,
             "safebrowsing.enabled": True}
    chrome_options.add_experimental_option('prefs', prefs)
    driver = webdriver.Chrome(chrome_options)

    # access the gene search page, which will have the refseq data
    driver.get("https://www.ncbi.nlm.nih.gov/datasets/genome/")

    # a dictionary that acts as a hash table to check if species have been
    # blasted before
    species_so_far = {}

    # another file, specifying the order of the taxa to be searched in
    # this is NON-trivial, since once a species has been BLASTED, it will not
    # be BLASTED again, even if a different query sequence is being used.
    ordered_taxa_path = r"C:\Users\tonyx\Downloads\ordered_taxa.txt"
    taxa_in_order = open(ordered_taxa_path, "r").readlines()

    refined_taxa_in_order = []
    for line in taxa_in_order:
        if line.strip() != "":
            refined_taxa_in_order.append(line.strip())

    for taxa in refined_taxa_in_order:

        # reference query path
        reference_filepath = rt_directory_path + "\\" + taxa + ".fas"

        # get the search bar element, click on it, and delete anything currently
        # in it
        search_bar = get_element(driver, By.ID, "taxonomy_autocomplete", 40)
        clicking_wrapper(driver, driver, By.ID, "taxonomy_autocomplete", 40)
        search_bar.send_keys(Keys.BACKSPACE)
        # sometimes the page changes
        search_bar = driver.find_element(By.ID, "taxonomy_autocomplete")
        # search with the taxa (name of the query file)
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
        # tracks which tabs correspond to which species (not directly evident on
        # the blast page
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
                # DO NOT BLAST AGAIN if a species has already been BLASTed
            else:
                # open the menu
                print("on: " + taxon_name)
                species_so_far[taxon_name] = True

                clicking_wrapper(driver, current_row, By.XPATH,
                                 ".//button[@data-ga-action='click_open_menu']", 40)
                # wait until the menu item shows up
                # if the species is blastable
                if blast_menu_button_clicker(driver, 40):

                    tabs_open += 1
                    species_open.append(taxon_name)
                    print(species_open)
                    driver.switch_to.window(driver.window_handles[tabs_open])

                    clicking_wrapper(driver, driver, By.XPATH, "//*[text()[. = 'Somewhat similar sequences (blastn)']]", 40)
                    file_input = get_element(driver, By.ID, "upl", 40)
                    file_input.send_keys(reference_filepath)
                    clicking_wrapper(driver, driver, By.CLASS_NAME, "blastbutton", 40)

                # if the maximum number of blast entries has been reached
                if tabs_open == MAX_NUM_TABS:

                    # wait until the results of a tab show up
                    tab_no = juggle_blast_tabs(driver, By.XPATH,
                                               "//a[text()[. = 'Download All']]",
                                               MAX_NUM_TABS, 2)
                    print("tab no: " + str(tab_no))

                    # download the xml file for the blast results, get its path
                    xml_download_clicker(driver)
                    file_to_analyze = get_downloaded_xml_file(path)

                    # parse the xml file, creating files that contain the
                    # results in fasta format
                    parse_blast_xml(file_to_analyze, path, taxa,
                                    species_open[tab_no - 1])
                    os.remove(file_to_analyze)

                    driver.close()
                    tabs_open -= 1
                    species_open.pop(tab_no - 1)
                    print(species_open)

                driver.switch_to.window(driver.window_handles[0])
                ActionChains(driver).move_by_offset(10, 10).click().perform()

        # after blast searches have been run on every available species in the
        # taxa, winds down and collects results from the remaining tabs that
        # are running
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





