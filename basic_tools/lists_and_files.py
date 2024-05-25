import os
import time
from typing import List, Dict, Tuple
from Bio import Entrez
import re


def file_to_list(filepath: str) -> List:
    """
    Converts the contents of a file into a list by line

    :param filepath: the path of the file to convert
    :return: a list with each index storing a line of the file
    """

    f = open(filepath)
    raw = f.readlines()
    output = []
    for i in range(0, len(raw)):

        # strips to remove any potential whitespace
        if raw[i].strip() != "":
            output.append(raw[i].strip())

    return output


def list_to_string_csv(lst: List[object]) -> str:
    """
    Converts a list to a string, with each index separated by a comma
    :param lst: the list to convert
    :return: a comma-separated string
    """
    first = True
    string = ""
    for item in lst:
        if first:
            string = item
            first = False
        else:
            string = string + "," + item

    return string
