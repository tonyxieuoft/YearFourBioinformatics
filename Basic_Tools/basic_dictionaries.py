from typing import Dict, Iterable, List
import json


def dict_get_values(dct: Dict) -> List:

    encountered_values = {}
    for key in dct.keys():
        if isinstance(dct[key], Iterable):
            for item in dct[key]:
                if item not in encountered_values:
                    encountered_values[item] = True
        else:
            if dct[key] not in encountered_values:
                encountered_values[dct[key]] = True

    return list(encountered_values.keys())


def print_dict(obj, spacer):

    if isinstance(obj, dict):
        for key in obj:
            print(spacer + str(key))
            print_dict(obj[key], spacer + "  ")

    elif isinstance(obj, list):
        print(spacer + "List:")
        counter = 1
        for item in obj:
            print(spacer + "Item " + str(counter))
            print_dict(item, spacer + "  ")
            counter += 1

    else:
        print(spacer + str(obj))


def json_to_dict(json_file):

    string = open(json_file).read().strip()
    return json.loads(string)
