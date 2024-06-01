from typing import Dict, Iterable, List


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
