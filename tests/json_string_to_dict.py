import json

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


if __name__ == "__main__":

    json_file = r"C:\Users\tonyx\Downloads\elasmo_genome_info.txt"
    string = open(json_file).read().strip()
    dct = json.loads(string)
    print_dict(dct, "")




