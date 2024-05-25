import xml.etree.ElementTree as ET


def xml_to_dictionary(file):

    tree = ET.parse(file)
    root = tree.getroot()
    return {root.tag: xml_to_dictionary_recursive(root)}


def xml_to_dictionary_recursive(root: object) -> object:

    use = "dict"

    counter = 0
    last_tag = ""
    for child in root:
        if child.tag == last_tag:
            use = "array"
            break
        last_tag = child.tag
        counter += 1

    if counter == 0:
        use = "value"

    if use == "value":
        return root.text
    elif use == "dict":
        new_dict = {}
        for child in root:
            new_dict[child.tag] = xml_to_dictionary_recursive(child)
        return new_dict
    else:
        new_array = []
        for child in root:
            new_array.append(xml_to_dictionary_recursive(child))
        return new_array
