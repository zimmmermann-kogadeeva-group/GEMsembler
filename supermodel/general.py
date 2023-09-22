def findKeysByValue(dict_to_check, item_to_find, compar_operator):
    keys = [k for (k, v) in dict_to_check.items() if compar_operator(v, item_to_find)]
    return sorted(keys)


def is_float(element: any) -> bool:
    # If you expect None to be passed:
    if element is None:
        return False
    try:
        float(element)
        return True
    except ValueError:
        return False
