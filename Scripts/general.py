def findKeysByValue(dict_to_check, item_to_find, compar_operator):
    keys = [k for (k, v) in dict_to_check.items() if compar_operator(v, item_to_find)]
    return sorted(keys)


def intersection(a, b):
    out = []
    for el in a:
        if (el in b) & (el not in out):
            out.append(el)
    return out


def substraction(a, *args):
    out = []
    for el in a:
        for b in args:
            if (el not in b) & (el not in out):
                out.append(el)
    return out
