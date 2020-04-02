from collections import OrderedDict
import re
import os

from cogent3 import LoadTable

ref = re.compile(r"\\ref\{[^\}]+\}")
label = re.compile(r"\\label{\S+}")


def get_tags(filepath, tag):
    """returns the set of tags from a latex file.

    tag can be either label or ref"""
    pattern = {"label": label, "ref": ref}[tag]
    with open(filepath) as infile:
        data = "".join(infile.readlines())

    all_tags = pattern.findall(data)
    all_tags = [t.split("{")[-1][:-1] for t in all_tags]
    #  we ignore equation tags
    all_tags = [
        (t.strip(), filepath.strip()) for t in all_tags if not t.startswith("eq:")
    ]
    return OrderedDict(all_tags)


def split_main_supp_tags(tags):
    """splits the labels into those for the main manuscript and
    for the supplementary"""
    supp_tags = OrderedDict()
    main_tags = OrderedDict()
    for tag in tags:
        store = supp_tags if tag.startswith("sup") else main_tags
        store[tag] = tags[tag]
    return main_tags, supp_tags


def filtertags(filterfunc, tags):
    new_tags = OrderedDict()
    for tag, value in tags.items():
        if tag in new_tags:
            continue

        if filterfunc(tag):
            new_tags[tag] = value
    return new_tags


def get_ms_supp_labels(float_type, texdir="../ENU-ms-genetics-v2", verbose=False):
    """returns ordered dicts of labels from the manuscript for
    supplementary and main manuscript body"""
    #     assert float_type in ('fig', 'tab')
    # hardcoding these, in the manuscript order of the sections
    texfns = [
        os.path.join(texdir, tfn)
        for tfn in (
            "MS-introduction.tex",
            "MS-results.tex",
            "MS-discussion.tex",
            "MS-methods.tex",
        )
    ]
    alllabels = None
    for tfn in texfns:
        tags = get_tags(tfn, "label")
        if alllabels is None:
            alllabels = tags
        else:
            alllabels.update(tags)

    print("\n\nWorking on labels")
    alllabels = filtertags(float_type, alllabels)

    allrefs = None
    for tfn in texfns:
        tags = get_tags(tfn, "ref")
        if allrefs is None:
            allrefs = tags
        else:
            allrefs.update(tags)
    print("\n\nWorking on refs")
    allrefs = filtertags(float_type, allrefs)
    mainrefs = filtertags(lambda x: not x.startswith("sup"), allrefs)
    suprefs = filtertags(lambda x: x.startswith("sup"), allrefs)

    missing = set(alllabels) - set(mainrefs)

    rows = [(missed, alllabels[missed]) for missed in missing]
    table = LoadTable(header=["label missing", "referenced in"], rows=rows)
    if verbose:
        print(table)
    return mainrefs, suprefs
