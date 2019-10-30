import numpy
from cogent3 import LoadTable


def format_group(group):
    group_dict = {'ENU_variants': 'ENU',
                  'germline_variants': 'Spontaneous'}
    group = group.split('/')[7]
    return group_dict[group]


def format_direction(direction):
    direction = direction.replace('to', r'$\rightarrow$')
    return direction


def format_positions(positions):
    positions = positions.split(":")
    old_new = {"0": "-2", "1": "-1", "2": "+1", "3": "+2"}
    redone = [old_new[p[-1]] for p in positions]
    if len(redone) == 1:
        template = "%s"
    else:
        template = "(%s)"

    return template % ", ".join(redone)


def latex_format_numbers(value):
    try:
        value = float(value)
    except ValueError:
        return value

    if value == 0:
        return "0"
    v = "%.1e" % value
    v = v.split("e")
    v[1] = str(int(v[1]))
    v = r"$%s\times 10^{%s}$" % tuple(v)
    return v


def format_pvalue(value):
    """0.0 if 0, float if <= 0.0001, scientific otherwise"""
    epsilon = 1e-6
    if value <= 1e-100:
        result = "0.0"
    elif 1 - epsilon <= value <= 1.0 + epsilon:
        result = "1.0"
    elif 1e-4 <= value:
        result = "%.4f" % value
    elif 0 < value < 1e-4:
        result = latex_format_numbers(value)
    else:
        raise ValueError(value)
    return result


def format_latex_table(table, justify, label=None):
    caption = table.title or ""
    table.title = ""
    label = label or ""
    result = table.to_string(format="latex", justify=justify)

    if caption:
        caption = r"\caption{%s}" % caption

    if label:
        caption += r"\label{%s}" % label

    if caption:
        result = result.splitlines()
        result.insert(-1, caption)
        result = "\n".join(result)

    return result


def classifier_summary_stats(table, stat, factors):
    """returns a table of summary statistics from the table"""
    dvs = table.distinct_values(factors)
    rows = []
    for dv in dvs:
        sub = table.filtered(lambda x: tuple(x) == dv, columns=factors)
        est = sub.tolist(columns=stat)
        est = numpy.array(est)
        row = list(dv) + [est.mean(), est.std(ddof=1), est.min(), est.max()]
        rows.append(row)

    t = LoadTable(header=list(factors) +
                  ['mean(%s)' % stat, 'std(%s)' % stat,
                   'min(%s)' % stat, 'max(%s)' % stat],
                  rows=rows, digits=3)
    t = t.sorted(factors)
    t = t.with_new_header(['size', 'name'], ['Training size', 'Feature set'])
    t.format_column('Training size', "{:,}".format)
    return t
