import numpy
import plotly.graph_objs as go
from plotly import subplots, figure_factory as ff

def MakeMatch(indexed_funcs, debug=False):
    def matches(vals):
        result = True
        for index, func in indexed_funcs.items():
            if not func(vals[index]):
                result = False
                break

        if debug:
            print(result, vals)
        return result
    return matches


def Coloriser(num, div='qual', comp='Paired'):
    CL = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    # dict of different colours
    coloured = {}

    def get_color(name, coloured=coloured):
        if name not in coloured:
            coloured[name] = CL[len(coloured)]
        c = coloured[name]
        return c
    return get_color


def get_layout(stat='AUC', width=600, height=400, title=None,
               label_font_size=20, margins=None, x_stacked=False,
               tickformat='.2f'):
    # xstacked means a subplots stacked along x-axis
    ylabel = f'<span style="text-decoration: overline">{stat}</span>'
    if margins is None:
        margins = dict(l=50, r=50, b=50, t=50)

    xattr = dict(title='Training size',
                 automargin=True,
                 titlefont=dict(size=label_font_size),
                 showgrid=True,
                 zeroline=False,
                 showline=True,
                 mirror=True,
                 tickformat='s')
    yattr = dict(title=ylabel,
                 automargin=True,
                 titlefont=dict(size=label_font_size),
                 showgrid=True,
                 zeroline=False,
                 showline=True,
                 mirror=True,
                 tickformat=tickformat)
    if x_stacked:
        xattr.pop('title')
        xattr['mirror'] = False
        yattr['mirror'] = False

    lo = go.Layout(autosize=False, width=width, height=height,
                   title=title,
                   yaxis=go.layout.YAxis(**yattr),
                   xaxis=go.layout.XAxis(**xattr),
                   margin=go.layout.Margin(**margins))
    return lo


def get_traces(table, stat='auc', category='name',
               order=("M", "M+I", "M+I+2D", "FS")):
    dvs = table.distinct_values(category)
    assert set(order) & dvs, dvs
    traces = []
    for element in order:
        subtable = table.filtered(lambda x: x == element, columns=category)
        subtable = subtable.sorted("size")
        if subtable.shape[0] == 0:
            continue
        data = subtable.tolist(columns=["size", f"mean_{stat}", f"std_{stat}"])
        x, y, std = zip(*data)
        # these should be 0 <= y + std <= 1.0
        arr_min = [std[i] if y[i] - std[i] > 0 else y[i] - 0
                   for i in range(len(y))]
        arr_max = [std[i] if y[i] + std[i] <= 1 else 1 - y[i]
                   for i in range(len(y))]
        trace = go.Scatter(
            x=x,
            y=y,
            error_y=dict(array=arr_max,
                         arrayminus=arr_min,
                         symmetric=False,
                         type='data'),
            name=element)
        traces.append(trace)

    return traces


def get_plotly_fig(traces, stat='AUC', title=None, layout=None):
    lo = layout or get_layout(title=title, stat=stat)
    fig = go.Figure(layout=lo)

    for trace in traces:
        fig.add_trace(trace)
    return fig


def stacked_fig_from_traces(traces, titles, label_font_size=20,
                            y_bounds=None, layout=None):
    fig = subplots.make_subplots(rows=1, cols=2, shared_xaxes=True,
                              shared_yaxes=True, vertical_spacing=0.1,
                              horizontal_spacing=0.05,
                              subplot_titles=titles,
                              print_grid=False)
    if layout:
        fig["layout"].update(layout)

    for col, traced in enumerate(traces):
        get_color = Coloriser(len(traced))
        for trace in traced:
            clr = get_color(trace["name"])
            trace["marker"].color = clr
            trace["error_y"].color = clr
            if col != 0:
                trace["showlegend"] = False
            fig.append_trace(trace, 1, col + 1)

    fig["layout"].update(
        yaxis=dict(
            range=y_bounds,
            title='<span style="text-decoration: overline">AUC</span>',
            tickformat=".2f",
        ),
        xaxis1=dict(showline=True),
        xaxis2=dict(showline=True),
        height=500,
        width=700)
    fig.layout.annotations += ({'font': {'color': '#000000',
                                         'size': label_font_size},
                                'text': 'Training size',
                                'x': 0.5,
                                'xanchor': 'center',
                                'xref': 'paper',
                                'y': -0.1,
                                'borderpad': 4,
                                'showarrow': False,
                                'yanchor': 'middle',
                                'yref': 'paper'},)
    return fig


def fig_grid_from_traces(traces, titles, label_font_size=20,
                         y_bounds=None, layout=None):
    fig = subplots.make_subplots(rows=2, cols=2, shared_xaxes=True,
                              shared_yaxes=True, vertical_spacing=0.1,
                              horizontal_spacing=0.05,
                              subplot_titles=titles,
                              print_grid=False)
    if layout:
        fig["layout"].update(layout)

    for (row, col), traced in traces.items():
        get_color = Coloriser(len(traced))
        for trace in traced:
            clr = get_color(trace["name"])
            trace["marker"].color = clr
            trace["error_y"].color = clr
            if (row, col) != (1, 1):
                trace["showlegend"] = False
            fig.append_trace(trace, row, col)

    fig["layout"].update(
        xaxis1=dict(showline=True),
        xaxis2=dict(showline=True),
        height=500,
        width=700)
    fig.layout.annotations += ({'font': {'color': '#000000',
                                         'size': label_font_size},
                                'text': 'Training size',
                                'x': 0.5,
                                'xanchor': 'center',
                                'xref': 'paper',
                                'y': -0.1,
                                'showarrow': False,
                                'yanchor': 'middle',
                                'yref': 'paper'},
                               {'font': {'color': '#000000',
                                         'size': label_font_size},
                                'text': '<span style="text-decoration: overline">AUC</span>',
                                'y': 0.5,
                                'xanchor': 'center',
                                'xref': 'paper',
                                'x': -0.1,
                                'showarrow': False,
                                'yanchor': 'middle',
                                'borderpad': 4,
                                'textangle': 270,
                                'yref': 'paper'})
    return fig


def faceted_fig(df, x, y, facet_col, y_bounds=None):
    fig = ff.create_facet_grid(df, x=x, y=y,
                               facet_col=facet_col,
                               facet_col_labels='name',
                               facet_row_labels='name')
    return fig


def get_genome_barchart(chrom_aucs, chrom_order, train_chrom,
                        title=None, y_bounds=None):
    assert train_chrom in chrom_order
    y_bounds = y_bounds or [0.0, 1.0]
    x = [str(c) for c in chrom_order if c != train_chrom]
    y = [chrom_aucs[c] for c in chrom_order if c != train_chrom]
    trace1 = go.Bar(
        x=x,
        y=y,
        name='Not trained',
        marker=dict(color='rgb(49,130,189)')
    )

    trace = go.Bar(
        x=[str(train_chrom)],
        y=[chrom_aucs[train_chrom]],
        name='Trained',
        marker=dict(color='rgb(204,204,204)')
    )
    fig = go.Figure(data=[trace, trace1])
    d = numpy.array([v for k, v in chrom_aucs.items() if k != 1])
    mean = d.mean()
    std = d.std(ddof=1)
    txt = dict(text='<span style="text-decoration: overline'
               '">AUC</span>=%.2f (±%.2f)' % (mean, std),
               showarrow=False, x=18, y=y_bounds[1] * 0.98)
    fig["layout"].update(yaxis=dict(range=y_bounds,
                                    tickformat=".2f",
                                    showline=True,
                                    title='AUC'
                                    ),
                         xaxis=dict(showline=True,
                                    title='Chromosome',
                                    dtick=1,
                                    type="category"),
                         height=500,
                         width=700,
                         title=title,
                         annotations=[txt])

    return fig


def get_fig_latex(path, caption, label):
    """returns latex figure section"""
    template = '''\\begin{figure}[!ht]
              \\centering
              \\includegraphics[width=1.0\\textwidth]{%(path)s}
              \\caption{%(caption)s}
              \\label{%(label)s}
            \\end{figure}
            '''
    data = dict(path=path, caption=caption, label=label)
    return template % data
