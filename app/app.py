#!/usr/bin/env python3
# -*- coding=utf-8 -*-

"""
## Documentation

This application aims to visualize local geometric informations about a
molecular structure using the 
[`pychemcurv`](https://pychemcurv.readthedocs.io/) package. In particular,
the application computes geometric quantities which provide an insight of the 
discrete curvature of a molecular or periodic structure. 

Custom data can be visualized by editing manually the table.

### Global overview

The dashed box on the top of the page
allows you to upload an xyz file. Click into the box or drag and
drop your file there. The page is splitted in three parts. On the left, you can 
visualize your structure, on the rigth, the selected data are plotted and below
a that gathers the data.

#### On the left

* The *"Select data"* dropdown, allows to select the data you want to map on 
the structure.
* The *"Select colormap"* dropdown, allows you to select a colormap. The `_r` 
label corresponds to colormap in reverse order.
* The *"bounds"* inputs change the min and max values used to 
compute the colors associated to the data.
* The *"Nan color"* input, can be used to set a color to atoms for which the
selected data does not exist.

#### On the right

By default, the right panel displays an histogram of the selected data, used
on the structure visualization. On top of this plot, a box plot presents an
overview of the distribution. Below the plot, a table gathers statistical
information of the data. In that case, the slider allows you to change
the number of bins of the histogram.

The dropdown menu *"Hisogram or abscissa"* allows you to plot either an
histogram of the selected data (default) or to plot the selected data as 
a function of another data. In that case, a trend line is also plotted.
Statistical information of both data are then 
displayed in the table below the plot.

#### Data table

Below the visualization part, a table displays all the data provided by the
`pychemcurv` package. Select the columns you want to see using the dropdown
menu. All the table is editable, manually, and the 
visualization is updated each time you modify a value.

If you want to add manualy custom data, you can add the `custom` column to the
table and fill it with your values. You can copy and paste data from a 
spreadsheet or a text file.

The whole data can be downloaded in csv format from the `export` button at the 
top.

**Warning:** If you edit the data in the table, you have first to refresh the
application before uploading a new molecule.

### Geometrical data

All the definitions of the atomic quantities available in this application are
defined in details in
[this publication](https://hal.archives-ouvertes.fr/hal-02490358/document) 
or are briefly described in the
[pychemcurv documentation](https://pychemcurv.readthedocs.io/en/latest/).

Hereafter is a quick list that gives the units of the quantities:

`atom_idx`
: index of the atom in the system, starting from 0

`species`
: chemical element as provided in the xyz file

`pyrA`
: pyramidalization angle in degrees

`angular_defect`
: angular defect in degrees

`n_star_A`
: number of atoms bonded to this atom

`spherical_curvature`
: spherical curvature, no unit

`improper`
: improper angle in degrees

`pyr_distance`
: distance of atom A from the plane defined from atoms bonded to A

`atom_A`
: coordinates of atom A

`star_A`
: coordinates of atom bonded to A

`hybridization`
: hybridization as define by Haddon et al., n tilde, amount of pz AO in 
the system sigma

`m`
: `m = (c_pi / lambda_pi)^2`

`n`
: `n = 3m + 2`

`c_pi^2`
: c_pi is the coefficient of the s AO in the h_pi hybrid orbital

`lambda_pi^2`
: lambda_pi is the coefficient of the p_pi AO in the h_pi hybrid orbital

`ave. neighb. dist.`
: Average distance with neighbors of atom A.


### File and data upload

The application accepts standard xyz files.
Such a file is suposed to display the number of atoms on the first line,
followed by a title line and followed by the structure in cartesian
coordinates. Each line contains the element as first column and the
cartesian coordinates as 2d, 3th and 4th columns, for example:

    3
    H2O molecule
    O   -0.111056  0.033897  0.043165
    H    0.966057  0.959148 -1.089095
    H    0.796629 -1.497157  0.403985

Coordinates have to be in angstrom to determine correctly the connectivity.

"""

import os
import base64
import re
import yaml

import dash
from dash import html, dcc, dash_table
from dash.dash_table.Format import Format, Scheme
from dash.dependencies import Input, Output, State
import plotly.graph_objs as go
import plotly.express as px
import dash_bio

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from pymatgen.core import Molecule
from pychemcurv import CurvatureAnalyzer

__author__ = "Germain Salvato Vallverdu"
__title__ = "Pychemcurv data viewer"
__subtitle__ = "Part of the Mosaica project"

HEX_COLOR_PATT = re.compile(r"^#[A-Fa-f0-9]{6}$")

# ---- Set up App ----
ext_css = ["https://use.fontawesome.com/releases/v5.8.1/css/all.css"]
app = dash.Dash(__name__,
                external_stylesheets=ext_css,
                # url_base_pathname="/mosaica/",
                suppress_callback_exceptions=True)
server = app.server

# with open(app.get_asset_url("data/elementColors.yml"), "r") as fyml:
with open("assets/data/elementColors.yml", "r") as fyml:
    ELEMENT_COLORS = yaml.load(fyml, Loader=yaml.SafeLoader)["jmol"]

#
# Layout
# ------------------------------------------------------------------------------

# --- define tab style
tab_style_header = {
    'backgroundColor': 'white',
    "padding": "5px",
    'fontWeight': 'bold',
    "textAlign": "center",
    "borderBottom": "2px solid rgb(60, 93, 130)",
    "borderTop": "2px solid rgb(60, 93, 130)",
    "fontFamily": "sans-serif"
}

style_data_conditional = [
    {'if': {'row_index': 'odd'}, 'backgroundColor': 'rgba(60, 93, 130, .05)'}
]

# --- header ---
header = html.Div(className="head", children=[
    html.Div(className="container", children=[
        html.H1(children=[html.Span(className="fas fa-atom"), " ", __title__]),
        # html.H2(__subtitle__)
        html.A(
            id="github-link",
            href="https://github.com/gVallverdu/pychemcurv",
            children=[
                "View on GitHub",
            ]
        ),
        html.Span(id="github-icon", className="fab fa-github fa-2x"),
    ])
])

# --- Footer ---
footer = html.Div(className="foot", children=[
    html.Div(className="container", children=[
        html.Div(className="row", children=[
            html.Div(className="eight columns", children=[
                html.H5("About:"),
                html.A("Germain Salvato Vallverdu",
                       href="https://gsalvatovallverdu.gitlab.io/"),
                html.Br(),
                html.A("University of Pau & Pays Adour",
                       href="https://www.univ-pau.fr")
            ]),
            html.Div(className="four columns", children=[
                html.A(href="https://www.univ-pau.fr", children=[
                    html.Img(
                        src=app.get_asset_url("img/LogoUPPAblanc.png"),
                    )
                ])
            ])
        ]),
    ])
])

# --- Body: main part of the app ---
body = html.Div(className="container", children=[

    # --- store components for the data
    dcc.Store(id="data-storage", storage_type="memory"),

    # --- upload part
    html.Div(className="row", id="top-panel", children=[
        # --- upload xyz file
        html.Div(className="four columns", children=[
            dcc.Upload(
                id='file-upload',
                children=html.Div(
                    className="upload-area control",
                    children="Upload xyz file here"
                ),
            ),
        ]),

        # --- intro text
        html.Div(className="eight columns", children=[
            dcc.Markdown("""
            The [documentation is available at the bottom of this page](#documentation). 

            Upload an xyz file on the left. The structure will appear
            on the left and the data can be plotted on the right.
            A [data table is available below](#data-table-title). 

            """)
        ])
    ]),

    html.Div(className="row", children=[
        # --- left panel: 3D visualization
        html.Div(id="left-panel", children=[

            # --- dash bio Molecule 3D Viewer
            html.Div(id="dash-bio-viewer"),

            # --- color bar
            dcc.Graph(id='colorbar', config=dict(displayModeBar=False)),

            # --- controls of mapped data
            html.Div(className="row", children=[

                # --- select data
                html.Div(className="six columns", children=[
                    html.Span("Select data", className="control-label"),
                    dcc.Dropdown(id='dropdown-data',
                                 placeholder="Select data"),
                ]),

                # --- colormap
                html.Div(className="six columns", children=[
                    html.Span("Colormap", className="control-label"),
                    dcc.Dropdown(
                        id='dropdown-colormap',
                        options=[{"label": cm, "value": cm}
                                 for cm in plt.colormaps()],
                        value="cividis"
                    ),

                    # --- colormap boundaries
                    html.Div(className="row", children=[
                        html.Div(className="four columns", children=[
                            html.Span(
                                "bounds", className="control-label",
                                style={"lineHeight": "38px"}
                            )
                        ]),
                        html.Div(className="four columns", children=[
                            dcc.Input(
                                id="cm-min-value", type="number", debounce=True,
                                placeholder="min", style={"width": "100%"}),
                        ]),
                        html.Div(className="four columns", children=[
                            dcc.Input(
                                id="cm-max-value", type="number", debounce=True,
                                placeholder="max", style={"width": "100%"}),
                        ]),
                    ], style={"marginTop": "10px"}),

                    # --- nan color selector
                    html.Div(className="row", children=[
                        html.Div(className="four columns", children=[
                            html.Span(
                                "Nan color", className="control-label",
                                style={"lineHeight": "38px"}
                            )
                        ]),
                        html.Div(className="eight columns", children=[
                            dcc.Input(
                                id="nan-color-value",
                                debounce=True,
                                placeholder="#000000",
                                type="text",
                                pattern=u"^#[A-Fa-f0-9]{6}$",
                                style={"width": "100%"},
                            ),
                        ]),
                    ], style={"marginTop": "10px"})
                ]),
            ]),
        ]),

        # --- right panel: plot data
        html.Div(id="right-panel", children=[

            # --- plot figure
            dcc.Graph(id='plot-data'),

            # --- a table of the selected data
            dash_table.DataTable(
                id="plot-data-table",
                style_header=tab_style_header,
                style_data_conditional=style_data_conditional,
            ),

            html.Div(className="row", children=[
                # --- select histogram or absicssa
                html.Div(className="six columns", children=[
                    html.Span("Histogram or abscissa",
                              className="control-label"),
                    dcc.Dropdown(
                        id="plot-data-selector",
                        options=[{"label": "histogram", "value": "histogram"}],
                        value="histogram",
                        placeholder="Plot more data",
                    ),
                ]),

                # --- number of bins for histogram
                html.Div(className="six columns", children=[
                    html.Span("# bins", className="control-label"),
                    dcc.Slider(
                        id="nbins-slider",
                        min=5, max=50, step=1,
                        value=30,
                        marks={i: "%d" % i for i in range(5, 51, 5)},
                    ),
                ]),
            ], style={"marginTop": "10px"}),
        ]),
    ]),

    # --- Data table
    html.Div(id="data-table-container", children=[
        html.H4("Data Table", id="data-table-title"),
        html.Div(className="column-selector-label",
                 children="Select the columns of the table:"),
        dcc.Dropdown(
            id="data-column-selector",
            multi=True,
            clearable=False,
        ),

        html.Div(children=[
            dash_table.DataTable(
                id="data-table",
                export_format='csv',
                export_columns="all",
                editable=True,
                style_header=tab_style_header,
                style_data_conditional=style_data_conditional,
            )
        ]),
    ]),

    # --- Documentation
    html.Div(className="documentation", children=[
        dcc.Markdown(__doc__, id="documentation")
    ])
])

app.layout = html.Div([header, body, footer])

#
# callbacks
# ------------------------------------------------------------------------------


@app.callback(
    [Output("data-storage", "data"),
     Output("dash-bio-viewer", "children"),
     Output("dropdown-data", "options"),
     Output("data-column-selector", "options"),
     Output("data-column-selector", "value"),
     Output("plot-data-selector", "options")],
    [Input("file-upload", "contents"),
     Input("file-upload", "filename"),
     Input('data-table', 'data_timestamp')],
    [State("data-storage", "data"),
     State("data-table", "data"),
     State("data-column-selector", "value"),
     State("dash-bio-viewer", "children")
     ]
)
def upload_data(content, filename, table_ts, stored_data, table_data,
                selected_columns, dbviewer):
    """
    Uploads the data from an xyz file and store them in the store component.
    Then set up the dropdowns, the table and the molecule viewer.
    """

    if table_ts is not None:
        # update stored data from current data in the table
        df = pd.DataFrame(stored_data)
        try:
            table_df = pd.DataFrame(table_data)
            table_df = table_df.astype({col: np.float for col in table_df
                                        if col != "species"})
            df.update(table_df)
        except ValueError:
            print("No update of data")

        all_data = df.to_dict("records")

    else:
        # Initial set up, read data from upload

        # read a default file
        # filename = app.get_asset_url("data/C28-D2.xyz")
        filename = "assets/data/C28-D2.xyz"
        mol = Molecule.from_file(filename)

        if content:
            content_type, content_str = content.split(",")
            _, ext = os.path.splitext(filename)
            decoded = base64.b64decode(content_str).decode("utf-8")
            # fdata = io.StringIO(decoded)
            try:
                mol = Molecule.from_str(decoded, fmt=ext[1:])
            except NameError:
                # TODO: Manage format error
                print("Unable to read format")

        # comute data
        ca = CurvatureAnalyzer(mol)

        # add a custom column for manual editing
        ca.data["custom"] = 0.0

        # all data for the store component
        all_data = ca.data.to_dict(orient="records")

        # Set the molecule 3D Viewer component
        dbviewer = dash_bio.Molecule3dViewer(
            id='molecule-viewer',
            backgroundColor="#FFFFFF",
            # backgroundOpacity='0',
            modelData=ca.get_molecular_data(),
            atomLabelsShown=True,
            selectionType='atom'
        )

        # options for the checklist in order to select the columns of the table
        selected_columns = ["atom_idx", "species", "angular_defect",
                            "pyrA", "n_star_A"]

    # options to select data mapped on atoms
    options = [{"label": name, "value": name} for name in ca.data
               if name not in ["atom_idx", "species", "atom_A", "star_A"]]
    options2 = [{"label": "histogram", "value": "histogram"}] + options

    # checklist options to select table columns
    tab_options = [{"label": name, "value": name} for name in ca.data]

    return all_data, dbviewer, options, tab_options, selected_columns, options2


@app.callback(
    [Output("data-table", "data"),
     Output("data-table", "columns")],
    [Input("data-storage", "modified_timestamp"),
     Input("data-column-selector", "value")],
    [State("data-storage", "data")]
)
def select_table_columns(ts, values, data):
    """
    Select columns displayed in the table. A custom column is available and 
    filled with zero by default.
    """

    # get data from the Store component
    df = pd.DataFrame(data)

    if values is None:
        # initial set up
        return [], []
    else:
        # fill the table with the selected columns
        tab_df = df[values]
        data = tab_df.to_dict("records")

        # add format
        columns = list()
        for column in tab_df:
            if column in {"atom_idx", "species", "neighbors", "custom"}:
                columns.append({"name": column, "id": column})
            elif column == "custom":
                columns.append(
                    {"name": column, "id": column, "editable": True})
            else:
                columns.append({
                    "name": column, "id": column, "type": "numeric",
                    "format": Format(
                        precision=4,
                        scheme=Scheme.fixed,
                    )
                })

        return data, columns


@app.callback(
    Output('molecule-viewer', 'styles'),
    [Input('dropdown-data', 'value'),
     Input('dropdown-colormap', "value"),
     Input("data-storage", "modified_timestamp"),
     Input("cm-min-value", "value"),
     Input("cm-max-value", "value"),
     Input("nan-color-value", "value")],
    [State("data-storage", "data")]
)
def map_data_on_atoms(selected_data, cm_name, ts, cm_min, cm_max, nan_color, data):
    """
    Map the selected data on the structure using a colormap.
    """

    df = pd.DataFrame(data)

    if selected_data:
        values = df[selected_data].values
        minval, maxval = np.nanmin(values), np.nanmax(values)

        # get cm boundaries values from inputs if they exist
        if cm_min is not None:
            minval = cm_min
        if cm_max is not None:
            maxval = cm_max

        # check nan_color value
        if nan_color is None or not HEX_COLOR_PATT.match(nan_color):
            nan_color = "#000000"

        normalize = mpl.colors.Normalize(minval, maxval)

        cm = plt.cm.get_cmap(cm_name)

        colors = list()
        for value in values:
            if np.isnan(value):
                colors.append(nan_color)
            else:
                colors.append(mpl.colors.rgb2hex(
                    cm(X=normalize(value), alpha=1)))

#        nan_idx = np.nonzero(np.isnan(values))[0]
#        norm_cm = cm(X=normalize(values), alpha=1)
#        colors = [mpl.colors.rgb2hex(color) for color in norm_cm]

        styles_data = {
            str(iat): {
                "color": colors[iat],
                "visualization_type": "stick"
            }
            for iat in range(len(df))
        }

    else:
        styles_data = {
            str(iat): {
                "color": ELEMENT_COLORS[df.species[iat]]
                if df.species[iat] in ELEMENT_COLORS else "#000000",
                "visualization_type": "stick"
            }
            for iat in range(len(df))
        }

    return styles_data


@app.callback(
    Output("colorbar", "figure"),
    [Input('dropdown-data', 'value'),
     Input('dropdown-colormap', 'value'),
     Input("data-storage", "modified_timestamp"),
     Input("cm-min-value", "value"),
     Input("cm-max-value", "value")],
    [State("data-storage", "data")]
)
def plot_colorbar(selected_data, cm_name, data_ts, cm_min, cm_max, data):
    """
    Display a colorbar according to the selected data mapped on to the structure.
    """

    if selected_data:
        # get data and boundaries
        values = pd.DataFrame(data)[selected_data].values
        minval, maxval = np.nanmin(values), np.nanmax(values)

        # get cm boundaries values from inputs if they exist
        if cm_min is not None:
            minval = cm_min

        if cm_max is not None:
            maxval = cm_max

        # set up fake data and compute corresponding colors
        npts = 100
        values = np.linspace(minval, maxval, npts)
        normalize = mpl.colors.Normalize(minval, maxval)

        cm = plt.cm.get_cmap(cm_name)
        cm_RGBA = cm(X=normalize(values), alpha=1) * 255
        cm_rgb = ["rgb(%d, %d, %d)" % (int(r), int(g), int(b))
                  for r, g, b, a in cm_RGBA]
        colors = [[x, c] for x, c in zip(np.linspace(0, 1, npts), cm_rgb)]

        trace = [
            go.Contour(
                z=[values, values],
                x0=values.min(),
                dx=(values.max() - values.min()) / (npts - 1),
                colorscale=colors,
                autocontour=False,
                showscale=False,
                contours=go.contour.Contours(coloring="heatmap"),
                line=go.contour.Line(width=0),
                hoverinfo="skip",
            ),
        ]
        figure = go.Figure(
            data=trace,
            layout=go.Layout(
                width=600, height=100,
                xaxis=dict(showgrid=False, title=selected_data),
                yaxis=dict(ticks="", showticklabels=False),
                margin=dict(t=0, b=0, l=0, r=0)
                # margin=dict(l=40, t=0, b=40, r=20, pad=0)
            )
        )
    else:
        figure = go.Figure(
            data=[],
            layout=go.Layout(
                width=600, height=100,
                xaxis=dict(ticks="", showticklabels=False, showgrid=False,
                           title=selected_data, zeroline=False),
                yaxis=dict(ticks="", showticklabels=False, showgrid=False,
                           title=selected_data, zeroline=False),
                margin=dict(l=0, t=0, b=0, r=0, pad=0)
            )
        )

    return figure


@app.callback(
    [Output("plot-data", "figure"),
     Output("plot-data-table", "data"),
     Output("plot-data-table", "columns")],
    [Input('dropdown-data', 'value'),
     Input("data-storage", "modified_timestamp"),
     Input("plot-data-selector", "value"),
     Input("nbins-slider", "value")],
    [State("data-storage", "data")]
)
def plot_data(selected_data1, data_ts, selected_data2, nbins, data):
    """
    Make a plot according to the data mapped on the structure. By default
    a histogram of the data is plotted. If another abscissa is chosen the
    data are plotted against this abscissa. 

    Statistical descriptors of these data are displayed on a table.
    """

    figure = go.Figure(data=[],
                       layout=go.Layout(template="plotly_white", height=600))
    tabdata = list()
    columns = list()

    if selected_data1:
        df = pd.DataFrame(data).dropna()

        # plot a histogram
        if selected_data2 == "histogram":
            figure = px.histogram(
                data_frame=df,
                x=selected_data1,
                histnorm="probability",
                marginal="box",
                nbins=nbins,
                height=600,
                color_discrete_sequence=["#2980b9"],
                title=selected_data1,
                template="plotly_white",
            )
            figure.layout.update(
                yaxis=dict(showgrid=False),
                xaxis=dict(showgrid=False),
            )

        # scatter plot with trend line
        else:
            figure = px.scatter(
                data_frame=df,
                x=selected_data2, y=selected_data1,
                symbol_sequence=["circle-open"],
                color_discrete_sequence=["#2980b9"],
                template="plotly_white",
            )
            figure.update_traces(marker=dict(size=10, line=dict(width=3)))

            # add a polynomial trend line to the plot
            xmin = np.nanmin(df[selected_data2])
            xmax = np.nanmax(df[selected_data2])
            xmin -= .05 * (xmax - xmin)
            xmax += .05 * (xmax - xmin)
            x = np.linspace(xmin, xmax, 100)
            p = np.poly1d(np.polyfit(
                df[selected_data2], df[selected_data1], deg=2))
            figure.add_trace(
                go.Scatter(
                    x=x, y=p(x),
                    mode="lines", showlegend=False,
                    line=dict(color="#2980b9", width=1),
                )
            )

        # set up table of plotted data with statistical descriptors
        if selected_data2 == "histogram":
            tabdata = df[selected_data1].describe().to_frame(
                name=selected_data1)
        else:
            tabdata = df[[selected_data1, selected_data2]].describe()

        tabdata = tabdata.transpose()
        tabdata.index.name = "data"
        tabdata.reset_index(inplace=True)

        columns = ["data", 'mean', 'std', 'min', '25%', '50%', '75%', 'max']
        tabdata = tabdata[columns].to_dict("records")
        fformat = Format(precision=4, scheme=Scheme.fixed)
        columns = [{"name": c, "id": c, "type": "numeric", "format": fformat}
                   for c in columns]

    return figure, tabdata, columns


if __name__ == '__main__':
    # app.run_server(debug=True, host='127.0.0.1')
    app.run_server(debug=False)
