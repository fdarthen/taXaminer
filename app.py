import dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output
import plotly.express as px
import pandas as pd
import numpy as np
import sys
import yaml

config_path = sys.argv[1]

# read parameters from config file
config_obj=yaml.safe_load(open(config_path,'r'))
output_path=config_obj.get('output_path') # complete path to output directory

# reading data and creating preliminary python version of 3D plot
data = pd.read_csv(output_path + "taxonomic_assignment/gene_table_taxon_assignment.csv")
fig = px.scatter_3d(data, x='Dim.1', y='Dim.2', z='Dim.3',
              color='plot_label')
fig.update_traces(marker=dict(size=6))


app = dash.Dash(__name__)
app.title= "MILTS"

app.layout = html.Div(
    children=[
        html.H1(children="MILTS",
                className="header-title",
),
        dcc.Graph(
                id="scatter3d",
                config={"displayModeBar": True},
                animate=True,
                figure=fig,
                className="plot",
        ),
    ]
)


if __name__ == "__main__":
    app.run_server(host='127.0.0.1', port='8050', debug=True)
