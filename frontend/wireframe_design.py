#%%

import dash
import dash_core_components as dcc
import dash_html_components as html #elements of wireframe
import dash_bootstrap_components as dbc #grid
import plotly.graph_objs as go
import pandas as pd
import matplotlib.pyplot as plt

#%%

app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])

 
#%%
device_family_dropdown = dcc.Dropdown(
    options=[
        {'label': 'Device Family', 'value': 'Device Family'} #list comprehension for actual list
    ],
    placeholder="Select a Device Family"
    
)

device_dropdown = dcc.Dropdown(
    options=[
        {'label': 'Device', 'value': 'Device'}
    ],
    placeholder="Select a Device"
)

#%%

Elec_one_dropdown = dcc.Dropdown(
    options=[
        {'label': 'Electrode', 'value': 'Elect1'} #list comprehension for actual list
    ],
    placeholder="Select an Electrode"
    
) #is this upload

carrier_dropdown = dcc.Dropdown(
    options=[
        {'label': 'Carrier', 'value': 'Carr'} #list comprehension for actual list
    ],
    placeholder="Select an Carrier"    
)

simulation_type_dropdown = dcc.Dropdown(
    options=[
        {'label': 'Sim_type', 'value': 'Sim'} #list comprehension for actual list
    ],
    placeholder="Simuation Type"    
)

using_axons_dropdown = dcc.Dropdown(
    options=[
        {'label': 'using_axon', 'value': 'axon'} #list comprehension for actual list
    ],
    placeholder="Using Axon"    
)

axon_upload = dcc.Upload(
        id='upload_xml_axon',
        children=html.Div([
            html.A('Upload Axon XML file')
        ]),
        style={
            'width': '50%',
            'height': '40px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'dashed',
            'borderRadius': '5px',
            'textAlign': 'center',
            'margin': '10px'
        })

xml_upload = dcc.Upload(
        id='upload_xml',
        children=html.Div([
            html.A('Upload Nerve XML file')
        ]),
        style={
            'width': '25%',
            'height': '40px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'dashed',
            'borderRadius': '5px',
            'textAlign': 'center',
            'margin': '10px'
        })



email_results = dbc.InputGroup(
            [
                dbc.InputGroupAddon("Email results to : ", addon_type="prepend"),
                dbc.Input(placeholder="Type email address"),
            ],
            className="mb-3",
        )

run_button = dbc.Button("RUN", size="lg", className="mr-1")

fig= plt.figure(figsize = (10,10))
plt.title('LAYOUT PREVIEW', fontsize = 20)
layout_preview = dcc.Graph(figure=fig)
#%%   
layout = html.Div(
    [
        dbc.Row(
                [   
                    dbc.Col(html.Div([device_family_dropdown, device_dropdown]), width=4), \
                    dbc.Col(html.Div([simulation_type_dropdown]), width=4) 
                    
                ], justify = "between"
        ),        
        dbc.Row(
                [
                    dbc.Col(html.Div([Elec_one_dropdown]), width = 4), \
                    dbc.Col(html.Div([using_axons_dropdown, axon_upload]), width = 4)
                ],  justify = "between"            
        ),
        dbc.Row(
            [
                    dbc.Col(html.Div( dbc.Input(placeholder="X", type="text")), width={"size": 1, "order": 12} ) ,
                    dbc.Col(html.Div( dbc.Input(placeholder="Y", type="text")), width={"size": 1, "order": 12} ),
                    dbc.Col(html.Div( dbc.Input(placeholder="Z", type="text")), width={"size": 1, "order": 12} ),
            ]
        ),
        dbc.Row([
                    dbc.Col(html.Div([carrier_dropdown]), width=4), \
                    dbc.Col(html.Div([email_results, run_button]), width=4)
                    
                ], justify = "between"            
        ),
        dbc.Row(
            [
                dbc.Col(html.Div( dbc.Input(placeholder="L", type="text")), width={"size": 1, "order": 12} ) ,
                dbc.Col(html.Div( dbc.Input(placeholder="W", type="text")), width={"size": 1, "order": 12} ),
                dbc.Col(html.Div( dbc.Input(placeholder="D", type="text")), width={"size": 1, "order": 12} ),
            ]
        ),
        dbc.Row(                
                dbc.Col(html.Div([xml_upload]) , width = 4) )               
                        
    ])

#%%
app.layout = layout
if __name__ == '__main__':
    app.run_server(debug=True)