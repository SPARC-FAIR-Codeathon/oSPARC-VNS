#%%

import dash
import dash_core_components as dcc
import dash_html_components as html #elements of wireframe
import dash_bootstrap_components as dbc #grid
import plotly.graph_objs as go
#import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px

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
dcc.Upload( html.Div( [html.A('Upload Nerve XML file')] ), \
                       style = {'borderStyle': 'dashed',\
                                'borderWidth': '1px'})
       )

xml_upload = dcc.Upload(
        html.Div([
            html.A('Upload Nerve XML file')
        ])        
        )


email_results = dbc.InputGroup(
            [
                dbc.InputGroupAddon("Email", addon_type="prepend"),
                dbc.Input(placeholder="Type email"),
            ],
            className="mb-3",
        )

run_button = dbc.Button("RUN", size="lg", color = "danger", className="mr-1")


layout_preview = dcc.Graph(figure=px.scatter(width=480, height=400))


interactive_plot = dcc.Graph(figure=px.scatter(width=220, height=200))


plot_preview = dcc.Graph(figure=px.scatter(width=220, height=200))

#%%
#get card elements
#1st ROW
#card containing device family, device
First_column = html.Div([
                    dbc.Card(
                            dbc.CardBody([
                                    html.Div([
                                    device_family_dropdown, device_dropdown
                                    ]) 
                                ])
                            ),\
                    dbc.Card(
                            dbc.CardBody([
                                    html.Div([
                                    Elec_one_dropdown
                                    ]) 
                                ])
                            ),\
                    dbc.Card(
                            dbc.CardBody([
                                    html.Div([
                                    dbc.Row(
                                        [
                                           dbc.Col([dbc.Input(placeholder="X", type="text")]),\
                                           dbc.Col([dbc.Input(placeholder="Y", type="text")]),\
                                           dbc.Col([dbc.Input(placeholder="Z", type="text")])
                                        ])
                        ])
                    ])
                    ),\

                   dbc.Card(
                            dbc.CardBody([
                                    html.Div([
                                    dbc.Row(
                                        [
                                           dbc.Col([dbc.Input(placeholder="L", type="text")]),\
                                           dbc.Col([dbc.Input(placeholder="W", type="text")]),\
                                           dbc.Col([dbc.Input(placeholder="D", type="text")])
                                        ])
                        ])
                    ])
                    ),\

                   dbc.Card(
                            dbc.CardBody([
                                    html.Div([
                                    carrier_dropdown
                                    ]) 
                                ])
                            ),\

                   dbc.Card(
                            dbc.CardBody([
                                    html.Div([
                                    dbc.Row(
                                        [
                                           dbc.Col([dbc.Input(placeholder="L", type="text")]),\
                                           dbc.Col([dbc.Input(placeholder="W", type="text")]),\
                                           dbc.Col([dbc.Input(placeholder="D", type="text")])
                                        ])
                        ])
                    ])
                    ),\

                   dbc.Card(
                            dbc.CardBody([
                                    html.Div([
                                    dcc.Upload( html.Div( [html.A('Upload Nerve XML file')] ), style = {'borderStyle': 'dashed',\
                                                                                                        'borderWidth': '1px'})
                                    ]) 
                                ])
                            )  

])

#%%
#LAST COLUMN
Last_column = html.Div([
                    dbc.Card(
                            dbc.CardBody([
                                    html.Div([
                                    simulation_type_dropdown
                                    ]) 
                                ])
                            ),\

                    dbc.Card(
                            dbc.CardBody([
                                    html.Div([
                                    using_axons_dropdown, axon_upload
                                    ]) 
                                ])
                            ), \

                    dbc.Card(
                            dbc.CardBody([
                                    html.Div([
                                    email_results
                                    ]) 
                                ])
                            ), \

                    html.Br(), \
                    html.Br(), \
                    
                    html.Div([run_button]) 
                                

    


])


#%%
layout_preview = html.Div([
                    dbc.Card(
                            dbc.CardBody([
                                    html.Div([
                                    layout_preview
                                    ]) 
                                ])
                            ),
                        ])

int_preview = html.Div([
                    dbc.Card(
                            dbc.CardBody([
                                    html.Div([
                                    interactive_plot
                                    ]) 
                                ])
                            ),
                        ])

plot_preview = html.Div([
                    dbc.Card(
                            dbc.CardBody([
                                    html.Div([
                                    plot_preview
                                    ]) 
                                ])
                            ),
                        ])



#%%



#%%   
layout = html.Div([
            dbc.Card(
            dbc.CardBody([
    
                dbc.Row(
                        [   
                            dbc.Col([First_column], width= 3), \
                            dbc.Col(dbc.Card(
                                    dbc.CardBody(
                                            dbc.Col(
                                                    [dbc.Row([
                                                            layout_preview,                                         
                                                            dbc.Col([int_preview], width = "auto"), 
                                                            dbc.Col([plot_preview], width = "auto")
                                                            ])], width = "auto"                              
                                                   ))
                                             ), width = 6),\
                            dbc.Col([Last_column], width= 3)
                            
                        ], justify = "between", form = True, 
                    no_gutters=True)
                ])
            )
        ])

#%%
app.layout = layout
if __name__ == '__main__':
    app.run_server(debug=True)
