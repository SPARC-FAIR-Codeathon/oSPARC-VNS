#%%

from flask import Flask, send_from_directory

import dash
import dash_core_components as dcc
import dash_html_components as html #elements of wireframe
import dash_bootstrap_components as dbc #grid
import dash_extensions as dec
import plotly.graph_objs as go
#import pandas as pd

import callbacks
import graphics

#%%
# Normally, Dash creates its own Flask server internally. By creating our own,
# we can create a route for downloading files directly:
server = Flask(__name__)
app = dash.Dash(server=server, external_stylesheets=[dbc.themes.BOOTSTRAP])
# ,             external_scripts=["//daybrush.com/moveable/release/latest/dist/moveable.min.js"])


#%%
# Define path so that we can serve downloads 
@server.route("/download/<path:path>")
def download(path):
    """Serve a file from the upload directory."""
    return send_from_directory(UPLOAD_DIRECTORY, path, as_attachment=True)


#%% 

# I've reorganised this into a condensed form so we can see more of it in a single screen. Excessive whitespace is an anti-pattern.

# first_col is the array device controls 
app.layout = dbc.Row([
  dbc.Col([
    dbc.Card([ 
      dbc.CardHeader("Select device"),
      dbc.CardBody([
        dcc.Dropdown(id="device-family-dropdown", 
                     options = callbacks.list_deviceFamilies(), value = callbacks.get_default_dfamily(app), 
                     placeholder="Select a Device Family", persistence=True), 
        dcc.Dropdown(id="device-dropdown", options = callbacks.list_devices(None), value = None, 
                     placeholder="Select a Device", disabled=False, persistence=True),
        dcc.Upload(id="upload-device", className="uploada",
                   children=html.Div(
                     ["Drag and drop or click to",html.Br(),"select an array design to upload"]
                    )), # dcc.upload
        dbc.Button("Save Device Configuration", id='btn-save-array', color="primary", outline=True,
                       className="mr-1",n_clicks=0), 
        dec.Download(id="download-file") # invisible component 
        ]) # CardBody: device
      ]), 
    dbc.Card([
      dbc.CardHeader("Electrode Configuration"),
      dbc.CardBody([
        dcc.Dropdown(id="elec-dropdown", options = [{'label':'E1','value':1}], 
                     placeholder="Edit Electrode(s)", value = None, multi=True), 
        html.Div([ dbc.Col([
          dbc.Button("Add Electrode",id="add-elec", color="secondary", outline=True, className="mr-1",
                                n_clicks=0, size="sm", style={'width':'45%'}), 
          dbc.Button("Remove",id="rem-elec", color="secondary", outline=True, className="mr-1",
                                n_clicks=0, size="sm", style={'width':'45%'})], style={'padding':'4px'}) 
        ], style={'align-items':'center'}), 
        dbc.Row([
          dbc.Col([dbc.Input(id="elec-x",placeholder="X", type="number", debounce=True, step="Any")]),
          dbc.Col([dbc.Input(id="elec-y",placeholder="Y", type="number", debounce=True, step="Any")]),
          dbc.Col([dbc.Input(id="elec-z",placeholder="Z", type="number", debounce=True, step="Any")]) ]),
        dbc.Row([
          dbc.Col([dbc.Input(id="elec-w",placeholder="W", type="number", debounce=True, step="Any")]), # "width" is parallel to nerve
          dbc.Col([dbc.Input(id="elec-h",placeholder="H", type="number", debounce=True, step="Any")]), # "height" is perpendicular to nerve
          dbc.Col([dbc.Input(id="elec-d",placeholder="D", type="number", debounce=True, step="Any")       ]) ])
        ])  # CardBody: Electrode Configuration
      ]),
    dbc.Card([
      dbc.CardHeader("Electrode Carrier"),
      dbc.CardBody([
        dbc.Row([
          dbc.Col([dbc.Input(id="outer-x",placeholder="X", type="number", debounce=True, step="Any")]),
          dbc.Col([dbc.Input(id="outer-y",placeholder="Y", type="number", debounce=True, step="Any")]),
          dbc.Col([dbc.Input(id="outer-z",placeholder="L", type="number", debounce=True, step="Any")]) ])
        ])  # CardBody: Electrode Carrier
    ])], width=3), # Left column
  #%%
  dbc.Col( 
    html.Div([
      html.Img(src=graphics.encode(graphics.array_SVG(None,None)),id="view-device"), 
    ],style={'margin':'auto','display':'block'}), id="viewport-col", width=6), # middle 
  #%%
  dbc.Col([
    dbc.Card([
      dbc.CardHeader("Select Nerve Anatomy"),
      dbc.CardBody([
        html.Div(["nerve.xml"],id="nerve-name",style={'display':'none'}),
        dcc.Upload(id="upload-nerve", className="uploada",
                   children=html.Div(["Drag and drop or click to",html.Br(),"select a nerve anatomy (MBF-XML) file"],id="upload-nerve-text")),
        html.Div([html.Div([
          dcc.Slider(id="nerve-x",min=-1000,max=1000,step=10,value=0,updatemode='drag')],style={"margin-top":"9px"}),
          dbc.InputGroupAddon("0 µm",addon_type="append",id='nerve-x-lbl',style={"height":"30px"})],
          style={'height':'34px',"display":"grid","grid-template-columns": "85% 15%"}),
        html.Div([html.Div([
          dcc.Slider(id="nerve-y",min=-1000,max=1000,step=10,value=0,updatemode='drag')],style={"margin-top":"9px"}),
          dbc.InputGroupAddon("0 µm",addon_type="append",id='nerve-y-lbl',style={"height":"30px"})],
          style={'height':'34px',"display":"grid","grid-template-columns": "85% 15%"}),
        html.Div([html.Div([
          dcc.Slider(id="nerve-r",min=-360,max=360,step=5,value=0,updatemode='drag')],style={"margin-top":"9px"}),
          dbc.InputGroupAddon("0°",addon_type="append",id='nerve-r-lbl',style={"height":"30px"})],
          style={'height':'34px',"display":"grid","grid-template-columns": "85% 15%"}),        
        dcc.Dropdown(id="nerve-dropdown", options = callbacks.list_nerveClasses(),persistence=True),
        dbc.Button("Save Nerve Configuration", id='btn-save-nerve', color="primary", outline=True,
                       className="mr-1",n_clicks=0),
      ]) # CardBody: Nerve Configuration
    ]), 
    dbc.Card([
      dbc.CardHeader("Run Controls"),
      dbc.CardBody([
        dcc.Dropdown(id="run-mode-dropdown", 
                     options = [{'label':'Nerve Recording','value':'full'}, 
                                {'label':'Fields only (faster)','value':'fast'}], 
                     value = callbacks.get_default_runMode(app), 
                     clearable = False), 
        html.Div([
            html.Div("simulated spike-rate (imp/s/axon):"),
            dbc.Input(id="spike-rate",placeholder="flat[0.2,2]", type="text")],
            id="div-spike-rate",style={"display":"block"}),
        html.Div("If you would like your results emailed to you, please enter your email below:"),
        dbc.InputGroup([
          dbc.InputGroupAddon("Email", addon_type="prepend"),
          dbc.Input(id="user-email",value=callbacks.get_default_email(app),type="email",placeholder="(optional)"),
            ],className="mb-3"),
        dbc.Button("Run Model", id="btn-run", size="lg", color = "success", className="mr-1")
     ]) # CardBody: run control
    ])], width=3), # Right column
    dcc.Store(id="device-json", storage_type='session'),
    dcc.Store(id="nerve-json", storage_type='session')
  ]) # layout

callbacks.add_callbacks(app)
graphics.add_callbacks(app)
# middle_layer.add_callbacks(app)

#%%
if __name__ == '__main__':
    app.run_server(debug=True)

