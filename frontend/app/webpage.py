#%%

from flask import Flask, send_from_directory

import dash
import dash_core_components as dcc
import dash_html_components as html #elements of wireframe
import dash_bootstrap_components as dbc #grid
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
                     placeholder="Select a Device Family"), 
        dcc.Dropdown(id="device-dropdown", options = callbacks.list_devices(None), value = None, 
                     placeholder="Select a Device", disabled=False),
        dcc.Upload(id="upload-device", className="uploada",
                   children=html.Div(
                     ["Drag and drop or click to",html.Br(),"select an array design to upload"]
                    )),
        dbc.Button("Save", id='btn-save-array', color="primary", outline=True,
                       className="mr-1",n_clicks=0),         
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
          dbc.Col([dbc.Input(id="elec-x",placeholder="X", type="number", min=0)]),
          dbc.Col([dbc.Input(id="elec-y",placeholder="Y", type="number", min=0)]),
          dbc.Col([dbc.Input(id="elec-z",placeholder="Z", type="number", min=0)]) ]),
        dbc.Row([
          dbc.Col([dbc.Input(id="elec-w",placeholder="W", type="number", min=0)]), # "width" is parallel to nerve
          dbc.Col([dbc.Input(id="elec-h",placeholder="H", type="number", min=0)]), # "height" is perpendicular to nerve
          dbc.Col([dbc.Input(id="elec-d",placeholder="D", type="number")       ]) ])
        ])  # CardBody: Electrode Configuration
      ]),
    dbc.Card([
      dbc.CardHeader("Electrode Carrier"),
      dbc.CardBody([
        dbc.Row([
          dbc.Col([dbc.Input(id="outer-x",placeholder="X", type="number", min=0)]),
          dbc.Col([dbc.Input(id="outer-y",placeholder="Y", type="number", min=0)]),
          dbc.Col([dbc.Input(id="outer-z",placeholder="L", type="number", min=0)]) ])
        ])  # CardBody: Electrode Carrier
    ])], width=3), # Left column
  #%%
  dbc.Col([
      html.Img(src=graphics.encode(graphics.array_SVG(None,None)),id="view-upper"), 
      html.Img(src=graphics.encode(graphics.nerve_SVG(None,None)),id="view-lower")
    ], id="viewport", width=6, style={'align-items':'center'}), # middle 
  #%%
  dbc.Col([
    dbc.Card([
      dbc.CardHeader("Select Nerve Anatomy"),
      dbc.CardBody([
        dcc.Upload(id="upload-nerve", className="uploada",
                   children=html.Div(["Drag and drop or click to",html.Br(),"select a nerve anatomy (MBF-XML) file"])),
        html.Div([dcc.Slider(id="nerve-x",min=-1000,max=1000,step=10,value=0)],style={'height':'30px'}),
        html.Div([dcc.Slider(id="nerve-y",min=-1000,max=1000,step=10,value=0)],style={'height':'30px'}),
        html.Div([dcc.Slider(id="nerve-r",min=-360,max=360,step=5,value=0)],style={'height':'30px'}),
        dbc.Button("Save Position", id='btn-save-nerve', color="primary", outline=True,
                       className="mr-1",n_clicks=0),
        dcc.Dropdown(id="nerve-dropdown", 
                     options = callbacks.list_nerveClasses(), value = callbacks.get_default_nerveClass(app))
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
    dcc.Store(id="device-json", storage_type='local')
  ]) # layout

callbacks.add_callbacks(app)
graphics.add_callbacks(app)
# middle_layer.add_callbacks(app)

#%%
if __name__ == '__main__':
    app.run_server(debug=True)

