
import base64
import os

from urllib.parse import quote as urlquote
from flask import Flask, send_from_directory

import dash
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html
from dash.dependencies import Input, Output, State

from jinja2 import Template
import data_handling 
import model_info

external_stylesheets = ['https://fonts.googleapis.com/css?family=Roboto:300,400,400i,700%7CRoboto+Mono&display=fallback']

# Normally, Dash creates its own Flask server internally. By creating our own,
# we can create a route for downloading files directly:
server = Flask(__name__)
app = dash.Dash(server=server)
# , external_stylesheets=external_stylesheets

@server.route("/download/<path:path>")
def download(path):
    """Serve a file from the upload directory."""
    return send_from_directory(UPLOAD_DIRECTORY, path, as_attachment=True)

app.layout = html.Div(
    [
      html.H2("Dataset:"),
      html.Div([
        dcc.Dropdown(id='ds-select',options=data_handling.get_ds_names(and_new=True),
                     clearable=False,placeholder="Select a dataset..."),
        # dbc.Button("...", id="ds-options",color="primary", className="mr-1",n_clicks=0)
          ],style={'width': '100%', 'display': 'inline-block'}
      ),
      html.Div([
        dcc.Input(id="ds-name", type="text", placeholder="enter dataset name..."),          
        ],style={'width': '100%', 'display': 'block'},id="ds-name-panel"
      ),
      html.Hr(),
      dcc.Upload(
          id="upload-data",
          children=html.Div(
              ["Drag and drop or click to select an image to upload."]
          ),
          style={"width": "100%","height": "40px","lineHeight": "40px",
                 "borderWidth": "1px","borderStyle": "dashed","borderRadius": "5px",
                 "textAlign": "center","margin": "10px",
          },
          multiple=True,
      ), 
      html.Div(id='ds-hidden', style={'display':'none'}),

      # List of images / data
      data_handling.initialise_table(), 

      html.Hr(),
      html.Div([
         dbc.Button("Configuration", id='btn-config', color="primary", outline=True,
                   className="mr-1",n_clicks=0),
         dbc.Button("Wrangle Data", id='btn-wrangle', color="primary", outline=True,
                   className="mr-1",n_clicks=0),
         dbc.Button("Run Pipeline", id='btn-run', color="success", 
                   className="mr-1",n_clicks=0),
        ],style={'width': '100%', 'display': 'block'}
      ),
      html.Div([
        html.H2("Configuration Settings"),
        # html.Div([
        #     dcc.Input(id="opts-px2um", type="number", disabled=True),
        #     html.Label(" µm / pixel",htmlFor="opts-px2um")
        # ]),
        dbc.Button(id="opts-xml",color="primary",active=True,n_clicks=0,block=True),
        dcc.Dropdown(id="opts-ver",options=model_info.get_menu(),
                 clearable=False,value=1),
         html.Div([
          dcc.Input(id="opts-m1", type="number", disabled=True), # tooltip="Segmentation pixel size (µm)"
          html.Label(" Segmentation pixel size (µm)",htmlFor="opts-m1")
        ]), 
        html.Div([
          dcc.Input(id="opts-m2", type="number", disabled=True), # tooltip="Segmentation patch size (pix)"
          html.Label(" Segmentation patch size (pix)",htmlFor="opts-m2")
        ]),
      ],style={'width': '100%', 'display': 'block'},id="config-panel")
    ],
    style={"max-width": "400px"},
)

# Show or hide the extra configuration panel (working)
@app.callback(
    [Output("btn-config", "active"),
     Output("config-panel", "style")],
     Input("btn-config", "n_clicks"), 
     State('btn-config','active'),
)
def toggle_configuration_display(n_clicks,active):

  if n_clicks == 0: 
    active = True
  if (active):
    return False,{'display':'none'}
  else:
    return True,{'display':'block'}

# Enable or disable dataset naming text input
@app.callback(
   [Output("ds-name","value"), 
    Output("ds-name-panel","style"),
    Output("ds-table","data")],
    Input("ds-select","value"),
    State("ds-select","options")   
)
def toggle_dataset_name(ds_name,ds_menu):
  if ds_name == ds_menu[-1]['value']:
    return "",{'display':'block'},data_handling.get_ds_contents(None)
  else:
    return ds_name,{'display':'none'},data_handling.get_ds_contents(ds_name)

# Toggle XML button in options menu
@app.callback(
    [Output("opts-xml", "active"),
     Output("opts-xml", "children")],
     Input("opts-xml", "n_clicks"), 
     State('opts-xml',"active"),
)
def toggle_xml(n_clicks,active):
  if n_clicks == 0: 
    active = False
  if (active):
    return False,"TIFF Output"
  else:
    return True,"XML Output"


# functions to run pipeline 
def make_config_files(ds_name,do_xml,mdl_ver,mdl_pps,mdl_ppy):
  if n_clicks == 0:
    return False

  with open(r"./conf/conf.yaml.template","r") as fi:
    conf = fi.read()

  with open(r"./conf/conf.yaml","w") as fi:
    fi.write(Template(conf).render( 
      name=ds_name,
      doXML=do_xml,
      version=mdl_ver,
      model_p2um=mdl_pps,
      model_npx=mdl_ppy ))

  with open(r"./conf/local.yaml.template","r") as fi:
    local = fi.read()

  with open(r"./conf/local.yaml","w") as fi:
    fi.write(Template(local).render( 
      installPath=os.path.dirname(os.path.realpath(__file__)),
      userDir=os.path.realpath(data_handling.UPLOAD_DIRECTORY) 
      ))



# Because several actions cause updates to the data table and dash doesn't
# allow a property to be set by more than one callback, we're stuck in this
# situation where there's one callback to rule them all, one callback to find
# them, one callback to bring them forth and in the darkness bind them

@app.callback(
   [Output("btn-run","active"),
    Output("ds-table", "data")],

    Input("btn-run","n_clicks"),     
    Input("btn-wrangle","n_clicks"), 
    Input("upload-data", "filename"), 
    Input("upload-data", "contents"),

    State("ds-table", "data")
    State("btn-run","active"),
    State("ds-name","value"),
    State("opts-xml","active"),
    State("opts-ver","value"),
    State("opts-m1", "value"),
    State("opts-m2","value"),
)
def one_true_callback( c_run, c_wrangle, file_up, data_up, # inputs
                       t_data, running, ds_name,        # state 
                       do_xml, mdl_ver, mdl_pps, mdl_ppy): # more state

  ctx = dash.callback_context
  if not ctx.triggered:
      return False, t_data
  else:
      which_input = ctx.triggered[0]['prop_id'].split('.')[0]

  print(which_input)

  # detect which input active

  if which_input is 'upload-data': 
    new_row = data_handling.upload_files(file_up,data_up,ds_name)
    t_data = data_handling.update_table(new_row, t_data, -1)
    return False, t_data

  elif which_input is 'btn-wrangle': 
    t_data = data_handling.wrangle_files(t_data, ds_name) # wrangle only
    return running, t_data

  elif which_input is 'btn-run': 
    t_data = data_handling.wrangle_files(t_data, ds_name) # wrangle
    make_config_files(ds_name,do_xml,mdl_ver,mdl_pps,mdl_ppy)
    return True, t_data

  return False, t_data

data_handling.add_callbacks(app)
model_info.add_callbacks(app)

if __name__ == "__main__":
    app.run_server(debug=True, port=8888)