
import dash
import dash_core_components as dcc
import dash_html_components as html #elements of wireframe
import dash_bootstrap_components as dbc #grid
import dash_extensions as dec

from urllib.parse import quote as urlquote
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from glob import glob

import base64
import json
import math
import os
import re

import user_files
import osparc_api


# some utilities 

def file_download_link(filename):
    """Create a Plotly Dash "A" element that downloads a file from the app."""
    location = "/download/{}".format(urlquote(filename))
    return html.A(filename, href=location)

def which_input(ctx = dash.callback_context):    
    if not ctx.triggered: raise PreventUpdate    
    else: return ctx.triggered[0]["prop_id"].split(".")[0]

def save_file(filepath,content):
    """Decode and store a file uploaded with Plotly Dash."""
    try: 
      data = content.encode("utf8").split(b";base64,")[1]
    except:
      print(content.encode("utf8").split(b";base64,"))
      raise PreventUpdate           
    base = os.path.dirname(filepath)
    if not os.path.exists(base): os.makedirs(base)
    with open(filepath, "wb") as fp:
        fp.write(base64.decodebytes(data))


#%% 
#
# Emulate database functionality 
#
#%%

# default value (user/session-based)
def db_query(tag,table="USER"):
    print("TODO: query column {} from table {}".format(tag,table))

# get user ID 
def get_user_ID(app=None):
    return 1

# get session ID 
def get_session_ID(app=None):
    return 1

def get_default_session(app=None):
    return 1


# what is user"s email?
def get_default_email(app):
    db_query("email","USER")
    return None



#%% 
# 
# Callbacks for interacting with json representation of array
#
#%%

def update_electrodePositions(mode,this,e_select,ex=None,ey=None,ez=None):

  try:
    eti = this["array"]["ElectrodeTypeIndex"]
  except:
    print("\nIssue getting ElectrodeTypeIndex:")
    print(this["array"])
    print("\n")

  if mode == "WHAT": return eti
  if mode == "PUT": # update "this" from ex ey ez
    if ex is None or ey is None or ez is None: raise PreventUpdate
    if not e_select: raise PreventUpdate

    s = e_select[0]
    try:
      this["array"]["ElectrodePositions"][s] = [ex,ey,ez]
    except:
      print("\nIssue setting ElectrodePositions:")
      print(this["array"]["ElectrodePositions"])
      print(e_select)
      print(s)
      print("\n")
      raise PreventUpdate

    return this
  if mode == "GET": # update ex ey ez from "this"

    print("GET elec-xyz:")
    print(e_select)

    if e_select is not None and len(e_select) == 1:
      s = e_select[0]
      ex = this["array"]["ElectrodePositions"][s][0]
      ey = this["array"]["ElectrodePositions"][s][1]
      ez = this["array"]["ElectrodePositions"][s][2]
    else:
      ex = None
      ey = None
      ez = None

    return ex,ey,ez

def update_electrodeDimension(mode,this,e_select,ew=None,eh=None,ed=None):

  these = lambda v : [v[u] for u in range(0,len(v)) if u in e_select ]
  those = lambda v : [v[u] for u in range(0,len(v)) if u not in e_select ]

  try:
    eti = this["array"]["ElectrodeTypeIndex"]
  except:
    print("\nIssue getting ElectrodeTypeIndex:")
    print(this["array"])
    print("\n")

  single_etype = ( min(eti) == 1 and max(eti) == 1 )
  single_edata = ( not isinstance(this["array"]["ElectrodeDimensions"][0],list) )

  if mode == "WHAT": return single_etype, single_edata
  if mode == "PUT": 
    if ew is None or eh is None or ed is None: raise PreventUpdate
    if not e_select: 
      if single_etype and single_edata:
        this["array"]["ElectrodeDimensions"][0] = ew
        this["array"]["ElectrodeDimensions"][2] = eh
        this["array"]["InsetDepth"] = ed
      elif single_etype:
        this["array"]["ElectrodeDimensions"][0][0] = ew
        this["array"]["ElectrodeDimensions"][0][2] = eh
        this["array"]["InsetDepth"] = ed
      else: raise PreventUpdate  
      return this,ew,eh,ed
    else: # one or more electrodes selected 
      print("PUT elec-dim:")
      print(e_select)

      if min(these(eti)) == max(these(eti)): # "these" all same ETI
        print("these all same type-index")
        new_eti = these(eti)[0] # proposed index 
        print("ETI not in selection:")
        print(those(eti))
        if new_eti in those(eti): # some other things also have the proposed ID
          print("looking for a new ETI...")
          for n in range(1,max(eti)+2): 
            if n not in those(eti): # find lowest ID not already used
              print("found a new ETI: {}".format(n))
              new_eti = n 
              break
        print("update type index {}".format(new_eti))
      else: # "these" do not all have the same index 
        print("mismatch of type-indices")
        for n in range(1,max(eti)+2): 
          if n not in those(eti): # find lowest ID not already used
            new_eti = n 

      if single_edata: # make into a list
        print("ElectrodeDimensions becomes list")
        this["array"]["ElectrodeDimensions"] = [this["array"]["ElectrodeDimensions"]]

      ek = this["array"]["ElectrodeDimensions"][0][1]
      if new_eti >= len(this["array"]["ElectrodeDimensions"]):          
          this["array"]["ElectrodeDimensions"].append([ew,ek,eh])
          new_eti = len(this["array"]["ElectrodeDimensions"])
          print("ElectrodeDimensions append")
      else: 
          this["array"]["ElectrodeDimensions"][new_eti] = [ew,ek,eh]
          print("ElectrodeDimensions assign")

      for e in e_select: # assign ETI
        this["array"]["ElectrodeTypeIndex"][e] = new_eti
        print("ElectrodeTypeIndex assign [{}] {}".format(e,new_eti))
    # PUT suceeded in updating "this" from ew eh ed
    return this,ew,eh,ed

  if mode == "GET":
    print("GET elec-dim:")
    print(e_select)

    ed = this["array"]["InsetDepth"]
    if single_edata:
      ew = this["array"]["ElectrodeDimensions"][0]
      eh = this["array"]["ElectrodeDimensions"][2]
    elif single_etype:
      ew = this["array"]["ElectrodeDimensions"][0][0]
      eh = this["array"]["ElectrodeDimensions"][0][2]
    elif not e_select: 
      ew = None
      eh = None
    else:
      s = e_select[0]
      ew = this["array"]["ElectrodeDimensions"][s][0]
      eh = this["array"]["ElectrodeDimensions"][s][2]
    # GET suceeded in updating ew eh ed from "this"
    return ew,eh,ed

def update_device_carrier(mode,this,cx=None,cy=None,cz=None):

  if "c_thickness" in this["array"]["carrier"]: 
        device_type = "flat"
  elif "cuff_IDr" in this["array"]["carrier"] and \
                     this["array"]["carrier"]["cuff_IDr"] > 0.45 * \
                     this["array"]["carrier"]["cuff_IDx"] :
        device_type = "cyl"
  else: device_type = "rect"

  if mode == "WHAT": return device_type
  if mode == "PUT":
    if device_type == "flat": 
      this["array"]["carrier"]["c_len"] = cx
      this["array"]["carrier"]["c_wid"] = cy
      this["array"]["carrier"]["c_thickness"] = cz
    elif device_type == "cyl":
      this["array"]["carrier"]["cuff_IDx"] = cx
      this["array"]["carrier"]["cuff_IDy"] = cx
      this["array"]["carrier"]["cuff_IDr"] = cx * 0.99/2
      this["array"]["carrier"]["cuff_length"] = cz
    else:
      this["array"]["carrier"]["cuff_IDx"] = cx
      this["array"]["carrier"]["cuff_IDy"] = cy
      this["array"]["carrier"]["cuff_length"] = cz

    return this,device_type

  if mode == "GET":
    try:
      if device_type == "flat": 
        cx = this["array"]["carrier"]["c_len"]
        cy = this["array"]["carrier"]["c_wid"]
        cz = this["array"]["carrier"]["c_thickness"]
      else:      
        cx = this["array"]["carrier"]["cuff_IDx"]
        cy = this["array"]["carrier"]["cuff_IDy"]
        cz = this["array"]["carrier"]["cuff_length"]
    except:
        print("\nPossibly bad device data:")
        print(this["array"]["carrier"])
        print("\n")        
    
    return cx,cy,cz,device_type


#%% Layout


# I reorganised this into a condensed form so we can see more of it in a single screen. Excessive whitespace is an anti-pattern.

layout = dbc.Row([ # model set-up controls 
    dbc.Col([ # first_col is the array device controls 
      dbc.Card([ 
        dbc.CardHeader("Select device"),
        dbc.CardBody([
          dcc.Dropdown(id="device-family-dropdown", 
                       options = user_files.list_deviceFamilies(), 
                       placeholder="Select a Device Family", persistence=True), 
          dcc.Dropdown(id="device-dropdown", options = user_files.list_devices(None), value = None, 
                       placeholder="Select a Device", disabled=False, persistence=True),
          dcc.Upload(id="upload-device", className="uploada",
                     children=html.Div(
                       ["Drag and drop or click to",html.Br(),"select an array design to upload"]
                      )), # dcc.upload
          dbc.Button("Save Device Configuration", id="btn-save-array", color="primary", outline=True,
                         className="mr-1",n_clicks=0), 
          dec.Download(id="download-device") # invisible component 
          ]) # CardBody: device
        ]), 
      dbc.Card([
        dbc.CardHeader("Electrode Configuration"),
        dbc.CardBody([
          dcc.Dropdown(id="elec-dropdown", options = [{"label":"E1","value":1}], 
                       placeholder="Edit Electrode(s)", value = None, multi=True), 
          html.Div([ dbc.Col([
            dbc.Button("Add Electrode",id="add-elec", color="secondary", outline=True, className="mr-1",
                                  n_clicks=0, size="sm", style={"width":"45%"}), 
            dbc.Button("Remove",id="rem-elec", color="secondary", outline=True, className="mr-1",
                                  n_clicks=0, size="sm", style={"width":"45%"})], style={"padding":"4px"}) 
          ], style={"align-items":"center"}), 
          dbc.Row([
            dbc.Col([dbc.Input(id="elec-x",placeholder="X", type="number", debounce=True, step="Any")]),
            dbc.Col([dbc.Input(id="elec-y",placeholder="Y", type="number", debounce=True, step="Any")]),
            dbc.Col([dbc.Input(id="elec-z",placeholder="Z", type="number", debounce=True, step="Any")])],form=True),
          dbc.Row([
            dbc.Col([dbc.Input(id="elec-w",placeholder="W", type="number", debounce=True, step="Any")]), # "width" is parallel to nerve
            dbc.Col([dbc.Input(id="elec-h",placeholder="H", type="number", debounce=True, step="Any")]), # "height" is perpendicular to nerve
            dbc.Col([dbc.Input(id="elec-d",placeholder="D", type="number", debounce=True, step="Any")])],form=True)
          ])  # CardBody: Electrode Configuration
        ]),
      dbc.Card([
        dbc.CardHeader("Electrode Carrier"),
        dbc.CardBody([
          dbc.Row([
            dbc.Col([dbc.Input(id="outer-x",placeholder="X", type="number", debounce=True, step="Any")]),
            dbc.Col([dbc.Input(id="outer-y",placeholder="Y", type="number", debounce=True, step="Any")]),
            dbc.Col([dbc.Input(id="outer-z",placeholder="L", type="number", debounce=True, step="Any")])],form=True)
          ])  # CardBody: Electrode Carrier
      ])], width=3), # Left column
    #%%
    dbc.Col( 
      html.Div([
        html.Img(src="",id="view-device"), 
      ],style={"margin":"auto","display":"block","text-align":"center"}), id="viewport-col", width=6), # middle 
    #%%
    dbc.Col([
      dbc.Card([
        dbc.CardHeader("Select Nerve Anatomy"),
        dbc.CardBody([
          html.Div(["nerve.xml"],id="nerve-name",style={"display":"none"}),
          dcc.Upload(id="upload-nerve", className="uploada",
                     children=html.Div(["Drag and drop or click to",html.Br(),"select a nerve anatomy (MBF-XML) file"],id="upload-nerve-text")),
          html.Div([html.Div([
            dcc.Slider(id="nerve-x",min=-1000,max=1000,step=10,value=0,updatemode="drag")],style={"margin-top":"9px"}),
            dbc.InputGroupAddon("0 µm",addon_type="append",id="nerve-x-lbl",style={"height":"30px"})],
            style={"height":"34px","display":"grid","grid-template-columns": "85% 15%"}),
          html.Div([html.Div([
            dcc.Slider(id="nerve-y",min=-1000,max=1000,step=10,value=0,updatemode="drag")],style={"margin-top":"9px"}),
            dbc.InputGroupAddon("0 µm",addon_type="append",id="nerve-y-lbl",style={"height":"30px"})],
            style={"height":"34px","display":"grid","grid-template-columns": "85% 15%"}),
          html.Div([html.Div([
            dcc.Slider(id="nerve-r",min=-180,max=180,step=5,value=0,updatemode="drag")],style={"margin-top":"9px"}),
            dbc.InputGroupAddon("0°",addon_type="append",id="nerve-r-lbl",style={"height":"30px"})],
            style={"height":"34px","display":"grid","grid-template-columns": "85% 15%"}),        
          dcc.Dropdown(id="axon-pop-dropdown", options = user_files.list_nerveClasses(),persistence=True),
          dbc.Button("Save Nerve Configuration", id="btn-save-nerve", color="primary", outline=True,
                         className="mr-1",n_clicks=0),
          dec.Download(id="download-nerve") # invisible component         
        ]) # CardBody: Nerve Configuration
      ]), 
      dbc.Card([
        dbc.CardHeader("Run Controls"),
        dbc.CardBody([
          dcc.Dropdown(id="run-mode-dropdown", 
                       options = [{"label":"Nerve Recording","value":"full"}, 
                                  {"label":"Fields only (faster)","value":"fast"}],
                       clearable = False, persistence=True), 
          html.Div([
              html.Div("simulated spike-rate (imp/s/axon):"),
              dbc.Input(id="spike-rate",placeholder="flat[0.2,2]", type="text")],
              id="div-spike-rate",style={"display":"block"}),
          html.Div("If you would like your results emailed to you, please enter your email below:"),
          dbc.InputGroup([
            dbc.InputGroupAddon("Email", addon_type="prepend"),
            dbc.Input(id="user-email",type="email",placeholder="(optional)"),
              ],className="mb-3"),
          dbc.Button("Run Model", id="btn-run", size="lg", color = "success", className="mr-1")
       ]) # CardBody: run control
      ])], width=3), # Right column
    ]) # layout SETUP page


#%% 
#
# DEFINE CALLBACKS
# 
#%%
def add_callbacks(app):

  # Callback to show/hide "firing rate" field
  @app.callback( Output("div-spike-rate","style"), Input("run-mode-dropdown","value") )
  def toggle_show_spikerate(run_mode):        
    if run_mode == "full":  return {"display":"block"}
    else:                   return {"display":"none"}

  # Callback to set up device dropdown "firing rate" field
  @app.callback( [Output("device-dropdown","options"), 
                  Output("device-dropdown","value"), 
                  Output("device-dropdown","disabled")], 
                  Input("device-family-dropdown","value"),
                  State("navbar-session","value") )
  def update_device_dropdown(family,session):
    opts = user_files.list_devices(family)
    default = [{"label":"...","value":""}]
    if opts is None: return default, None, True
    if not opts:     return default, None, True
    if family == "user" and os.path.exists("../data/u/{}/{}/array.json".format(get_user_ID(),session)):
                     return opts, opts[session-1]["value"], False
    return opts, opts[0]["value"], False

    # TODO here: if url changed refresh and trigger refresh cascade


  # MAJOR callback: update device representation based on user IO 
  @app.callback(Output("device-json","data"),
                [Input("device-dropdown","value"),
                 Input("add-elec","n_clicks"),
                 Input("rem-elec","n_clicks"),
                 Input("elec-x","value"),
                 Input("elec-y","value"),
                 Input("elec-z","value"),
                 Input("elec-w","value"),
                 Input("elec-h","value"),
                 Input("elec-d","value"),
                 Input("outer-x","value"),
                 Input("outer-y","value"),
                 Input("outer-z","value"),
                 State("device-json","data"),
                 State("device-family-dropdown","value"),
                 State("elec-dropdown","value")] )
  def update_device(device,ua,ur,ex,ey,ez,ew,eh,ed,cx,cy,cz,
                    this, family, e_select):

    was_loaded = False
    if device is None: raise PreventUpdate
    if this is None: 
      this = user_files.get_device_json(device,family)
      was_loaded = True

    if isinstance(this,list) and this: this = this[0]
    clicked = which_input()

    print("update_device: "+clicked)

    if clicked=="device-dropdown" and not was_loaded: # got a new device?
        this = user_files.get_device_json(device,family)
        was_loaded = True
        return this

    these = lambda v : [v[u] for u in range(0,len(v)) if u in e_select ]
    those = lambda v : [v[u] for u in range(0,len(v)) if u not in e_select ]

    if clicked=="add-elec": 
      if e_select:
            eti = this["array"]["ElectrodeTypeIndex"][e_select[0]]
      else: eti = this["array"]["ElectrodeTypeIndex"][-1]

      this["array"]["ElectrodePositions"].append([0,0,0])
      this["array"]["ElectrodeTypeIndex"].append(eti)

      if "ElectrodeAngle" in this["array"]: 

        if e_select:
              eti = this["array"]["ElectrodeTypeIndex"][e_select[0]]
        else: eti = this["array"]["ElectrodeTypeIndex"][-1]
        this["array"]["ElectrodeAngle"].append(this["array"]["ElectrodeAngle"][eti])
        return this

    elif clicked=="rem-elec": 
      if not e_select: raise PreventUpdate

      this["array"]["ElectrodePositions"] = those(this["array"]["ElectrodePositions"])
      this["array"]["ElectrodeTypeIndex"] = those(this["array"]["ElectrodePositions"])
      if "ElectrodeAngle" in this["array"]: 
        this["array"]["ElectrodeAngle"] = those(this["array"]["ElectrodeAngle"])

      return this

    elif clicked=="elec-x" or clicked=="elec-z" or clicked=="elec-z":
      if not e_select: raise PreventUpdate
      if ex is None or ey is None or ez is None: raise PreventUpdate
      this = update_electrodePositions("PUT",this,e_select,ex,ey,ez)
      return this

    elif clicked=="elec-w" or clicked=="elec-h" or clicked=="elec-d":
      this = update_electrodeDimension("PUT",this,e_select,ew,eh,ed)
      return this # Ignoring ElectrodeKind for now  

    elif clicked=="outer-x" or clicked=="outer-y" or clicked=="outer-z":
      this = update_device_carrier("PUT",this,cx,cy,cz)

    return this


  # Simple update to electrode drop-down
  @app.callback([Output("elec-dropdown","value"),
                 Output("elec-dropdown","options")],
                 Input("device-dropdown","value"),
                 Input("add-elec","n_clicks"),
                 Input("rem-elec","n_clicks"),
                 State("device-json","data"),
                 State("elec-dropdown","value") )
  def update_elec_dropdown(ud,ua,ur,this,e_select):

    if this is None: raise PreventUpdate
    clicked = which_input()

    if isinstance(this,list) and this: this = this[0]
    try: numel = range(0,len(this["array"]["ElectrodeTypeIndex"]))
    except:
      print("\nPossibly bad device data:")
      print(this)
      print("\n")
      raise PreventUpdate

    these = [{"label":"Elec{}".format(n+1),"value":n} for n in numel]

    if clicked == "add-elec": 
        return [len(this["array"]["ElectrodeTypeIndex"])-1], these
    if clicked == "rem-elec": 
        return [],these
    if not these: return [],these
    return [],these


  @app.callback([Output("elec-x","value"),
                 Output("elec-y","value"),
                 Output("elec-z","value"),
                 Output("elec-w","value"),
                 Output("elec-h","value"),
                 Output("elec-d","value"),
                 Output("outer-x","value"),
                 Output("outer-y","value"),
                 Output("outer-z","value"),
                 Output("outer-y","disabled")],
                [Input("device-dropdown","value"),
                 Input("elec-dropdown","value"),
                 State("device-json","data"),
                 State("device-family-dropdown","value"),
                 State("elec-x","value"),
                 State("elec-y","value"),
                 State("elec-z","value"),
                 State("elec-w","value"),
                 State("elec-h","value"),
                 State("elec-d","value"),
                 State("outer-x","value"),
                 State("outer-y","value"),
                 State("outer-z","value")] )
  def update_elec_props(device,e_select,this,family,ex,ey,ez,ew,eh,ed,cx,cy,cz):

    if device is None: raise PreventUpdate
    
    if this is None: 
      this = user_files.get_device_json(device,family)
      # using memoise to reduce waste here 
    elif "ui" not in this:
      this = user_files.get_device_json(device,family)
    
    if isinstance(this,list) and this: this = this[0]

    if not this["ui"]["device"] == device or \
       not this["ui"]["family"] == family:
      this = user_files.get_device_json(device,family)

    cy_disable = False

    ex,ey,ez      = update_electrodePositions("GET",this,e_select)
    ew,eh,ed      = update_electrodeDimension("GET",this,e_select)
    cx,cy,cz, cy_type = update_device_carrier("GET",this)
    if cy_type == "cyl": cy_disable = True

    return ex,ey,ez,ew,eh,ed,cx,cy,cz,cy_disable

  
  @app.callback([Output("device-family-dropdown","options"), 
                 Output("device-family-dropdown","value")],
                  Input("upload-device","filename"),
                  Input("upload-device","contents"), 
                  Input("download-device","data"),                  
                  State("navbar-session","value"),   # change to INPUT
                  prevent_initial_call=True)
  def upload_device(name,data,dl_data,session=1):
    if name is None or data is None: raise PreventUpdate    
    clicked = which_input()

    if clicked == "download-device":
      u_list = list_deviceFamilies()
      return u_list,u_list[-1]["value"]

    if clicked == "navbar-session":
      u_list = list_deviceFamilies()
      return u_list,u_list[-1]["value"]

    path = "../data/u/{}/{}/array.json".format(get_user_ID(),session)
    print("Uploading " + name + " --> " + path)
    save_file(path,data)

    u_list = list_deviceFamilies()
    return u_list,u_list[-1]["value"]
  

  @app.callback(Output("download-device","data"),
                Input("btn-save-array","n_clicks"), 
                State("device-json","data"),
                State("navbar-session","value"),
                prevent_initial_call=True)
  def save_array(n_clicks,data,session=1):

    array,json_string = user_files.make_ARRAY_json(data)
    path = "../data/u/{}/{}/array.json".format(get_user_ID(),session)

    with open(path,"wt") as f:
      f.write(json_string)

    return dict(content=json_string, filename="array.json")


  @app.callback([Output("nerve-x","value"), 
                 Output("nerve-y","value"), 
                 Output("nerve-r","value"), 
                 Output("nerve-name","children"),
                 Output("nerve-name","style"),
                 Output("axon-pop-dropdown","value"),
                 Output("anatomy-json","data")], 
                 Input("upload-nerve","filename"),
                 Input("upload-nerve","contents"),
                 State("url", "pathname"),   # change to INPUT
                 State("nerve-x","value"),
                 State("nerve-y","value"),
                 State("nerve-r","value"),
                 State("nerve-name","children"),
                 State("nerve-name","style"),
                 State("axon-pop-dropdown","value"),
                 State("anatomy-json","data"),
                 State("navbar-session","value"))
  def upload_nerve(name,data,url,nx,ny,nr,nn,nn_style,ap,anat,session=1):

    USER = get_user_ID()

    if name is None or data is None:

      anat,config = user_files.check_existing_nerve_files(USER,session)

      nn = "last used nerve"
      if anat is None: raise PreventUpdate
      if config is not None and "nerve" in config:
        config = config["nerve"]

        if "xRotate" in config: nr = config["xRotate"]
        else:                   nr = 0
        if "xMove" in config: 
          nx = config["xMove"][0]
          ny = config["xMove"][1]
        else: 
          nx = 0
          ny = 0          
        if "uiAxonPop" in config: ap = config["uiAxonPop"]
        if "uifileName" in config: nn = config["uifileName"]
      
      nn_style = {"display":"block"}

      return nx,ny,nr,nn,nn_style,ap,anat

    # name and data defined
    clicked = which_input()

    if clicked == "url":
      anat,config = user_files.check_existing_nerve_files(USER,session)

      nn = "last used nerve"
      if anat is None: raise PreventUpdate
      if config is not None and "nerve" in config:
        config = config["nerve"]

        if "xRotate" in config: nr = config["xRotate"]
        else:                   nr = 0
        if "xMove" in config: 
          nx = config["xMove"][0]
          ny = config["xMove"][1]
        else: 
          nx = 0
          ny = 0          
        if "uiAxonPop" in config: ap = config["uiAxonPop"]
        if "uifileName" in config: nn = config["uifileName"]


    if clicked == "download-file":
      u_list = list_deviceFamilies()
      return u_list,u_list[-1]["value"]

    file_ext = os.path.splitext(name)

    if file_ext[1] in [".json"]:
      path = "../data/u/{}/{}/nerve.json".format(USER,session)
    elif file_ext[1] in [".xml"]: 
      path = "../data/u/{}/{}/nerve.xml".format(USER,session)      
      nx=0
      ny=0
      nr=0
      nn=os.path.basename(file_ext[0])
      nn_style = {"display":"block"}
    else:
      print("{}: not a valid file (expected .xml, .json)".format(name))
      raise PreventUpdate

    print("Uploading " + name + " --> " + path)
    save_file(path,data)
    
    if file_ext[1] in [".json"]:

      config = user_files.parse(path)        
      config = config["nerve"]

      if "xRotate" in config: nr = config["xRotate"]
      else:                   nr = 0
      if "xMove" in config: 
        nx = config["xMove"][0]
        ny = config["xMove"][1]
      else: 
        nx = 0
        ny = 0
      if "uiAxonPop" in config: ap = config["uiAxonPop"]
    else:
      anat = user_files.get_user_XML_contours(USER,session)

    return nx,ny,nr,nn,nn_style,ap,anat


  @app.callback([Output("nerve-json","data"),
                 Output("nerve-x-lbl","children"),
                 Output("nerve-y-lbl","children"),
                 Output("nerve-r-lbl","children")],
                 Input("nerve-x","value"),
                 Input("nerve-y","value"),
                 Input("nerve-r","value"),
                 State("nerve-name","children"))
  def update_nerve_sliders(nx,ny,nr,name):
    dat = {"nerve":{"source":name,"xRotate":nr,"xMove":[nx,ny]}} # other parameters generated at compile
    nxl = "{} µm".format(nx)
    nyl = "{} µm".format(ny)
    nrl = "{}°".format(nr)
    return dat,nxl,nyl,nrl


  @app.callback([Output("nerve-x","min"),
                 Output("nerve-x","max"),
                 Output("nerve-y","min"),
                 Output("nerve-y","max")],
                 Input("anatomy-json","data"))
  def update_nerve_slider_range(anat):
    
    if anat is None: raise PreventUpdate
    if "anat" in anat: anat = anat["anat"]

    x0 = min([obj["xr"][0] for obj in anat])
    x1 = max([obj["xr"][1] for obj in anat])
    x0 = max([2*x1-x0,x1-2*x0])

    return -x0,x0,-x0,x0


  @app.callback(Output("download-nerve","data"),
                Input("btn-save-nerve","n_clicks"), 
                State("nerve-json","data"),
                State("device-json","data"),
                State("axon-pop-dropdown","value"),
                State("navbar-session","value"),
                prevent_initial_call=True)
  def save_nerve_json(n_clicks,nerve,array,axons,session=1):

    json_string = user_files.save_json_files(array,nerve,axons,session)
    return dict(content=json_string, filename="nerve.json")

    



  @app.callback(Output("btn-run","color"),
                Input("btn-run","n_clicks"),
                State("nerve-json","data"),
                State("device-json","data"),
                State("axon-pop-dropdown","value"),
                State("navbar-session","value"),
                prevent_initial_call=True)
  def run_model(n_clicks,nerve,array,axons,session=1):
    if not n_clicks: raise PreventUpdate

    # Step 1: save needed files 
    user_files.save_json_files(array,nerve,axons,session)
    USER = 1
    
    # filename = ['','','','']
    # filename[0] = "../data/u/{}/{}/array.json".format(USER,session)
    # filename[1] = "../data/u/{}/{}/nerve.xml".format(USER,session)
    # filename[2] = "../data/u/{}/{}/nerve.json".format(USER,session)
    # filename[3] = "../data/share/axon/{}.mat".format(axons)
    
    cfg = osparc_api.cfg
    input_file_1, input_file_2, input_file_3, input_file_4 = osparc_api.upload_files(cfg, axons, session, USER)
    

    # step 2: 
    # run nerve_mesh
    nm_results, nm_download_path = osparc_api.run_node(cfg, \
                     [input_file_1, input_file_2, input_file_3],\
                     "simcore/services/comp/nerve-mesh", \
                     "1.0.2")
            
    # run axon-population
    ap_results, ap_download_path = osparc_api.run_node(cfg, \
                     [input_file_4, input_file_2, input_file_3],\
                     "simcore/services/comp/axon-population", \
                     "1.0.2")
            
    # run fields-solver
    es_results, es_download_path = osparc_api.run_node(cfg, \
                     [nm_results],\
                     "simcore/services/comp/eidors-solver", \
                     "2.0.2")
            
    # nerve recording
    nr_results, nr_download_path = osparc_api.run_node(cfg, \
                     [ap_results, es_results],\
                     "simcore/services/comp/nerve-recording", \
                     "1.0.1")
                
    return "danger" 

#%%

if __name__ == "__main__":
  # run in debug modew
  import webpage
  webpage.app.run_server(debug=True)

