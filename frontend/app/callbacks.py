
import dash_core_components as dcc
import dash_html_components as html
import dash

from urllib.parse import quote as urlquote
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from glob import glob
import functools 

import json
import os
import math

import generate_user_files


def parse(text):
    try:
        return json.load(text)
    except ValueError as e:
        print('invalid json: %s' % e)
        return None # or: raise

def file_download_link(filename):
    """Create a Plotly Dash 'A' element that downloads a file from the app."""
    location = "/download/{}".format(urlquote(filename))
    return html.A(filename, href=location)

def which_input(ctx = dash.callback_context):    
    if not ctx.triggered: raise PreventUpdate    
    else: return ctx.triggered[0]['prop_id'].split('.')[0]

def save_file(filepath,content):
    """Decode and store a file uploaded with Plotly Dash."""
    data = content.encode("utf8").split(b";base64,")[1]
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
def db_query(tag,table='USER'):
    print('TODO: query column {} from table {}'.format(tag,table))

# get user ID 
def get_user_ID(app):
    return 1

# get session ID 
def get_session_ID(app):
    return 1

# what device family did the user last select
def get_default_dfamily(app):
    db_query('dfamily','DEFAULTS')
    return None

# what nerve class did the user last select
def get_default_nerveClass(app):
    db_query('nerveclass','DEFAULTS')
    return None

# what nerve class did the user last select
def get_default_runMode(app):
    db_query('runmode','DEFAULTS')
    return 'full'

# what is user's email?
def get_default_email(app):
    db_query('email','USER')
    return None


#%%
#
# Get menu options for device and deviceFamily lists 
# 
#%%

# get list of user-defined devices
def list_userDevices():
  my_list = glob('../data/u/{}/*/array.json'.format(get_user_ID()))
  
  dev_list = []

  if not my_list: return []
  for filepath in my_list:
    with open(filepath) as f:
      array = parse(f)

    if isinstance(array,list): array = array[0]
    try:
      dev_list.append({"label":me_array['array']['name'],"value":filepath})
    except:
      print("error parsing user JSON array: "+filepath)

  return dev_list

# Load device families from /data/share
def list_deviceFamilies():
  with open(r'../data/share/array/index.json') as f:
    arrays = parse(f)

  dflist = set([a['family'] for a in arrays['list']])
  dflist = ([{"label":a,"value":a} for a in dflist])

  mylist = list_userDevices()
  if mylist: dflist.append({"label":"My devices","value":"user"})

  # todo append my-devices    
  return dflist

# Load devices given family from /data/share
def list_devices(family):
  # print(arrays)
  if family is None: return []
  if family == 'user': 
    return list_userDevices()

  with open(r'../data/share/array/index.json') as f:
    arrays = parse(f)

  return [{"label":a['name'],"value":a['file']} 
          for a in arrays['list'] if a['family'] == family]

# Load nerve classes from /data/share
def list_nerveClasses():
    with open(r'../data/share/axon/index.json') as f:
      axons = parse(f)
    # print(axons)
    return([{"label":a["label"],"value":a["value"]} for a in axons['list']])


#%% 
# 
# Callbacks for interacting with json representation of array
#
#%%

def update_electrodePositions(mode,this,e_select,ex=None,ey=None,ez=None):

  try:
    eti = this['array']['ElectrodeTypeIndex']
  except:
    print("\nIssue getting ElectrodeTypeIndex:")
    print(this['array'])
    print("\n")

  if mode == 'WHAT': return eti
  if mode == 'PUT': # update 'this' from ex ey ez
    if ex is None or ey is None or ez is None: raise PreventUpdate
    if not e_select: raise PreventUpdate

    s = e_select[0]
    try:
      this['array']['ElectrodePositions'][s] = [ex,ey,ez]
    except:
      print("\nIssue setting ElectrodePositions:")
      print(this['array']['ElectrodePositions'])
      print(e_select)
      print(s)
      print("\n")
      raise PreventUpdate

    return this
  if mode == 'GET': # update ex ey ez from 'this'

    print('GET elec-xyz:')
    print(e_select)      

    if len(e_select) == 1:
      s = e_select[0]
      ex = this['array']['ElectrodePositions'][s][0]
      ey = this['array']['ElectrodePositions'][s][1]
      ez = this['array']['ElectrodePositions'][s][2]
    else:
      ex = None
      ey = None
      ez = None

    return ex,ey,ez

def update_electrodeDimension(mode,this,e_select,ew=None,eh=None,ed=None):

  these = lambda v : [v[u] for u in range(0,len(v)) if u in e_select ]
  those = lambda v : [v[u] for u in range(0,len(v)) if u not in e_select ]

  try:
    eti = this['array']['ElectrodeTypeIndex']
  except:
    print("\nIssue getting ElectrodeTypeIndex:")
    print(this['array'])
    print("\n")

  single_etype = ( min(eti) == 1 and max(eti) == 1 )
  single_edata = ( not isinstance(this['array']['ElectrodeDimensions'][0],list) )

  if mode == 'WHAT': return single_etype, single_edata
  if mode == 'PUT': 
    if ew is None or eh is None or ed is None: raise PreventUpdate
    if not e_select: 
      if single_etype and single_edata:
        this['array']['ElectrodeDimensions'][0] = ew
        this['array']['ElectrodeDimensions'][2] = eh
        this['array']['InsetDepth'] = ed
      elif single_etype:
        this['array']['ElectrodeDimensions'][0][0] = ew
        this['array']['ElectrodeDimensions'][0][2] = eh
        this['array']['InsetDepth'] = ed
      else: raise PreventUpdate  
      return this,ew,eh,ed
    else: # one or more electrodes selected 
      print('PUT elec-dim:')
      print(e_select)

      if min(these(eti)) == max(these(eti)): # 'these' all same ETI
        new_eti = these(eti)[0] # proposed index 
        if new_eti in those(eti): # some other things also have the proposed ID
          for n in range(1,max(eti)+1): 
            if n not in those(eti): # find lowest ID not already used
              new_eti = n 
              break
      else: # "these" do not all have the same index 
        for n in range(1,max(eti)+1): 
          if n not in those(eti): # find lowest ID not already used
            new_eti = n 

      if single_edata: # make into a list
        this['array']['ElectrodeDimensions'] = [this['array']['ElectrodeDimensions']]

      ek = this['array']['ElectrodeDimensions'][0][1]
      if new_eti >= len(this['array']['ElectrodeDimensions']):
          new_eti = len(this['array']['ElectrodeDimensions'])
          this['array']['ElectrodeDimensions'].append([ew,ek,eh])
      else: 
          this['array']['ElectrodeDimensions'][new_eti] = [ew,ek,eh]

      for e in e_select: # assign ETI
        this['array']['ElectrodeTypeIndex'][e] = new_eti
    # PUT suceeded in updating 'this' from ew eh ed
    return this,ew,eh,ed

  if mode == 'GET':
    print('GET elec-dim:')
    print(e_select)

    ed = this['array']['InsetDepth']
    if single_edata:
      ew = this['array']['ElectrodeDimensions'][0]
      eh = this['array']['ElectrodeDimensions'][2]
    elif single_etype:
      ew = this['array']['ElectrodeDimensions'][0][0]
      eh = this['array']['ElectrodeDimensions'][0][2]
    elif not e_select: 
      ew = None
      eh = None
    else:
      s = e_select[0]
      ew = this['array']['ElectrodeDimensions'][s][0]
      eh = this['array']['ElectrodeDimensions'][s][2]
    # GET suceeded in updating ew eh ed from 'this'
    return ew,eh,ed

def update_device_carrier(mode,this,cx=None,cy=None,cz=None):

  if "c_thickness" in this['array']['carrier']: 
        device_type = 'flat'
  elif "cuff_IDr" in this['array']['carrier'] and \
                     this['array']['carrier']['cuff_IDr'] > 0.45 * \
                     this['array']['carrier']['cuff_IDx'] :
        device_type = 'cyl'
  else: device_type = 'rect'

  if mode == 'WHAT': return device_type
  if mode == 'PUT':
    if device_type == 'flat': 
      this['array']['carrier']['c_len'] = cx
      this['array']['carrier']['c_wid'] = cy
      this['array']['carrier']['c_thickness'] = cz
    elif device_type == 'cyl':
      this['array']['carrier']['cuff_IDx'] = cx
      this['array']['carrier']['cuff_IDy'] = cx
      this['array']['carrier']['cuff_IDr'] = cx * 0.99/2
      this['array']['carrier']['cuff_length'] = cz
    else:
      this['array']['carrier']['cuff_IDx'] = cx
      this['array']['carrier']['cuff_IDy'] = cy
      this['array']['carrier']['cuff_length'] = cz

    return this,device_type

  if mode == 'GET':
    try:
      if device_type == 'flat': 
        cx = this['array']['carrier']['c_len']
        cy = this['array']['carrier']['c_wid']
        cz = this['array']['carrier']['c_thickness']
      else:      
        cx = this['array']['carrier']['cuff_IDx']
        cy = this['array']['carrier']['cuff_IDy']
        cz = this['array']['carrier']['cuff_length']
    except:
        print("\nPossibly bad device data:")
        print(this['array']['carrier'])
        print("\n")        
    
    return cx,cy,cz,device_type

#%%
# 
# Load device from filesystem 
@functools.lru_cache(maxsize=32)
def get_device_json(name,family):

    is_user = (family == 'user')
    print("Loading '{}'".format(name))

    if is_user: # user specified device path 
      if not os.path.isfile(name):
        with open(name,'wt') as f: 
          json.dump(name, outfile)        
      else:
        with open(name,'rt') as f:
          this = parse(f)
    else: 
      file = r'../data/share/array/{}.json'.format(name)      
      with open(file,'rt') as f:
        this = parse(f)

    this['ui'] = {'device':name,'family':family}
    return this


#%% 
#
# DEFINE CALLBACKS
# 
#%%
def add_callbacks(app):

  # Callback to show/hide "firing rate" field
  @app.callback( Output("div-spike-rate","style"), Input("run-mode-dropdown","value") )
  def toggle_show_spikerate(run_mode):        
    if run_mode == 'full':  return {'display':'block'}
    else:                   return {'display':'none'}

  # Callback to set up device dropdown "firing rate" field
  @app.callback( [Output("device-dropdown","options"), 
                  Output("device-dropdown","value"), 
                  Output("device-dropdown","disabled")], 
                  Input("device-family-dropdown","value") )
  def update_device_dropdown(family):
    opts = list_devices(family)
    default = [{"label":"...","value":""}]
    if opts is None: return default, None, True
    if not opts: return default, None, True
    return opts, opts[0]['value'], False


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
      this = get_device_json(device,family)
      was_loaded = True

    if isinstance(this,list) and this: this = this[0]
    clicked = which_input()

    print("update_device: "+clicked)

    if clicked=='device-dropdown' and not was_loaded: # got a new device?
        this = get_device_json(device,family)
        was_loaded = True
        return this

    these = lambda v : [v[u] for u in range(0,len(v)) if u in e_select ]
    those = lambda v : [v[u] for u in range(0,len(v)) if u not in e_select ]

    if clicked=='add-elec': 
      if e_select:
            eti = this['array']['ElectrodeTypeIndex'][e_select[0]]
      else: eti = this['array']['ElectrodeTypeIndex'][-1]

      this['array']['ElectrodePositions'].append([0,0,0])
      this['array']['ElectrodeTypeIndex'].append(eti)

      if 'ElectrodeAngle' in this['array']: 

        if e_select:
              eti = this['array']['ElectrodeTypeIndex'][e_select[0]]
        else: eti = this['array']['ElectrodeTypeIndex'][-1]
        this['array']['ElectrodeTypeIndex'].append(eti)
        return this

    elif clicked=='rem-elec': 
      if not e_select: raise PreventUpdate

      this['array']['ElectrodePositions'] = those(this['array']['ElectrodePositions'])
      this['array']['ElectrodeTypeIndex'] = those(this['array']['ElectrodePositions'])
      if 'ElectrodeAngle' in this['array']: 
        this['array']['ElectrodeAngle'] = those(this['array']['ElectrodeAngle'])

      return this

    elif clicked=='elec-x' or clicked=='elec-z' or clicked=='elec-z':
      if not e_select: raise PreventUpdate
      if ex is None or ey is None or ez is None: raise PreventUpdate
      this = update_electrodePositions('PUT',this,e_select,ex,ey,ez)
      return this

    elif clicked=='elec-w' or clicked=='elec-h' or clicked=='elec-d':
      this = update_electrodeDimension('PUT',this,e_select,ew,eh,ed)
      return this # Ignoring ElectrodeKind for now  

    elif clicked=="outer-x" or clicked=="outer-y" or clicked=="outer-z":
      this = update_device_carrier('PUT',this,cx,cy,cz)

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
    try: numel = range(0,len(this['array']['ElectrodeTypeIndex']))
    except:
      print("\nPossibly bad device data:")
      print(this)
      print("\n")
      raise PreventUpdate

    these = [{"label":"Elec{}".format(n+1),"value":n+1} for n in numel]

    if clicked == 'add-elec': 
        return [len(this['array']['ElectrodeTypeIndex'])-1], these
    if clicked == 'rem-elec': 
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
      this = get_device_json(device,family)
      # using memoise to reduce waste here 
    elif "ui" not in this:
      this = get_device_json(device,family)
    
    if isinstance(this,list) and this: this = this[0]

    if not this['ui']['device'] == device or \
       not this['ui']['family'] == family:
      this = get_device_json(device,family)

    cy_disable = False

    ex,ey,ez      = update_electrodePositions('GET',this,e_select)
    ew,eh,ed      = update_electrodeDimension('GET',this,e_select)
    cx,cy,cz, cy_type = update_device_carrier('GET',this)
    if cy_type == 'cyl': cy_disable = True

    return ex,ey,ez,ew,eh,ed,cx,cy,cz,cy_disable


  @app.callback([Output("device-family-dropdown","options"), 
                 Output("device-family-dropdown","value")],
                  Input("upload-device","filename"),
                  Input("upload-device","contents"), 
                  Input("download-file","data"), 
                  prevent_initial_call=True)
  def upload_device(upload_filename,upload_data,download_data):
    if upload_filename is None or upload_data is None: raise PreventUpdate    
    clicked = which_input()

    if clicked == "download-file":
      u_list = list_deviceFamilies()
      return u_list,u_list[-1]['value']      


    for name, data in zip(upload_filename, upload_data):

      path = '../data/u/{}/{}/array.json'.format(get_user_ID(),get_session_ID())

      print('Uploading ' + name + ' --> ' + path)
      save_file(path,data)

    u_list = list_deviceFamilies()
    return u_list,u_list[-1]['value']

  #
  @app.callback(Output("download-file","data"),                 
                Input("btn-save-array","n_clicks"), 
                State("device-json","data"),
                prevent_initial_call=True)
  def make_local_array(n_clicks,data):

    array,json_string = generate_user_files.mk_array(data)
    path = '../data/u/{}/{}/array.json'.format(get_user_ID(),get_session_ID())

    with open(path,'wt') as f:
      f.write(json_string)

    return dict(content=json_string, filename="array.json")



  print('Inserted UI callbacks')







#%%
if __name__ == '__main__':

  import webpage
  webpage.app.run_server(debug=True)

