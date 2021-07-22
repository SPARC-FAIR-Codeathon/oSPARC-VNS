
import dash_core_components as dcc
import dash_html_components as html
import dash

from urllib.parse import quote as urlquote
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import json
import os
import math


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




# Load device families from /data/share
def list_deviceFamilies():
    with open(r'../data/share/array/index.json') as f:
      arrays = parse(f)

    dflist = set([a['family'] for a in arrays['list']])
    # todo append my-devices    
    return([{"label":a,"value":a} for a in dflist])

# Load devices given family from /data/share
def list_devices(family):
    with open(r'../data/share/array/index.json') as f:
      arrays = parse(f)
    # print(arrays)
    if family is None: return []
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
# Emulate database functionality 

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
# Load device from filesystem 

def get_device_json(name,is_user,this):

    print("Loading '{}'".format(name))

    if is_user:
      file = r'../data/u/{}/{}/array.json'.format(get_user_ID(None), name)

      if not os.path.isfile(file):
        with open(file,'wt') as f: 
          json.dump(this, outfile)        
      else:
        with open(file,'rt') as f:
          this = parse(f)
    else: 
      file = r'../data/share/array/{}.json'.format(name)      
      with open(file,'rt') as f:
        this = parse(f)

    return this



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
      this = get_device_json(device,family=='user',this)
      was_loaded = True

    clicked = which_input()

    print("update_device: "+clicked)

    if clicked=='device-dropdown' and not was_loaded: # got a new device?
        this = get_device_json(device,family=='user',this)
        return this

    these = lambda v : [v[u] for u in range(0,len(v)) if u in e_select ]
    those = lambda v : [v[u] for u in range(0,len(v)) if u not in e_select ]

    if clicked=='add-elec': 
      if e_select:
            eti = this['array']['ElectrodeTypeIndex'][e_select[0]]
      else: eti = this['array']['ElectrodeTypeIndex'][-1]

      this['array']['ElectrodePositions'].add([0,0,0])
      this['array']['ElectrodeTypeIndex'].add(eti)

      if 'ElectrodeAngle' in this['array']: 

        if e_select:
              eti = this['array']['ElectrodeTypeIndex'][e_select[0]]
        else: eti = this['array']['ElectrodeTypeIndex'][-1]
        this['array']['ElectrodeTypeIndex'].add(eti)
        return this

    if clicked=='rem-elec': 
      if not e_select: raise PreventUpdate

      this['array']['ElectrodePositions'] = those(this['array']['ElectrodePositions'])
      this['array']['ElectrodeTypeIndex'] = those(this['array']['ElectrodePositions'])

      if 'ElectrodeAngle' in this['array']: 
        this['array']['ElectrodeAngle'] = those(this['array']['ElectrodeAngle'])
        return this

    if clicked=='elec-x' or clicked=='elec-z' or clicked=='elec-z':
      if not e_select: raise PreventUpdate
      s = e_select[0]
      this['array']['ElectrodePositions'][s] = [ex,ey,ez]
      return this

    if clicked=='elec-w' or clicked=='elec-h' or clicked=='elec-d':
      
      eti = this['array']['ElectrodeTypeIndex']

      if not e_select: 
        if math.min(eti) == 1 and math.max(eti) == 1 :
              this['array']['ElectrodeDimensions'][0] = ew
              this['array']['ElectrodeDimensions'][2] = eh
              this['array']['InsetDepth'] = ed
        else: raise PreventUpdate
      else:

        if math.min(these(eti)) == math.max(these(eti)): # 'these' all same ETI

          tti = these(eti)[0] # proposed index 
          if tti in those(eti): # some other things also have the proposed ID
            for n in range(1,math.max(eti)+1): 
              if n not in those(eti): # find lowest ID not already used
                tti = n 
                break
        else: 
          for n in range(1,math.max(eti)+1): 
            if n not in those(eti): # find lowest ID not already used
              tti = n 

        for e in e_select: # assign ETI
          this['array']['ElectrodeTypeIndex'][e] = tti

        if not isinstance(this['array']['ElectrodeDimensions'][0],list):
          this['array']['ElectrodeDimensions'] = [this['array']['ElectrodeDimensions']]

        this['array']['ElectrodeDimensions'][tti] = [ew,
                                                     this['array']['ElectrodeDimensions'][0][1], 
                                                     eh]
        # Ignoring ElectrodeKind and InsetDepth here 
        return this 

    if clicked=="outer-x" or clicked=="outer-y" or clicked=="outer-z":

        if "c_thickness" in this['array']['carrier']:

            this['array']['carrier']['c_len'] = cx
            this['array']['carrier']['c_wid'] = cy
            this['array']['carrier']['c_thickness'] = cz
        elif "cuff_IDr" in this['array']['carrier'] and this['array']['carrier']['cuff_IDr'] > 0.45*this['array']['carrier']['cuff_IDx'] :

            this['array']['carrier']['cuff_IDx'] = cx
            this['array']['carrier']['cuff_IDy'] = cx
            this['array']['carrier']['cuff_IDr'] = cx * 0.99/2
            this['array']['carrier']['cuff_length'] = cz

        else:
            this['array']['carrier']['cuff_IDx'] = cx
            this['array']['carrier']['cuff_IDy'] = cy
            this['array']['carrier']['cuff_length'] = cz

    return this


  # Simple update to electrode drop-down
  @app.callback([Output("elec-dropdown","value"),
                 Output("elec-dropdown","options")],
                 Input("device-dropdown","value"),
                 Input("add-elec","n_clicks"),
                 Input("rem-elec","n_clicks"),
                 State("device-json","data"),
                 State("elec-dropdown","value") )
  def update_elec_select(ud,ua,ur,this,e_select):

    if this is None: raise PreventUpdate
    clicked = which_input()

    numel = range(0,len(this['array']['ElectrodeTypeIndex']))
    these = [{"label":"Elec{}".format(n+1),"value":n+1} for n in numel]

    if clicked == 'add-elec': 
        return [len(this['array']['ElectrodeTypeIndex'])-1], these
    if clicked == 'rem-elec': 
        return [],these
    if not these: return [],these
    return [0],these


  print('Inserted basic callbacks')