
import dash_core_components as dcc
import dash_html_components as html
import dash

from dash.dependencies import Input, Output
import json


def parse(text):
    try:
        return json.load(text)
    except ValueError as e:
        print('invalid json: %s' % e)
        return None # or: raise

def list_deviceFamilies():
    with open(r'../data/share/array/index.json') as f:
      arrays = parse(f)

    dflist = set([a['family'] for a in arrays['list']])
    # todo append my-devices    
    return([{"label":a,"value":a} for a in dflist])

def list_devices(family):
    with open(r'../data/share/array/index.json') as f:
      arrays = parse(f)

    if family is None: return []

    return [{"label":a['name'],"value":a['file']} 
            for a in arrays['list'] if a['family'] is family]
    

def list_nerveClasses():

    with open(r'../data/share/axon/index.json') as f:
      axons = parse(f)

    print(axons)
    return [{"label":"Rat Vagus, Cervical","value":"rat-cvn"}]



# default value (user/session-based)
def db_query(tag,table='USER'):
    print('TODO: query column {} from table {}'.format(tag,table))

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



def add_callbacks(app):

  print('TODO add add_callbacks')

  @app.callback( Output("div-spike-rate","style"), Input("run-mode-dropdown","value") )
  def toggle_show_spikerate(run_mode):        
    print('toggle_show_sr')
    if run_mode == 'full':
      return {'display':'block'}
    else: 
      return {'display':'none'}

  @app.callback( [Output("device-dropdown","options"), 
                  Output("device-dropdown","value"), 
                  Output("device-dropdown","enabled")], 
                  Input("device-family-dropdown","value") )
  def update_device_dropdown(family):
    opts = list_devices(family)
    if opts is None: return None, None, False
    return opts, opts[0]['value'], True