
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

def list_devices():

    with open(r'../data/share/array/index.json') as f:
      array = parse(f)

    print(array)
    return [{"label":"D1","value":"filename"}]


def list_nerveClasses():

    with open(r'../data/share/axon/index.json') as f:
      array = parse(f)

    print(array)
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



