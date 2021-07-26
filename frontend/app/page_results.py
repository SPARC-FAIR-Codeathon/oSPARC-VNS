
import os

from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

import dash
import dash_core_components as dcc
import dash_html_components as html #elements of wireframe
import dash_bootstrap_components as dbc #grid
import dash_extensions as dec
import plotly.express as px

import numpy
import pandas as pd
import scipy.io as spio
from glob import glob

import user_files


def which_input(ctx = dash.callback_context):    
    if not ctx.triggered: raise PreventUpdate    
    else: return ctx.triggered[0]['prop_id'].split('.')[0]

#%% 
#
# Get attributes from imported mat data 
#
#%%

def list_traces(data):
    if data is None: return [{"label":"no traces","value":-1}]
    print('get initial traces from data')


def list_epochs(data):
    if data is None: return [{"label":"no epochs","value":-1}]
    print('get initial epochs from data')

def list_axon_populations(data):
    if data is None: return [{"label":"no axons","value":-1}]
    print('get initial axons from data')


def list_fascicles(data):
    if data is None: return [{"label":"no fascicles","value":-1}]
    print('get initial fascicles from data')

def list_electrodes(data):
    # default = [{"label":"no electrodes","value":-1}]
    default = list()
    if data is None: return default
    if 'electrodes' not in data: return default
    return [{"label":e,"value":e} for e in data['electrodes']]






        # implement page routing
def add_callbacks(app):

  @app.callback(Output("results-json","data"),                
                Input("result-file-dropdown","value"))
  def update_results_json(filename): 
    print('Callback update_results_json')
    if filename is None: return {'filename':None,'type':None}
    return user_files.get_results_file(filename)


  @app.callback(Output("show-elec-controls","style"),
                Output("elec-active","options"),
                Input("results-json","data"))
  def set_electrode_menu(data):
    print('Callback set_electrode_menu')
    opts = list_electrodes(data)
    if data is None: return {"display":"none"}, opts
    return {"display":"block"}, opts


  @app.callback(Output("results-wave","figure"),                
                Input("results-json","data"), 
                Input("elec-active","value"))
  def update_results_figure(data,elec):

    print('Callback update_results_figure')
    if data is None: return dict()    
    if "type" not in data: return dict()
    if data['type'] == 've':
      if not elec: return dict()
      f_list = [f for f in data.keys() if f.startswith("Fascicle")]
      if not f_list: return dict()
      if isinstance(elec,list): elec = elec[0]

      df = pd.DataFrame()
      print(f_list)

      for f in f_list:
        dx = pd.DataFrame.from_dict(data[f])        
        df = df.append(dx, ignore_index=True)

      if df.empty: return dict()
      fig = px.scatter_3d(df, x='x', y='y', z='z',color=elec)

      scene=dict(aspectmode='cube', #this string can be 'data', 'cube', 'auto', 'manual'
           # a custom aspectratio is defined as follows:
           aspectratio=dict(x=1, y=1, z=1)
           )

      fig.update_layout(scene = scene)

      return fig

    else: 
      print("TYPE:'{}' not yet implemented'".format(data['type']))
      return dict()

    

    print(f_list)

    # filter f_list ? 

    return dict()

  @app.callback(Output("elec-active","value"),                
                Input("elec-active","options"),
                State("elec-active","value"))
  def default_elec_select(opts,curr_val):
    if not opts: raise PreventUpdate
    vals = [v['value'] for v in opts]
    if curr_val in vals: return curr_val
    return vals[0]
    

#%%

# first_col is the array device controls 
layout = html.Div([ dbc.Row([ # first row is graphics
  dcc.Graph(figure={},id="results-wave",style={'height':600})], style={"margin":"auto","display":"block","height":"75vh"}),
  # dbc.Row([]) for graphics controls ?
  dbc.Row([
    dbc.Col([dbc.Card([
      dcc.Dropdown(id="result-file-dropdown", 
                   options = user_files.list_resultsFiles(), 
                     placeholder="Select Results File", persistence=True), 
        # html.Div([
        #   dcc.Dropdown(id="result-trace", 
        #                options = list_traces(None), 
        #                placeholder="Select a Trace", persistence=True,multi=False), 
        #   dcc.Dropdown(id="result-epoch", 
        #                options = list_epochs(None), 
        #                placeholder="Select a Trace", persistence=True,multi=False), 
        # ],id="show-epoch-controls",style={"display":"none"}),
      ],body=True)],width=4,), 
    dbc.Col([dbc.Card([html.Div([ "nothing here for now"
        # dcc.Checklist(id="result-axons",options = list_axon_populations(None),value=None),
        # html.Hr(),
        # dcc.Checklist(id="result-fascicles",options = list_fascicles(None),value=None)
        ],id="show-axon-controls",style={"display":"none"}),
      ],body=True)],width=4,), 
    dbc.Col([dbc.Card([html.Div([
        dcc.Dropdown(id="elec-active",options = list_electrodes(None),value=list()),
        # dcc.Dropdown(id="elec-reference",options = list_electrodes(None),value=None)
        ],id="show-elec-controls",style={"display":"none"}),
      ],body=True)],width=4,)
    ],style={"height":"20vh"}) 
  ])



#%%

if __name__ == "__main__":
    # something
    filename = glob('../data/u/1/1/ve*.mat')

    # spio.get_matfile_version(filename[0])

    data = spio.loadmat(filename[0])
    




