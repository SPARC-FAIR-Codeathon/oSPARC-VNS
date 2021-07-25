
'''
import dash_core_components as dcc
import dash_html_components as html
import dash

from urllib.parse import quote as urlquote
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate


import functools 
import base64
import json
import math
import os
import re

import user_files
'''

from dash.dependencies import Input, Output, State

import dash
import dash_core_components as dcc
import dash_html_components as html #elements of wireframe
import dash_bootstrap_components as dbc #grid
import dash_extensions as dec

import os
from glob import glob
import scipy.io as spio

# first_col is the array device controls 
layout = html.Div([ dbc.Row([ # first row is graphics
  dcc.Graph(figure=None,id="results-wave")], style={"height":"75vh"}),   
  # dbc.Row([]) for graphics controls ?
  dbc.Row([
    dbc.Col([dbc.Card([
      dcc.Dropdown(id="result-file-dropdown", 
                     options = user_files.list_resultsFiles(), 
                     placeholder="Select a Device Family", persistence=True), 
        html.Div([
          dcc.Dropdown(id="result-trace", 
                       options = view_results.list_traces(None), 
                       placeholder="Select a Trace", persistence=True,multi=False), 
          dcc.Dropdown(id="result-epoch", 
                       options = view_results.list_epochs(None), 
                       placeholder="Select a Trace", persistence=True,multi=False), 
        ],id="show-epoch-controls",style={"display":"none"}),
      ],body=True)],width=4,), 
    dbc.Col([dbc.Card([html.Div([
        dcc.Checklist(id="result-axons",options = view_results.list_axon_populations(None),value=None),
        html.Hr(),
        dcc.Checklist(id="result-fascicles",options = view_results.list_fascicles(None),value=None)
        ],id="show-axon-controls",style={"display":"none"}),
      ],body=True)],width=4,), 
    dbc.Col([dbc.Card([html.Div([
        dcc.Dropdown(id="wave-elec-active",options = view_results.list_electrodes(None),value=None),      
        dcc.Dropdown(id="wave-elec-reference",options = view_results.list_electrodes(None),value=None)
        ],id="show-elec-controls",style={"display":"none"}),
      ],body=True)],width=4,)
    ],style={"height":"20vh"}) 
  ])



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
    if data is None: return [{"label":"no electrodes","value":-1}]
    print('get electrodes from data')













        # implement page routing
def add_callbacks(app):

 


    pass








if __name__ == "__main__":
    # something
    filename = glob('../data/u/1/1/nerve-recording*.mat')

    spio.get_matfile_version(filename[0])

    data = spio.loadmat(filename[0])
    

