#%%

from flask import Flask, send_from_directory
from dash.dependencies import Input, Output, State

import os
import re
import dash
import dash_core_components as dcc
import dash_html_components as html #elements of wireframe
import dash_bootstrap_components as dbc #grid
import dash_extensions as dec
#import plotly.graph_objs as go
#import pandas as pd


#%%
# Normally, Dash creates its own Flask server internally. By creating our own,
# we can create a route for downloading files directly:
server = Flask(__name__)
app = dash.Dash(server=server, external_stylesheets=[dbc.themes.BOOTSTRAP], 
                suppress_callback_exceptions=True)
# ,             external_scripts=["//daybrush.com/moveable/release/latest/dist/moveable.min.js"])

import page_setup
import page_results
import graphics
import user_files

# Since we're adding callbacks to elements that don't exist in the app.layout,
# Dash will raise an exception to warn us that we might be
# doing something wrong.
# In this case, we're adding the elements through a callback, so we can ignore
# the exception.

#%%
# Define path so that we can serve downloads 
@server.route("/download/<path:path>")
def download(path):
    """Serve a file from the upload directory."""
    return send_from_directory(UPLOAD_DIRECTORY, path, as_attachment=True)


#%% 


nav_bar = dbc.Card([
  dbc.Row([dbc.Col([
    dbc.Button("MODEL SETUP", id='navbar-setup', color="secondary", outline=True,
                       className="mr-1",n_clicks=0,size="sm"), 
    dbc.Button("MODEL RESULTS", id='navbar-results', color="secondary", outline=True,
                       className="mr-1",n_clicks=0,size="sm",disabled=True),
   ],width=9,style={'text-align':'left'}),
   dbc.Col([
        dcc.Dropdown(id="navbar-session", 
                     options = user_files.list_userSessions(), value = 1, 
                     placeholder="Select Session", persistence=True, clearable=False), 

    ],width=3)
  ])],body=True, style={"background-color": "#eeeeee","padding":"6px","margin-top":"4px","margin-bottom":"4px"})


app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    html.Div(id='page-content'),
    dcc.Store(id="user-data",storage_type='local'),
    dcc.Store(id="device-json", storage_type='session'),
    dcc.Store(id="nerve-json", storage_type='session'),
    dcc.Store(id="anatomy-json", storage_type='memory'),
    dcc.Store(id="results-json", storage_type='memory'),
])

#%%

page_setup.add_callbacks(app)
page_results.add_callbacks(app)
graphics.add_callbacks(app)

# middle_layer.add_callbacks(app)

#%%
# implement routing
@app.callback(Output('page-content', 'children'), 
              Input('url', 'pathname'))
def display_page(url):

  if url.startswith('/setup'):
      return [nav_bar, page_setup.layout]
  elif url.startswith('/results'):
      return [nav_bar, page_results.layout]
  else:
      return nav_bar
    # could also return a 404 "URL not found" page here

@app.callback([Output('url','pathname'),
               Output('navbar-session','options'),
               Output('navbar-results','disabled')],
               Input('navbar-setup','n_clicks'),
               Input('navbar-results','n_clicks'),
               Input('navbar-session','value'),
               State('navbar-session','options'),
               State('url','pathname'))
def on_nav_click(bs,br,session,ses_list,url):

  clicked = page_setup.which_input()
  # print("update_navbar: "+clicked)

  if session == ses_list[-1]['value']: # "new session"

    s_list = user_files.list_userSessions()
    ses_list[-1]['label'] = 'Session {}'.format(session)
    ses_list.append({'label':'New Session','value':session+1})
    session_path = '../data/u/{}/{}/'.format(1,session)
    if not os.path.exists(session_path):
      os.mkdir(session_path)

  r_ok = user_files.has_results(1,session)    

  if clicked == 'navbar-setup': url = "/setup/{}".format(session)
  elif clicked == 'navbar-results': 
    if r_ok: url = "/results/{}".format(session)
    else:    url = "/setup/{}".format(session)
  else: 
     stub = re.match(r"/[^/]/",url)
     if not stub: stub = ['/setup/']
     url = "{}{}".format(stub[0],session)
     

  return url, ses_list, not r_ok


#%%
if __name__ == '__main__':
    app.run_server(debug=False)



