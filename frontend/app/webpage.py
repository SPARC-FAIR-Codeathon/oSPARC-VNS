#%%

from flask import Flask, send_from_directory
from dash.dependencies import Input, Output, State

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
                     options = user_files.list_userSessions(), value = [], 
                     placeholder="Select Session", persistence=True, clearable=False), 

    ],width=3)
  ])],body=True, style={"background-color": "#eeeeee","padding":"6px","margin-top":"4px","margin-bottom":"4px"})


app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    dbc.Container(html.Div(id='page-content')),
    dcc.Store(id="user-data",storage_type='local'),
    dcc.Store(id="device-json", storage_type='session'),
    dcc.Store(id="nerve-json", storage_type='session'),
    dcc.Store(id="anatomy-json", storage_type='memory'),
    dcc.Store(id="results-json", storage_type='memory'),
])

#%%

callbacks.add_routing(app)
callbacks.add_callbacks(app)
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

#%%
if __name__ == '__main__':
    app.run_server(debug=False)



