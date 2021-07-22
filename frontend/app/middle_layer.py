
# the middle layer is responsible for keeping the 'task' alive
# after the user presses run, moving the SPARC api through the 
# necessary intermediate steps even if the user navigates away 
# from the front-end. The last celery task is to send an email
# to the user with their results. 

from celery import Celery
from dash.dependencies import Input, Output, State

app = Celery('vnsModel',broker='redis://localhost')


@app.task
def add(x, y):
    return x + y



import dash_core_components as dcc

# I'm also putting database-emulating code here 



def add_callbacks(app):

  # Simple update to electrode drop-down
  @app.callback(Output("...","value"),
                Input("btn-run","n_clicks"),
                State("device-json","data"),
                State("nerve-json","data"), 
                State("run-mode-dropdown","value")
                State("user-email","value"))
  def on-run_button_click(nc, array, nerve, mode, email):