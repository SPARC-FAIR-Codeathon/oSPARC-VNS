

# raghunath
# Jul '17
# This should work.

import dash
import dash_core_components as dcc
import dash_html_components as html

from dash.dependencies import Input, Output, State


import os
import glob

UPLOAD_DIRECTORY = ".\\data"

if not os.path.exists(UPLOAD_DIRECTORY):
       os.makedirs(UPLOAD_DIRECTORY)

ds_names = glob.glob('.\\data\\*') + ["... new dataset"]

app = dash.Dash()

fnameDict = {'chriddy': ['opt1_c', 'opt2_c', 'opt3_c'], 
               'jackp': ['opt1_j', 'opt2_j']}

names = list(fnameDict.keys())
nestedOptions = fnameDict[names[0]]

app.layout = html.Div(
    [
        html.Div([
        dcc.Dropdown(
            id='name-dropdown',
            options=[{'label':name, 'value':name} for name in names],
            value = list(fnameDict.keys())[0]
            ),
            ],style={'width': '20%', 'display': 'inline-block'}),
        html.Div([
        dcc.Dropdown(
            id='opt-dropdown',
            ),                                                                                                                                                                                                                                                             
            ],style={'width': '20%', 'display': 'inline-block'}
        ),
        html.Hr(),        
        html.Div(id='display-selected-values'),        
        html.Hr(),
        html.Div([
        dcc.Dropdown(
            id='ds-name-dropdown',
            options= [{'label':n, 'value':n} for n in ds_names],
            ),
            ],style={'width': '20%', 'display': 'block'}
        ),
        html.Div([
        dcc.Input(id='ds-name-input', value=ds_names[0], type='text',debounce=True,pattern=r'[^\s]*')
            ], style={'width': '20%', 'display': 'none'}, hidden=True
        ),
    ]
)

@app.callback(
    Output('opt-dropdown', 'options'),
    [Input('name-dropdown', 'value')]
)
def update_date_dropdown(name):
    return [{'label': i, 'value': i} for i in fnameDict[name]]

@app.callback(
    Output('display-selected-values', 'children'),
    [Input('opt-dropdown', 'value')])
def set_display_children(selected_value):
    return 'you have selected {} option'.format(selected_value)

# I added this one 

@app.callback(
    [Output('ds-name-dropdown', 'hidden'),
     Output('ds-name-input', 'hidden')],
     Input('ds-name-dropdown', 'value'),
     State('ds-name-dropdown','options'))
def set_show_new_ds_name(selected,opts):

    print("SEL = " + str(selected))
    print("OPTS = " + str(opts[-1]['value']))
    if selected is None:
        return False,True

    if selected == opts[-1]['value']:
        print("I want to show a textbox please")
        return True,False
    else: 
        return False,True


@app.callback(
    [Output('ds-name-input', 'hidden'),
     Output('ds-name-dropdown', 'hidden'),
     Output('ds-name-dropdown', 'options')],
     Input('ds-name-input', 'value'), 
     Input('ds-name-input','n_clicks'),
     State('ds-name-dropdown','options'),      
     prevent_initial_call=True)
def set_hide_new_ds_name(input_value,n_clk,old_opts):

    if n_clk == 0:
        print ('prevent_initial_call failed')
        return True,False,old_opts

    print("NEW NAME = " + str(input_value))


    # if input_value is already in list 
    ds_list = [u['value'] for u in old_opts]

    if input_value in ds_list: 
        return True,False, old_opts

    new_upload_dir = os.path.join(UPLOAD_DIRECTORY,'/'+input_value+'/')
    if not os.path.exists(new_upload_dir):
       os.makedirs(new_upload_dir)

    ds_list = ds_list[0:-1] + [input_value] + ds_list[-1:]
    return True,False, [{'label':u, 'value':u} for u in ds_list]






if __name__ == '__main__':
    app.run_server()

