
# model_info.py

import yaml

# import dash
# import dash_core_components as dcc
# import dash_bootstrap_components as dbc
# import dash_html_components as html
from dash.dependencies import Input, Output, State

with open("./conf/models.yaml") as fi: 
	model_list = list(yaml.safe_load(fi).values())

def get_menu():
	return [{'label':"Version {}".format(m['index']), 
			 'value':m['index']} for m in model_list]

def run_unit_test():

	print('===============================')
	print(' UNIT TEST: model_info.py')
	print('===============================')

	print('Model List:')

	print(model_list)
	print(model_list.__class__)

	print('TODO')


def add_callbacks(app):
	# print('@app.callback add table_callbacks')

	@app.callback(
	    [Output("opts-m1", "value"),
	     Output("opts-m2","value")],
	     Input("opts-ver", "value")
	)
	def update_mdl_vars(model_version):
	    """Save uploaded files and regenerate the file list."""

	    m = model_list[model_version-1]
	    return m['patch_pixel_size'], m['patch_size']

	    # files = uploaded_files()
	    # if len(files) == 0:
	    #     return [html.Li("No files yet!")]
	    # else:
	    #     return [html.Li(file_download_link(filename)) for filename in files]

if __name__ == "__main__":
    run_unit_test()