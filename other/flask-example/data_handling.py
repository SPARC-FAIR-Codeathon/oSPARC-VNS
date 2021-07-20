
# data_handling.py

from glob import glob
import pandas as pd
import os
import regex as re # regular expressions 

import base64

import dash
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

import dash_table

# Define where data is stored 
UPLOAD_DIRECTORY = os.path.join("data","user")
# This will contain ~/dataset/primary/sub-xx/sam-xx/...
# This will contain ~/dataset/staging/...

if not os.path.exists(UPLOAD_DIRECTORY):
       os.makedirs(UPLOAD_DIRECTORY)

def get_path(dsName=None,sub=None,sam=None,ext=None): 

    if dsName is None: 
        return UPLOAD_DIRECTORY

    p = os.path.join(UPLOAD_DIRECTORY,dsName)

    if sub is None:
        return os.path.join(p,'staging')
    else:
    	p = os.path.join(p,'primary')
    if not isinstance(sub,str): 
        sub = "sub-{}".format(sub)

    p = os.path.join(p,sub)

    if sam is None:
        return p
    if not isinstance(sam,str): 
        sam = "sam-{}".format(sam)

    p = os.path.join(p,sam)

    if ext is None:
        return p

    ext = "{}_{}_{}".format(sub,sam,ext)

    return os.path.join(p,ext)

def parse_filename(file=None,dsName=None):
    # return subject sample
    
    digits_of_ = lambda f : re.split(r'[^\d]+',f+' 1')
    sub = None
    sam = None

    if file is None: 
        return
    # Option 1 - first value sub*/d
    # Option 2 - first value > 10
    vals = re.search(r"(?<=sub[^\d]*)\d+",file)
    if vals is None: # re.search found nothing
        vals = [ int(s) for s in digits_of_(file) if len(s) > 0 ]
    else: 
        vals = [ int(vals[0]) ]
    try: 
        sub = next(x for x in vals if x > 20) 
    except StopIteration: # if no vals > 20, get first ##
        sub = next(x for x in vals)
    
    if sub is None: sub = 1 # default value
    
    vals = re.search(r"(?<=sam[^\d]*)\d+",file)
    if vals is None:
        vals = glob(get_path(dsName,sub,"*"))
        sam = len(vals) + 1
    else:
        sam = int(vals[0])
    
    if sam is None: sam = 1 # default value
    return sub,sam

def get_ds_names(and_new=False):

    ds_names = glob(os.path.join(UPLOAD_DIRECTORY,"*")) 
    ds_names = [n.split(os.path.sep)[-1] for n in ds_names]

    if and_new: 
        ds_names = ds_names + ["... new dataset"]

    return [{'label':n, 'value':n} for n in ds_names]

def get_ds_contents(dsName=None):

    default = [{'img':'','sub':0,'sam':0,'px-um':0.010,'path':''}]
    if dsName is None: # default table row
        return default

    print('DATASET /'+dsName+":")

    contents = glob(get_path(dsName,'*','*')) + glob(get_path(dsName)) 

	print('ERROR get-extra-contents')


    for ii in range(len(contents)): 

        img = glob(os.path.join(contents[ii],'*.tif'))

        if len(img) == 0:
            print('no TIF found in '+contents[ii])
            contents[ii] = ''
            continue
        else: print(img[0])

        img = os.path.split(img[0])
        sub,sam = parse_filename(contents[ii])
        size_file = os.path.join(img[0],'sub-{}_sam-{}_pixel_size.txt'.format(sub,sam))
        pixel = 0.010

        if os.path.exists(size_file):
          with open(size_file,'r') as fi:
            pixel = float(fi.readline())
        elif ii > 0:
            pixel = contents[ii-1]['px-um']

        contents[ii] = {'img':img[1],'sub':sub,'sam':sam,
                        'px-um':pixel,'path':os.path.join(img[0],img[1])}

    contents = [c for c in contents if len(c) > 0]
    if len(contents) > 0: 
        return contents
    else:
        return default

def initialise_table(id='ds-table',name=None): 
    return dash_table.DataTable(id=id,columns=[
            {"id":"img","name":"File","editable":False}, # "presentation":"markdown"
            {"id":"sub","name":"Subject ID","type":"numeric","editable":True},
            {"id":"sam","name":"Sample ID","type":"numeric","editable":True},
            {"id":"px2um","name":"pixel size (µm)","type":"numeric","editable":True}],
           data=get_ds_contents(name),
           style_cell=dict(textAlign='left'),
           row_deletable=True)

# Alicia
def save_file(filepath,content):
    """Decode and store a file uploaded with Plotly Dash."""
    data = content.encode("utf8").split(b";base64,")[1]
    base = os.path.dirname(filepath)
    if not os.path.exists(base): os.makedirs(base)
    with open(filepath, "wb") as fp:
        fp.write(base64.decodebytes(data))

"""@app.callback(
    Output("ds-table", "data"),
    Input("upload-data", "filename"), 
    Input("upload-data", "contents"),
    State("ds-name","value") )"""
def upload_files(upload_filename,upload_data,dsName):
    """Save uploaded files and regenerate the file list."""

    if upload_filename is None and upload_data is None:
        return []

    if dsName is not None: 
        dsName = dsName.strip()

    rows = []

    for name, data in zip(upload_filename, upload_data):
        # sub,sam = parse_filename(name,dsName)
        path = get_path(dsName) # want to send to /staging/
        file_ext = os.path.splitext(name)
        name = file_ext[0]

        if file_ext[1] in ['.tif','.tiff','.png']:
            # name = "sub-{}_sam-{}_Image-EM".format(sub,sam)
            pass
        elif file_ext[1] in ['.txt']: 
            # name = "sub-{}_sam-{}_pixel_size".format(sub,sam)
            pass
        else:
            print("WARNING: {} unknown to mamua".format(name))
        path = os.path.join(path,name) + file_ext[1]

        print('Uploading ' + file_ext[0] + file_ext[1] +
                         ' --> ' + path)
        
        save_file(path,data)

        # get corresponding dataset row
        sub,sam = parse_filename(name)
        rows.append({'img':name,'sub':sub,'sam':sam,'px-um':-1,'path':path})

    return rows


def update_table(new_rows,dat,row_idx)

  if dat: # isempty check  
    for r in new_rows: # each row is a dict 
      if r['px-um'] is -1: 
        r['px-um'] = dat[-1]['px-um']

  if row_idx is None: 
    dat.extend(new_rows)
  else:
    dat[row_idx] = new_rows[0] 



def wrangle_files(dat, dsName):

    if dsName is None: 
        return
    
    dsName = dsName.strip()

    print(dat)
    print('TODO - move files from /staging/ to /primary/')

    for d in dat:

        file_ext = os.path.splitext(d['name'])
        new_path = get_path(dsName,d['sub'],d['sam'],'image-EM.'+file_ext)

        if not os.path.exists(d['path']):
            print('WARNING: Missing File: '+d['path'])
            continue

        if new_path is d['path']: 
            continue

        os.replace(d['path'],new_path) 

    print('CHECK DONE - move files from /staging/ to /primary/')

def add_callbacks(app):

        # files = uploaded_files()
        # if len(files) == 0:
        #     return [html.Li("No files yet!")]
        # else:
        #     return [html.Li(file_download_link(filename)) for filename in files]

    @app.callback(
        Output('ds-hidden','children'),
        Input('ds-table', 'data'),        
        State('ds-table','data_previous')
    )
    def update_dataset(data,data_old):

    	print('TODO update_dataset_from_table')
    	print('DATA = ' + str(data))
    	print('OLD  = ' + str(data_old))

    	if data is None: raise PreventUpdate

    	return data
    


def run_unit_test():

    print('===============================')
    print(' UNIT TEST: data_handling.py')
    print('===============================')
    print('ERROR - TODO')

    parse_filename('sub-401-sam-3-xxx.ext2')
    parse_filename('sub-401-xxx.ext2','example')


if __name__ == "__main__":
    run_unit_test()




''' Extra code from Alicia: 

def save_file(name, content):
    """Decode and store a file uploaded with Plotly Dash."""
    data = content.encode("utf8").split(b";base64,")[1]
    with open(os.path.join(UPLOAD_DIRECTORY, name), "wb") as fp:
        fp.write(base64.decodebytes(data))

def uploaded_files():
    """List the files in the upload directory."""
    files = []
    for filename in os.listdir(UPLOAD_DIRECTORY):
        path = os.path.join(UPLOAD_DIRECTORY, filename)
        if os.path.isfile(path):
            files.append(filename)
    return files

def file_download_link(filename):
    """Create a Plotly Dash 'A' element that downloads a file from the app."""
    location = "/download/{}".format(urlquote(filename))
    return html.A(filename, href=location)

'''    



'''
app.clientside_callback(
    """
    function (input,oldinput) {
        if(JSON.stringify(input) != JSON.stringify(oldinput)) {
            for (i in Object.keys(input)) {
                newArray = Object.values(input[i])
                oldArray = Object.values(oldinput[i])
                if (JSON.stringify(newArray) != JSON.stringify(oldArray)) {
                    entNew = Object.entries(input[i])
                    entOld = Object.entries(oldinput[i])
                    for (const j in entNew) {
                        if (entNew[j][1] != entOld[j][1]) {
                            changeRef = [i, entNew[j][0]] 
                            break        
                        }
                    }
                }
            }
        }
        return changeRef
    }    
    """,
    Output('output-container', 'children'),
    [Input('TableID', 'data')],
    [State('TableID', 'data_previous')]
)
'''