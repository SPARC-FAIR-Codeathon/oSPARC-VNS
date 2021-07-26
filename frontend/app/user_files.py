# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 02:35:52 2021

@author: Calvin
"""

import os
import json
import xml.etree.ElementTree as et
import functools

from glob import glob
import scipy.io as spio

import numpy as np
import pandas as pd


#%%
# 
# Load device from filesystem 
#
#%%

def parse(file=None,text=None): # json parse
  
  if text is not None:
    try:
      return json.load(text)
    except ValueError as e:
      print('invalid json: %s' % e)
  else:
    with open(file,'rt') as text:
      try:
        return json.load(text)
      except ValueError as e:
        print('invalid json in {}: {}'.format(file,e))
  return None # or: raise

@functools.lru_cache(maxsize=32)
def get_device_json(name,family):

  is_user = (family == 'user')
  print("Loading '{}'".format(name))

  if is_user: # user specified device path 
    if not os.path.isfile(name):
      with open(name,'wt') as f: 
        json.dump(name, outfile)        
    else:
        this = parse(name)
  else: 
    file = r'../data/share/array/{}.json'.format(name)
    this = parse(file)

  this['ui'] = {'device':name,'family':family}
  return this

@functools.lru_cache(maxsize=32)
def get_MBF_XML_contours(xml_file):
  
  root = et.parse(xml_file).getroot()
  loop = []
  
  for c in root.findall('{*}contour'):
    
    nom = c.attrib['name'].lower()
    if 'blood' in nom: continue
    if 'outer' in nom: continue
    
    xy = [(float(p.attrib['x']),-float(p.attrib['y'])) for p in c.findall('{*}point')]

    xmin = min([x[0] for x in xy])
    xmax = max([x[0] for x in xy])

    loop.append({'name': c.attrib['name'], 'xy': xy, 'xr': (xmin,xmax) })
  
  return {"anat": loop} # return as JSONable dict


#%%
def get_Ve_matfile(filename,olddata=None):

    if filename is None: return olddata
    if isinstance(filename,list): filename = filename[0]
    if olddata is not None:
      if 'filename' in olddata and olddata['filename'] == filename: return olddata
    
    data = spio.loadmat(filename)
    
    output = list()
        
    obj_name = data['model'][0,0]['object_name']
    obj_name = [n[0] for n in obj_name[0]]
    
    obj_eidx = data['model'][0][0]['object_id']
    obj_eidx = [u for n,u in zip(obj_name,obj_eidx[0]) if "Fascicle" in n]
    obj_name = [n for n in obj_name if "Fascicle" in n]
    
    elem_indices = data['model'][0][0]['elems']
    
    for obj in obj_eidx:
        
        nn = elem_indices[obj-1,:]
        nn = np.unique(nn) - 1
        ve = data['v_extracellular'][nn,:]
        xyz = data['model'][0][0]['nodes'][nn,:]
    
        electrodes = ['elec{}'.format(e+1) for e in range(0,len(ve[0]))]    
        df = pd.DataFrame(data = xyz, columns = ['x','y','z'], index=None)        
        d2 = pd.DataFrame(data = ve, columns = electrodes, index=None)
        df = df.join(d2)
        
        output.append(df.to_dict())
    
    output = dict(zip(obj_name,output))
    output['filename'] = filename
    output['electrodes'] = electrodes
    output['type'] = 've'
    # print(output)
    return output

#%%
#
# Get menu options for device and deviceFamily lists 
# 
#%%


# get list of user-defined devices
def list_userSessions(user=1):

  session_folders = glob(r'../data/u/{}/*'.format(user))
  sessions = list()

  if not session_folders: sessions = [{"label":"Session 1","value":1}]
  for sid in [os.path.basename(p) for p in session_folders]:    
    sessions.append({"label":"Session {}".format(sid),"value":int(sid)})

  sessions.append({"label":"New Session","value":len(sessions)+1})

  return sessions




# get list of user-defined devices
def list_userDevices(user=1):
  my_list = glob(r'../data/u/{}/*/array.json'.format(user))

  dev_list = []

  if not my_list: return []
  for filepath in my_list:
    
    array = parse(filepath)

    if isinstance(array,list): array = array[0]
    try:
      dev_list.append({"label":array['array']['Name'],"value":filepath})
    except:
      print("error parsing user JSON array: "+filepath)
      print(array)

  return dev_list

# Load device families from /data/share
def list_deviceFamilies():
  
  arrays = parse(r'../data/share/array/index.json')

  dflist = set([a['family'] for a in arrays['list']])
  dflist = ([{"label":a,"value":a} for a in dflist])

  mylist = list_userDevices()
  if mylist: dflist.append({"label":"My devices","value":"user"})

  # todo append my-devices    
  return dflist

# Load devices given family from /data/share
def list_devices(family,user=1):
  if family is None: return []
  if family == 'user': 
    return list_userDevices(user)
  
  arrays = parse(r'../data/share/array/index.json')

  return [{"label":a['name'],"value":a['file']} 
          for a in arrays['list'] if a['family'] == family]

# Load nerve classes from /data/share
def list_nerveClasses():    
    axons = parse(r'../data/share/axon/index.json')
    return([{"label":a["label"],"value":a["value"]} for a in axons['list']])


def list_resultsFiles(user=1,session=1):

  # if not session == 1: print('LIST results-files: {}'.format(session))
  mat_files = glob('../data/u/{}/{}/*.mat'.format(user,session))
  print(mat_files)
  if not mat_files: return None  
  return([{"label":os.path.basename(p),"value":p} for p in mat_files])




#%%
#
# Make "complete" JSON to pass to oSPARC 
#
#%%

def make_NERVE_json(rot=0,dx=0,dy=0,sz=12,lc=0.1):
  
  s = {"nerve":{"xRotate":rot,   # degrees
                "xMove":[dx,dy], # µm
                "zRange":sz,     # mm
                "meshList":["Fascicle","Epineurium"],
                "MeshLengthFascicle":lc }}

  return s,json.dumps(s,indent=2)

def make_ARRAY_json(array):

  this = array['array']

  if isinstance(this['ElectrodeDimensions'][0],list):
        ew = this['ElectrodeDimensions'][0][1]
  else: ew = this['ElectrodeDimensions'][1]

  if this['InsetDepth'] < 0:
        min_tk = max(this['InsetDepth']+ew,0) + 0.1
  else: min_tk = 2*this['InsetDepth']+ew

  c = this['carrier']
  if "c_thickness" in this['carrier']:

    array['array']['carrier']['c_thickness'] = max(c['c_thickness'],min_tk)
    d_xyz = [ max(c['c_len']/2 + 1, 3), max(c['c_wid']/2 + 1, 3), 3 ]

  else:
    array['array']['carrier']['cuff_thickness'] = max(c['cuff_thickness'],min_tk)
    min_tk = array['array']['carrier']['cuff_thickness']

    d_xyz = [ max(c['cuff_length']+1, 6), max(c['cuff_IDy'] + min_tk + 1, 3), 
                                          max(c['cuff_IDx'] + min_tk + 1, 3)]
  
  array['mesh']['DomainSize'] = [max(a,b) for a,b in zip(array['mesh']['DomainSize'],d_xyz)]
  array['mesh']['MeshLengthMax'] =  min(array['mesh']['MeshLengthMax'], 0.2)
  array.pop("ui",[]) # remove 'ui' from dict 

  return array,json.dumps(array,indent=2)




def check_existing_nerve_files(user=1,session=1):

  filename = '../data/u/{}/{}/nerve.xml'.format(user,session)
  if not os.path.exists(filename): return None,None

  anat = get_MBF_XML_contours(filename)

  filename = '../data/u/{}/{}/nerve.json'.format(user,session)
  if not os.path.exists(filename): 
    return anat,make_NERVE_json() # default values

  config = parse(filename)

  return anat,config 




def get_user_XML_contours(user=1,session=1):

  path = '../data/u/{}/{}/nerve.xml'.format(user,session)
  return get_MBF_XML_contours(path)



def has_results(user=1,session=1):

  file_list = list_resultsFiles(user,session)
  if not file_list: return False
  else: return True


def get_results_file(filename):

  if filename is None: return None
  if "nerve" in filename: 
    print('load NERVE file not implemented yet')
  if "v-extra" in filename:
    return get_Ve_matfile(filename)



def save_json_files(array,nerve,axon_select,session=1):


    path = '../data/u/{}/{}/array.json'.format(get_user_ID(),session)
    array,json_string = user_files.make_ARRAY_json(array)

    print('Saving '+path)
    with open(path,'wt') as f:
      f.write(json_string)

    r = nerve['nerve']['xRotate']
    m = nerve['nerve']['xMove']
    z = array['mesh']['DomainSize']

    n2,json_string = user_files.make_NERVE_json(r,m[0],m[1],z[0],lc=0.1)

    n2['nerve']['uifileName'] = nerve['nerve']['source']
    n2['nerve']['uiAxonPop']  = axon_select

    json_string = json.dumps(nerve,indent=2)
    path = '../data/u/{}/{}/nerve.json'.format(get_user_ID(),session)
    print('Saving '+path)
    with open(path,'wt') as f:
      f.write(json_string)

    return json_string # dict(content=json_string, filename="nerve.json")




