# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 02:35:52 2021

@author: Calvin
"""

import json



def mk_nerve(rot=0,dx=0,dy=0,sz=12):
    
    s = {"nerve":{"xRotate":rot, 
                  "xMove":[dx,dy],
                  "zRange":sz,
                  "meshList":["Fascicle","Epineurium"],
                  "MeshLengthFascicle":0.1 }}

    return s,json.dumps(s,indent=2)




def mk_array(array):

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




def check_existing_nerve_files():




  print('todo: check existing nerve files / json')