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
                  "MeshLengthFascicle":0.1               }}

    return json.dumps(s)




def mk_array(rot=0,dx=0,dy=0,sz=12):
    
    s = {"nerve":{"xRotate":rot, 
                  "xMove":[dx,dy],
                  "zRange":sz,
                  "meshList":["Fascicle","Epineurium"],
                  "MeshLengthFascicle":0.1               }}

    return json.dumps(s)
