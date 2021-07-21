# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import json
import matplotlib.pyplot as plt
import math
import io
  
def parse(text):
    try:
        return json.load(text)
    except ValueError as e:
        print('invalid json: %s' % e)
        return None # or: raise
    
def elec_xy(d, eid):    
    
    etype = d['ElectrodeTypeIndex'][eid]-1    
    epos = d['ElectrodePositions'][eid]
    elwh = d['ElectrodeDimensions']
    
    if isinstance(elwh[0],list):
        elwh = elwh[etype]
    
    x = [epos[2]+elwh[2]*u/2 for u in [1,-1,-1,1,1]]
    y = [epos[0]+elwh[0]*u/2 for u in [1,1,-1,-1,1]]
    l = '-'
    
    if 'ElectrodeAngle' in d:
        a = d['ElectrodeAngle'][eid]*math.pi/180
        if a:
            z = [epos[0]+elwh[0]*u/2 for u in [1,1,-1,-1,1]]
            x = [x*math.cos(a)+z*math.sin(a) for x,z in zip(x,z)]
            l = '--'
    
    return x, y, l

def array_SVG(filename,nerve_json=None):
    
    with open(filename) as f:
      array = parse(f)
    array = array['array'] # get rid of {mesh}    
    assert('ElectrodeDimensions' in array) # QC
      
    # make matplotlib axis and fill
    fig, ax = plt.subplots()  # a figure with a single Axes

    if 'carrier' in array:        
        print('TODO: show carrier outline')


    for u in range(0,len(array['ElectrodeTypeIndex'])):
        
        x,y,s = elec_xy(array,u) # get electrode
        h = plt.plot(x,y,s)    
        
        if 'ElectrodeAngle' in array:
            if s == '--': va = 'top'
            else:         va = 'bottom'
        else:             va = 'center'    
        
        plt.text((x[0]+x[1])/2, (y[0]+y[2])/2, 'E%d'%(u+1), 
                     color=h[0].get_color(),
                        ha='center', va=va, 
                        fontweight='bold',fontsize=14)

    if nerve_json is not None:
        print('TODO: show nerve outline')
        

    # render to SVG
    f = io.BytesIO()
    plt.savefig(f, format = "svg")
    return f.getvalue() # svg string


if __name__ == '__main__':
    
    print(array_SVG('misc/C-FINE.json'))