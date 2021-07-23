# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 06:38:15 2021

@author: Calvin
"""

import json
import matplotlib.pyplot as plt
import xml.etree.ElementTree as et
from numpy import linspace
from plotly.tools import mpl_to_plotly
import math
import io
import base64

from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate


from statistics import median

def parse(text):
    try:
        return json.load(text)
    except ValueError as e:
        print('invalid json: %s' % e)
        return None # or: raise

def encode(svg_string):
    encoded = base64.b64encode(svg_string) 
    return 'data:image/svg+xml;base64,{}'.format(encoded.decode()) 

# Electrode view 1  
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
            z = [epos[1]+elwh[1]*u/2 for u in [1,-1,-1,1,1]]
            x = [x*math.cos(a)+z*math.sin(a) for x,z in zip(x,z)]
            l = '--'
    
    return x, y, l

def array_SVG(filename,nerve_json=None):
    
    if filename is None:
        return b''

    if isinstance(filename, dict ): 
      array = filename
    else: 
      with open(filename) as f: 
        array = parse(f)
      
      
    array = array['array'] # get rid of {mesh}    
    assert('ElectrodeDimensions' in array) # QC
      
    # make matplotlib axis and fill
    fig, ax = plt.subplots()  # a figure with a single Axes
    ax.set_aspect('equal', 'box')

    if 'carrier' in array:
      if 'c_len' in array['carrier']: # flat
        if 'c_radius' in array['carrier']: cr = array['carrier']['c_radius']
        else:                              cr = 0.3 # default value defined in +mesh.insert_gmsh_electrodes
        a = linspace(0,math.pi/2,16)
        xc = array['carrier']['c_len']/2-cr
        yc = array['carrier']['c_wid']/2-cr

        x = [cr*math.cos(a)+xc for a in a]
        x.extend([-cr*math.sin(a)-xc for a in a])
        x.extend([-cr*math.cos(a)-xc for a in a])
        x.extend([cr*math.sin(a)+xc for a in a])
        x.append(x[0])
        y = [cr*math.sin(a)+yc for a in a]
        y.extend([cr*math.cos(a)+yc for a in a])
        y.extend([-cr*math.sin(a)-yc for a in a])
        y.extend([-cr*math.cos(a)-yc for a in a])
        y.append(y[0])

        plt.plot(x,y,'-',color='k',linewidth=0.9)
      else:

        x1 = array['carrier']['cuff_IDx']/2
        x2 = array['carrier']['cuff_thickness']
        y1 = array['carrier']['cuff_length']/2

        x = [x1+a*x2 for a in [0,1,1,0,0]]
        y = [ y1*b    for b in [1,1,-1,-1,1]]

        ax.fill(y,x,'#ccc', hatch='//',edgecolor='k',linewidth=1.15)

        x = [-x for x in x]        
        ax.fill(y,x,'#ccc', hatch='//',edgecolor='k',linewidth=1.15)

        x = [(x1+x2)*a for a in [1,-1,-1,1,1]]
        plt.plot(y,x,'-',color='k',linewidth=0.9)


    for u in range(0,len(array['ElectrodeTypeIndex'])):
        
        x,y,s = elec_xy(array,u) # get electrode
        h = plt.plot(y,x,s)
        
        if 'ElectrodeAngle' in array:
            if s == '--': va = 'top'
            else:         va = 'bottom'
        else:             va = 'center'
        
        plt.text((y[0]+y[2])/2, (x[0]+x[2])/2, 'E%d'%(u+1), 
                     color=h[0].get_color(),
                        ha='center', va=va, 
                        fontweight='bold',fontsize=14)

    if nerve_json is not None:
        print('TODO: show nerve outline')
        

    # render to SVG
    f = io.BytesIO()
    plt.savefig(f, format = "svg")
    return f.getvalue() # svg string

# Nerve + cross-section view

def get_contours(xml_file):
    
    root = et.parse(xml_file).getroot()
    loop = []
    
    for c in root.findall('{*}contour'):
        
        nom = c.attrib['name'].lower()
        if 'blood' in nom: continue
        if 'outer' in nom: continue
        
        xy = [(float(p.attrib['x']),-float(p.attrib['y'])) for p in c.findall('{*}point')]

        loop.append({'name': c.attrib['name'], 
                     'xy': xy })
    return loop
    
def find_centroid(c):
    
    x0 = []
    y0 = []
    n = 0
    for u in range(0,len(c)):
        x0.append(sum([xy[0] for xy in c[u]['xy']]))
        y0.append(sum([xy[1] for xy in c[u]['xy']]))
        n = n+len(c[u]['xy'])
        
    return (sum(x0)/n, sum(y0)/n)


def nerve_SVG(xml_file, json_file = None):

    if xml_file is None:
        return b''

    c = get_contours(xml_file) # load data
    
    if json_file is not None:
        print('loading' + json_file)
        with open(json_file) as f:
            xf = parse(f)
        print(xf.__class__)
        xf = xf['nerve']
        print(xf.__class__)
        
        xy0 = find_centroid(c) # find xy0, needed for rotate

    # make matplotlib axis and fill
    fig, ax = plt.subplots()  # a figure with a single Axes
    ax.set_aspect('equal', 'box')
    # plt.axis('off')
    
    
    for u in range(0,len(c)):
        
        x = [xy[0]/1e3 for xy in c[u]['xy']]
        y = [xy[1]/1e3 for xy in c[u]['xy']]

        if json_file is not None:
          if 'xRotate' in xf:
            a = xf['xRotate']/180*math.pi
            x = [math.cos(a)*(xy[0]-xy0[0])/1e3 - 
                 math.sin(a)*(xy[1]-xy0[1])/1e3 + xy0[0]/1e3 for xy in c[u]['xy']]
            y = [math.sin(a)*(xy[0]-xy0[0])/1e3 + 
                 math.cos(a)*(xy[1]-xy0[1])/1e3 + xy0[1]/1e3 for xy in c[u]['xy']]

          if 'xMove' in xf:
            x = [x+xf['xMove'][0]/1e3 for x in x]
            y = [y+xf['xMove'][1]/1e3 for y in y]

        x.append(x[0])
        y.append(y[0])
        
        h = plt.plot(x,y,'-')    
        
        if 'inner' in c[u]['name'].lower():
            
            plt.text(median(x), median(y), 'F%d'%(u+1), 
                         color=h[0].get_color(),
                            ha='center',fontweight='bold',fontsize=14)
        else:
            h[0].set_color('#bbbbbb')

    plt.plot(xy0[0]/1e3,xy0[1]/1e3,'r+')

    # render to SVG
    f = io.BytesIO()
    plt.savefig(f, format = "svg")
    return f.getvalue() # svg string


if __name__ == "__main__":
  
    uid = 1 # user ID
    sid = 1 # session ID 
      
    # xml_file = 's3://pennsieve-prod-discover-publish-use1/65/6/files/derivative/sub-47/sam-2/sub-47_sam-2_C47-2MergeMask.xml'
    # with open("nerve.svg",'wb') as f:
    #     f.write(nerve_SVG(xml_file))

    xml_file = r'..\data\u\{0}\{1}\nerve.xml'.format(uid,sid)
    json_file = r'..\data\u\{0}\{1}\nerve.json'.format(uid,sid)

    print( nerve_SVG(xml_file, json_file) )
    
    
    # print(array_SVG('misc/C-FINE.json'))








def add_callbacks(app):


  # Simple update to electrode drop-down
  @app.callback(Output("view-device","src"), 
                 Input("device-json","data"), 
                 Input("nerve-json","data"), 
                 )
  def update_array(device,nerve):
    if device is None: 
      raise PreventUpdate

    if isinstance(device,list) and device: device = device[0]
    return encode(array_SVG(device))



  print('TODO add graphics (nerve) callback')


