# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 06:38:15 2021

@author: Calvin
"""

import io
import math
import json
import base64

from numpy import linspace
import matplotlib.pyplot as plt
import xml.etree.ElementTree as et
from plotly.tools import mpl_to_plotly

from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from statistics import median

import user_files

def parse(text):
    try:
        return json.load(text)
    except ValueError as e:
        print('invalid json: %s' % e)
        return None # or: raise

def encode(svg_string):
    encoded = base64.b64encode(svg_string) 
    return 'data:image/svg+xml;base64,{}'.format(encoded.decode()) 

def gcf_to_svg(): # convert active figure (gcf) to svg

  # render to SVG
  f = io.BytesIO()
  plt.savefig(f, format = "svg",pad_inches=0)
  return f.getvalue() # svg string

# Electrode view 1 (xy)
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

# Electrode view 2 (yz)
def elec_yz(d, eid):


    print('TODO elec_yz')



def rounded_rect(w,h,r):
        
    a = linspace(0,math.pi/2,16)
    xc = w/2 - r
    yc = w/2 - r

    x =      [ r*math.cos(a)+xc for a in a]
    x.extend([-r*math.sin(a)-xc for a in a])
    x.extend([-r*math.cos(a)-xc for a in a])
    x.extend([ r*math.sin(a)+xc for a in a])
    x.append(x[0])
    y =      [ r*math.sin(a)+yc for a in a]
    y.extend([ r*math.cos(a)+yc for a in a])
    y.extend([-r*math.sin(a)-yc for a in a])
    y.extend([-r*math.cos(a)-yc for a in a])
    y.append(y[0])

    return x,y


# Nerve + cross-section view

def find_centroid(c):
    
    x0 = []
    y0 = []
    n = 0
    for u in range(0,len(c)):
        x0.append(sum([xy[0] for xy in c[u]['xy']]))
        y0.append(sum([xy[1] for xy in c[u]['xy']]))
        n = n+len(c[u]['xy'])
        
    return (sum(x0)/n, sum(y0)/n)

def apply_xForm(anat,xf)

  if xf is None: return anat
  if 'nerve' in xf: xf = xf['nerve']

  xy0 = find_centroid(anat['anat']) # find xy0, needed for rotate

  print(anat)
  print(xf)

  for u in range(0,len(anat['anat'])):

    x = [xy[0] for xy in anat['anat'][u]['xy']]
    y = [xy[1] for xy in anat['anat'][u]['xy']]

    if 'xRotate' in xf:
      a = xf['xRotate']/180*math.pi
      x = [math.cos(a)*(xy[0]-xy0[0]) - 
           math.sin(a)*(xy[1]-xy0[1]) + xy0[0] for xy in c[u]['xy']]
      y = [math.sin(a)*(xy[0]-xy0[0]) + 
           math.cos(a)*(xy[1]-xy0[1]) + xy0[1] for xy in c[u]['xy']]

    if 'xMove' in xf:
      x = [x+xf['xMove'][0] for x in x]
      y = [y+xf['xMove'][1] for y in y]

    xmin = min([x(0) for x in xy])
    xmax = max([x(1) for x in xy])

    anat['anat'][u]['xy'] = [(x,y) for x,y in zip(x,y)]
    anat['anat'][u]['xr'] = (xmin,xmax)

def mk_ARRAY_subplot(array,nerve,ax):
    
    if array is None: return None
    if 'array' in array: array = array['array']

    ax.set_aspect('equal', 'box')

    # make outline of PDMS carrier

    if 'carrier' in array:
      if 'c_len' in array['carrier']: # flat
        if 'c_radius' in array['carrier']: cr = array['carrier']['c_radius']
        else:                              cr = 0.3 # default value defined in +mesh.insert_gmsh_electrodes
        
        x,y = rounded_rect(array['carrier']['c_len'], 
                           array['carrier']['c_wid'], cr)

        ax.plot(x,y,'-',color='k',linewidth=0.9)
      else:

        x1 = array['carrier']['cuff_IDx']/2
        x2 = array['carrier']['cuff_thickness']
        y1 = array['carrier']['cuff_length']/2

        x = [x1+a*x2 for a in [0,1,1,0,0]]
        y = [ y1*b   for b in [1,1,-1,-1,1]]

        ax.fill(y,x,'#ccc', hatch='//',edgecolor='k',linewidth=1.15)

        x = [-x for x in x]        
        ax.fill(y,x,'#ccc', hatch='//',edgecolor='k',linewidth=1.15)

        x = [(x1+x2)*a for a in [1,-1,-1,1,1]]
        ax.plot(y,x,'-',color='k',linewidth=0.9)

    # show each electrode
    for u in range(0,len(array['ElectrodeTypeIndex'])):
        
        x,y,s = elec_xy(array,u) # get electrode
        h = ax.plot(y,x,s)
        
        if 'ElectrodeAngle' in array:
            if s == '--': va = 'top'
            else:         va = 'bottom'
        else:             va = 'center'
        
        ax.text((y[0]+y[2])/2, (x[0]+x[2])/2, 'E%d'%(u+1), 
                     color=h[0].get_color(),
                        ha='center', va=va, 
                        fontweight='bold',fontsize=14)

    if nerve is not None:

        print('TODO: show nerve outline')
        print(nerve)

def mk_NERVE_subplot(array,nerve,ax):





    print('TODO')
    print(array)
    print(nerve)
    
    ax.set_aspect('equal', 'box')
    
    
    for u in range(0,len(c)):
        
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



def view_model_inputs(array,nerve):

    if array is None and nerve is None: raise PreventUpdate

    if not isinstance(array, dict ): 
      filename = array
      with open(filename) as f: 
        array = parse(f)
      
    if not isinstance(nerve, dict ): 
      filename = array
      with open(filename) as f: 
        array = parse(f)
    
      
    if 'array' in array: array = array['array'] # get rid of {mesh}    
    assert('ElectrodeDimensions' in array) # QC
      
    # make matplotlib axis and fill


    fig, ax = plt.subplots(2,1,num=1, clear=True,'sharex'=True)
    fig.set_size_inches(8.0, 4.0)

    ax = subplot(2, 1, 1, figure=fig)
    mk_ARRAY_subplot(array,nerve,ax[0])

    ax = subplot(2, 1, 1, figure=fig)
    mk_NERVE_subplot(array,nerve,ax[1])

    




def add_callbacks(app):


  # Simple update to electrode drop-down
  @app.callback(Output("view-device","src"), 
                 Input("device-json","data"), 
                 Input("nerve-json","data"), 
                 Input("anatomy-json","data"))
  def update_array(array,nerve,anatomy):
    if array is None: 
      raise PreventUpdate

    if isinstance(array,list) and array: array = array[0]
    if isinstance(nerve,list) and nerve: nerve = nerve[0]

    nerve = apply_xForm(anatomy,nerve)


    view_model_inputs(array,nerve) # make mpl figure

    return encode(gcf_to_svg())


