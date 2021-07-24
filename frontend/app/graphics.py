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
            x = [x*math.cos(a)-z*math.sin(a) for x,z in zip(x,z)]
            l = '--'
    
    return x, y, l

# Electrode view 2 (yz)
def elec_yz(d, eid):

  etype = d['ElectrodeTypeIndex'][eid]-1
  epos = d['ElectrodePositions'][eid]
  elwh = d['ElectrodeDimensions']
  if 'ElectrodeAngle' in d:
        erot = d['ElectrodeAngle'][eid]*math.pi/180
  else: erot = 0
  
  if isinstance(elwh[0],list):
    elwh = elwh[etype]

  if 'ElectrodeKind' in d: 
    if isinstance(d['ElectrodeKind'],list):
      if etype < len(d['ElectrodeKind']):
        ekind = d['ElectrodeKind'][etype]
      else:
        ekind = d['ElectrodeKind'][0]
    else: # not a list
      ekind = d['ElectrodeKind']
  else: ekind = 'rect'
  
  if ekind.startswith('circum'):
    
     radius = (d['carrier']['cuff_IDx']/2+d['InsetDepth'])  
     
     a = linspace(-elwh[2]/radius/2, 
                   elwh[2]/radius/2, 33).tolist()
     
     x = [radius*math.sin(a-erot) for a in a]
     y = [radius*math.cos(a-erot) for a in a]
     
     radius = radius + elwh[1]
     
     x.extend([radius*math.sin(-a-erot) for a in a])
     y.extend([radius*math.cos(-a-erot) for a in a])
     
  else:    
      
    x0 = [epos[2]+a*elwh[2]/2 for a in [-1,-1,1,1]]
    y0 = [epos[1]-d['InsetDepth']-a*elwh[1] for a in [0,1,0,1]]
    
    x = [x*math.cos(erot)-y*math.sin(erot) for x,y in zip(x0,y0)]
    y = [y*math.cos(erot)+x*math.sin(erot) for x,y in zip(x0,y0)]
    
    x = [x[u] for u in [0,1,3,2,0]]
    y = [y[u] for u in [0,1,3,2,0]]
      
  return x,y


def rounded_rect(w,h,r):
        
    a = linspace(0,math.pi/2,16).tolist()
    xc = w/2 - r
    yc = h/2 - r

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

def apply_xForm(anat,xf):

  if xf is None or anat is None: return anat
  if 'nerve' in xf: xf = xf['nerve']

  try:
    xy0 = find_centroid(anat['anat']) # find xy0, needed for rotate
  except:

    print('BEGIN DEBUG apply_xForm info')
    print(anat)
    print(xf)
    print('END DEBUG apply_xForm info')


  for u in range(0,len(anat['anat'])):

    xy = anat['anat'][u]['xy']
    

    if 'xRotate' in xf:
      a = xf['xRotate']/180*math.pi
      x = [math.cos(a)*(xy[0]-xy0[0]) - 
           math.sin(a)*(xy[1]-xy0[1]) + xy0[0] for xy in xy]
      y = [math.sin(a)*(xy[0]-xy0[0]) + 
           math.cos(a)*(xy[1]-xy0[1]) + xy0[1] for xy in xy]
    else: 
      x = [u[0] for u in xy]
      y = [u[1] for u in xy]

    if 'xMove' in xf:
      x = [x+xf['xMove'][0] for x in x]
      y = [y+xf['xMove'][1] for y in y]


    anat['anat'][u]['xy'] = [(x,y) for x,y in zip(x,y)]
    anat['anat'][u]['xr'] = (min(x),max(x))

  if 'zRange' in xf: 
    if 'info' in anat: 
          anat['info']['zRange'] = xf['zRange']
    else: anat['info'] = {'zRange':xf['zRange']}
    
  return anat

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
            
      zr = ax.get_xlim()
      x = [zr[u] for u in [0,1,1,0,0]]
      
      inner = [obj for obj in nerve['anat'] if 'inner' in obj['name'].lower()]
      outer = [obj for obj in nerve['anat'] if 'inner' not in obj['name'].lower()]
      
      for obj in outer:
        y = [obj['xr'][u]/1e3 for u in [0,0,1,1,0]]
        ax.fill(x,y,'#333', '-',alpha=0.2,edgecolor='#333')        
      
      ax.set_prop_cycle(None) # reset color sequence
      
      for fid in range(0,len(inner)):      
        y = [inner[fid]['xr'][u]/1e3 for u in [0,0,1,1,0]]
        h = ax.fill(x,y,alpha=0.3)
        
        ax.text(x[1]+0.1, y[0]/2+y[2]/2, 'F%d'%(fid+1), 
                            color=h[0].get_fc(),
                            fontsize=14,alpha=1.0)
        
        
        
        
def mk_NERVE_subplot(array,nerve,ax):

    ax.set_aspect('equal', 'box')        
    if 'array' in array: array = array['array']
    
    if 'carrier' in array: 
      c = array['carrier']
      if 'c_len' in c: # flat

        x = [0,0,-c['c_thickness'],-c['c_thickness'],0]
        y = [   a*c['c_wid']/2  for a in [1,-1,-1,1,1]]
        ax.fill(y,x,'#ccc', hatch='//',edgecolor='k',linewidth=1.15)

      else:

        x,y = rounded_rect(c['cuff_IDx']+2*c['cuff_thickness'], 
                           c['cuff_IDy']+2*c['cuff_thickness'], 
                           c['cuff_IDr']+c['cuff_thickness'])
        
        ax.fill(x,y,'#ccc', hatch='//',edgecolor='k',linewidth=1.15)  
        
        x,y = rounded_rect(c['cuff_IDx'], c['cuff_IDy'], c['cuff_IDr'])
        ax.fill(x,y,'w',edgecolor='k',linewidth=1.15)

    
    for u in range(0,len(array['ElectrodeTypeIndex'])):
      x,y = elec_yz(array,u) # get electrode
      ax.fill(x,y,linewidth=1.15)
    
    inner = [obj for obj in nerve['anat'] if 'inner' in obj['name'].lower()]
    outer = [obj for obj in nerve['anat'] if 'inner' not in obj['name'].lower()]
    
    for obj in outer:
        
      x = [x[0]/1e3 for x in obj['xy']]
      y = [x[1]/1e3 for x in obj['xy']]        
      ax.fill(x,y,'#333', '-',alpha=0.2,edgecolor='#333')

    ax.set_prop_cycle(None) # reset color sequence

    for u in range(0,len(inner)):
        
      x = [x[0]/1e3 for x in inner[u]['xy']]
      y = [x[1]/1e3 for x in inner[u]['xy']]
        
      h = ax.plot(x,y,'-')    
        
      ax.text(median(x), median(y), 'F%d'%(u+1), 
                          color=h[0].get_color(),
                          ha='center',fontweight='bold',fontsize=14)

    xl = ax.get_xlim()
    yl = ax.get_ylim()
    if xl[1] < 1: 
      xl = [xl[1],xl[1]-0.1]
      lbl = '100 µm'
    else:
      xl = [xl[1],xl[1]-1]
      lbl = '1 mm'

    ax.plot(xl,[yl[0],yl[0]],'-',color='#333',linewidth=3)
    ax.text(xl[0]/2+xl[1]/2, yl[0]*1.05-yl[1]*0.05, lbl, 
                                color='#333',ha='center',va='top' )


def view_model_inputs(array,nerve,xform=None):

    if array is None and nerve is None: raise PreventUpdate

    if not isinstance(array, dict): 
      filename = array
      with open(filename) as f: 
        array = parse(f)
    
    if nerve is None:
      nerve = {'anat':[]}
    elif not isinstance(nerve, dict): 
      filename = nerve
      with open(filename) as f: 
        nerve = parse(f)
        
    if xform is not None:
      if not isinstance(xform, dict):
        filename = xform
        with open(filename) as f: 
          xform = parse(f)
      nerve = apply_xForm(nerve,xform)
      
    if 'array' in array: array = array['array'] # get rid of {mesh}    
      
    # make matplotlib figure axes
    fig, ax = plt.subplots(2,1,num=1, clear=True, 
                               sharex=True, figsize=(4.8,9.2))
    
    mk_ARRAY_subplot(array,nerve,ax[0])
    mk_NERVE_subplot(array,nerve,ax[1])

    ax[0].set_axis_off()
    ax[1].set_axis_off()

    fig.tight_layout()



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



if __name__ == "__main__":
    
    
    array = "../data/u/1/1/array.json"
    nerve = "../data/u/1/1/nerve.xml"
    xform = "../data/u/1/1/nerve.json"
    
    nerve = user_files.get_MBF_XML_contours(nerve)
    
    view_model_inputs(array,nerve,xform)
  
    # x,y = elec_yz(array,u) # get electrode



    
    