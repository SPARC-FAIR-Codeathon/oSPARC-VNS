# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 06:38:15 2021

@author: Calvin
"""

import json
import matplotlib.pyplot as plt
import xml.etree.ElementTree as et
import math
import io

from statistics import median


def get_contours(xml_file):
    
    root = et.parse(xml_file).getroot()
    loop = []
    
    for c in root.findall('{*}contour'):
        
        nom = c.attrib['name'].lower()        
        if 'blood' in nom: continue
        if 'outer' in nom: continue
        
        xy = [(float(p.attrib['x']),float(p.attrib['y'])) for p in c.findall('{*}point')]

        loop.append({'name': c.attrib['name'], 
                     'xy': xy })
    return loop
    
    
def nerve_SVG(xml_file, json_xform = None):

    c = get_contours(xml_file) # load data 

    # make matplotlib axis and fill
    fig, ax = plt.subplots()  # a figure with a single Axes
    ax.set_aspect('equal', 'box')

    for u in range(0,len(c)):
        
        x = [xy[0]/1e3 for xy in c[u]['xy']]
        y = [xy[1]/1e3 for xy in c[u]['xy']]

        if json_xform is not None:
          xf = json_xform['nerve']
          if 'xRotate' in xf
            a = xf['xRotate']/180*math.pi
            x = [math.cos(a)*xy[0]/1e3 - math.sin(a)*xy[1]/1e3 for xy in c[u]['xy']]
            y = [math.sin(a)*xy[0]/1e3 + math.cos(a)*xy[1]/1e3 for xy in c[u]['xy']]

          if 'xMove' in xf:
            x = [x+xf['xMove'][0] for x in x]
            y = [y+xf['xMove'][1] for y in y]

        x.append(x[0])
        y.append(y[0])
        
        h = plt.plot(x,y,'-')    
        
        if 'inner' in c[u]['name'].lower():
            
            plt.text(median(x), median(y), 'F%d'%(u+1), 
                         color=h[0].get_color(),
                            ha='center',fontweight='bold',fontsize=14)
        else:
            h[0].set_color('#bbbbbb')

    # render to SVG
    f = io.BytesIO()
    plt.savefig(f, format = "svg")
    return f.getvalue() # svg string


if __name__ == "__main__":
  
    # xml_file = r'C:\Users\Calvin\Documents\MATLAB\Keast-lab\oSPARC-VNS\backend\module-mesher\input\sub-57_sam-1.xml'
    xml_file = 's3://pennsieve-prod-discover-publish-use1/65/6/files/derivative/sub-47/sam-2/sub-47_sam-2_C47-2MergeMask.xml'
    with open("nerve.svg",'wb') as f:
        f.write(nerve_SVG(xml_file))
    
