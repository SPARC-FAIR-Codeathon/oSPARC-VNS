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
    
    
    
    
    
xml_file = r'C:\Users\Calvin\Documents\MATLAB\Keast-lab\oSPARC-VNS\backend\module-mesher\input\sub-57_sam-1.xml'

c = get_contours(xml_file)

# make matplotlib axis and fill
fig, ax = plt.subplots()  # a figure with a single Axes
ax.set_aspect('equal', 'box')

for u in range(0,len(c)):
    
    x = [xy[0]/1e3 for xy in c[u]['xy']]
    y = [xy[1]/1e3 for xy in c[u]['xy']]
    
    x.append(x[0])
    y.append(y[0])
    
    h = plt.plot(x,y,'-')    
    
    if 'inner' in c[u]['name'].lower():
        
        plt.text(median(x), median(y), 'F%d'%(u+1), 
                     color=h[0].get_color(),
                        ha='center',fontweight='bold',fontsize=14)
    else:
        h[0].set_color('#bbbbbb')

# if __name__ == "__main__":
#   
#    r = jsonify(r'C:\Users\Calvin\Documents\MATLAB\Keast-lab\oSPARC-VNS\backend\module-mesher\input\sub-57_sam-1.xml')


'''svg
<svg width="750" height="500" style="background: gray">
<svg x="100" y="100" style="fill: yellow; stroke: red">
<rect x="0" y="0" width="142" height="142" />
<rect x="100" y="-50" width="100" height="100" style="transform: rotate(45deg)" />
</svg>
</svg>
'''


'''html

<div class='draggable'></div>

'''


'''css 

.draggable {
  position: absolute;
  background: #EEE;
  border: 1px solid #AAA;  
  top: 50%;
  left: 50%;
  width: 300px;
  height: 200px;
  cursor: move;
}

.resize-handle {
  position: absolute;
  width: 0px;
  height: 0px;
  margin-top: -17px;
  margin-left: -17px;
  border-style: solid;
  border-width: 0 0 16px 16px;
  border-color: transparent transparent #AAA transparent;
  cursor: grab;
  cursor: -moz-grab;
  cursor: -webkit-grab;
}

.resize-handle:active {
  cursor: grabbing;
  cursor: -moz-grabbing;
  cursor: -webkit-grabbing;
}

'''

'''.js

var drag = $(".draggable");
var handle = $("<div class='resize-handle'></div>").appendTo(drag);
TweenLite.set(drag, { xPercent: -50, yPercent: -50 });
TweenLite.set(handle, { x: drag.width(), y: drag.height()});

// Start rotation dragging
$(handle).on("mousedown", function(event) {

  event.stopPropagation(); // cancel x,y drag
  Draggable.create(drag, { type: "rotation" })[0].startDrag(event);
});

// Start x,y dragging
$(drag).on("mousedown", function(event) {

  Draggable.create(drag)[0].startDrag(event);
});

''' 