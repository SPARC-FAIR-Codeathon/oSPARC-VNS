
import dash
import dash_core_components as dcc

# layout_preview = dcc.Graph(figure=px.scatter(width=480, height=400))

import json


# The local store will take the initial data
# only the first time the page is loaded
# and keep it until it is cleared.
dcc.Store(id='array', storage_type='local'),
dcc.Store(id='nerve', storage_type='local'),
dcc.Store(id='config', storage_type='local'),


with open('path_to_file/person.json') as f:
  data = json.load(f)


person = '{"name": "Bob", "languages": ["English", "Fench"]}'
person_dict = json.loads(person)




import matplotlib.pyplot as plt
import numpy as np
import io

f = io.BytesIO()
a = np.random.rand(10)
plt.bar(range(len(a)), a)
plt.savefig(f, format = "svg")

print(f.getvalue()) # svg string
