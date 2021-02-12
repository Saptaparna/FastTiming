import plotly.express as px

import chart_studio.plotly as py
import plotly.graph_objs as go
import plotly.graph_objs.layout as gos
import pandas as pd
import chart_studio
chart_studio.tools.set_credentials_file(username='Sapta', api_key='l8vqYeEid8uR8xiRQ7in')

import numpy as np

from numpy import genfromtxt

recHitPos = genfromtxt('info.txt', delimiter=',')
x = []
y = []
z = []
t = []
for i in range(0, len(recHitPos)):
  x.append(recHitPos[i][0])
  y.append(recHitPos[i][1])
  z.append(recHitPos[i][2])
  t.append(recHitPos[i][3])



#df = px.data.gapminder()
df = pd.read_csv('info.txt')
df = df.sort_values(by=['t'], ascending=True)
print (df)
fig = px.scatter_3d(df, x="x", y="y", z="z", animation_frame="t",
           log_x=False, size_max=5, opacity=0.7, range_x=[-300,300], range_y=[-300,300])

fig["layout"].pop("updatemenus") # optional, drop animation buttons
fig.show()
