import plotly.plotly as py
import plotly.graph_objs as go

import numpy as np

from numpy import genfromtxt
recHitPos = genfromtxt('test_All_50.csv', delimiter=',')
x = []
y = []
z = []
for i in range(0, len(recHitPos)):
  x.append(recHitPos[i][0])
  y.append(recHitPos[i][1])
  z.append(recHitPos[i][2])

trace1 = go.Scatter3d(
    x=x,
    y=y,
    z=z,
    mode='markers',
    marker=dict(
        size=6,
        line=dict(
            color='rgba(217, 217, 217, 0.14)',
            width=0.5
        ),
        opacity=0.8
    )
)

data = [trace1]
layout = go.Layout(
    scene=go.Scene(
        xaxis=dict(title='Time [ns]'),
        yaxis=dict(title='z-coordinate [cm]'),
        zaxis=dict(title='rho [cm]')
    ),
    margin=dict(
        l=0,
        r=0,
        b=0,
        t=0
    )
)
fig = go.Figure(data=data, layout=layout)
py.plot(fig, filename='TimeRhoZClustering')
