import plotly.graph_objects as go
import numpy as np

DATA = np.loadtxt(open("surface.dat", "rb"))

Xs = DATA[:,0]
Ys = DATA[:,1]
Zs = DATA[:,2]

value = np.loadtxt(open("0.dat", "rb"))

fig = go.Figure(data=[
    go.Mesh3d(
        x=Xs,
        y=Ys,
        z=Zs,
        colorbar_title='value',
        colorscale="hot",
        # Intensity of each vertex, which will be interpolated and color-coded
        intensity=value,
        showscale=True
    )
])

fig.show()