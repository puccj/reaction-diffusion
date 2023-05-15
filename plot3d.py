import plotly.graph_objects as go
import numpy as np


# DATA = np.array([
#     [-0.807237702464, 0.904373229492, 111.428744443],
#     [-0.802470821517, 0.832159465335, 98.572957317],
#     [-0.801052795982, 0.744231916692, 86.485869328],
#     [-0.802505546206, 0.642324228721, 75.279804677],
#     [-0.804158144115, 0.52882485495, 65.112895758],
#     [-0.806418040943, 0.405733109371, 56.1627277595],
#     [-0.808515314192, 0.275100227689, 48.508994388],
#     [-0.809879521648, 0.139140394575, 42.1027499025],
#     [-0.810645106092, -7.48279012695e-06, 36.8668106345],
#     [-0.810676720161, -0.139773175337, 32.714580273],
#     [-0.811308686707, -0.277276065449, 29.5977405865],
#     [-0.812331692291, -0.40975978382, 27.6210856615],
#     [-0.816075037319, -0.535615685086, 27.2420699235],
#     [-0.823691366944, -0.654350489595, 29.1823292975],
#     [-0.836688691603, -0.765630198427, 34.2275056775],
#     [-0.854984518665, -0.86845932028, 43.029581434],
#     [-0.879261949054, -0.961799684483, 55.9594146815],
#     [-0.740499820944, 0.901631050387, 97.0261463995],
#     [-0.735011699497, 0.82881933383, 84.971061395],
#     [-0.733021568161, 0.740454485354, 73.733621269],
#     [-0.732821755233, 0.638770044767, 63.3815970475],
#     [-0.733876941678, 0.525818698874, 54.0655910105],
#     [-0.735055978521, 0.403303715698, 45.90859502],
#     [-0.736448900325, 0.273425879041, 38.935709456],
#     [-0.737556181137, 0.13826504904, 33.096106049],
#     [-0.738278724065, -9.73058423274e-06, 28.359664343],
#     [-0.738507612286, -0.138781586244, 24.627237837],
#     [-0.738539663773, -0.275090412979, 21.857410904],
#     [-0.739099040189, -0.406068448513, 20.1110519655],
#     [-0.741152200369, -0.529726022182, 19.7019157715],
# ])

DATA = np.loadtxt(open("surface.dat", "rb"))

Xs = DATA[:,0]
Ys = DATA[:,1]
Zs = DATA[:,2]
value = np.loadtxt(open("value.dat", "rb"))
# value = DATA[:,3]

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