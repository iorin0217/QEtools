# %%
import plotly.offline as pltly
import plotly.graph_objs as go
import numpy as np

# %%
with open("/home/CMD35/cmd35stud07/experiments/SRO_fine/vfermi.frmsf", "r") as f:
    all_data = f.readlines()
    frmsf = [i.rstrip().split() for i in all_data]
nk = frmsf[0]
points = nk[0]*nk[1]*nk[2]
nbands = frmsf[2]
bvec = frmsf[3:6]

# %%
X, Y, Z = np.mgrid[0:nk[0]:16j, 0:nk[1]:16j, 0:nk[2]:16j]
datas = []
for i in range(nbands):
    datas.append(go.Isosurface(
        x=(X@np.array(bvec[0])).flatten(),
        y=(Y@np.array(bvec[1])).flatten(),
        z=(Z@np.array(bvec[2])).flatten(),
        value=frmsf[6+points*i:6+points*(i+1)],
        isomin=0,
        isomax=0,
        caps=dict(x_show=False, y_show=False)
    ))
pltly.iplot(go.Figure(data=datas))

# %%
