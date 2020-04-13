# %%
import plotly.offline as pltly
import plotly.graph_objs as go
import numpy as np

# %%
with open("./vfermi.frmsf", "r") as f:
    all_data = f.readlines()
    frmsf = [i.rstrip().split() for i in all_data]
nk = frmsf[0]
points = int(nk[0])*int(nk[1])*int(nk[2])
nbands = int(frmsf[2][0])
bvec = np.array(frmsf[3:6]).astype(np.float32)

# %%
X, Y, Z = np.mgrid[0:1:int(nk[0])*1j, 0:1:int(nk[1])*1j, 0:1:int(nk[2])*1j]

# %%
datas = []
for i in range(nbands):
    datas.append(go.Isosurface(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=frmsf[6+points*i:6+points*(i+1)],
        isomin=0,
        isomax=0,
        caps=dict(x_show=False, y_show=False)
    ))
pltly.plot(go.Figure(data=datas))

# %%
'''
k = []
for a in range(24):
    for b in range(24):
        for c in range(24):
            k.append(bvec[0]*a/24+bvec[1]*b/24+bvec[2]*c/24)
# %%
datas = []
for i in range(nbands):
    datas.append(go.Isosurface(
        x=np.array(k).T[0],
        y=np.array(k).T[1],
        z=np.array(k).T[2],
        value=frmsf[6+points*i:6+points*(i+1)],
        isomin=0,
        isomax=0,
        caps=dict(x_show=False, y_show=False)
    ))
pltly.plot(go.Figure(data=datas))

# %%
tmp = np.dot(np.array([X, Y, Z]).T, bvec)
x = tmp[:, :, :, 0]
y = tmp[:, :, :, 1]
z = tmp[:, :, :, 2]
'''
