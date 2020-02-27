import numpy as np
from scipy.spatial import Voronoi
import plotly.offline as pltly
import plotly.graph_objs as go
import re

'''
python structure_plot.py /path/to/env/file
This code automatically search the input files:
- env.hdf5
This code plots:
- crystal structure at the input stage
- crystal structure at the submission stage
- BZ of it
* if you want to just view the input structure, you can see it in VESTA
'''


def draw_structure(structure_title, crystal_structure, spin_structure):
    structure_layout = go.Layout(title=structure_title, xaxis=,)
    structure_figure = go.Figure(data=structure, layout=structure_layout)
    return structure_figure


def draw_BZ(BZ_title, bvec, klabel):
    # special_point_annotation
    kpoints = []
    for i in range(-1, 2):
        for j in range(-1, 2):
            for k in range(-1, 2):
                vec = i*bvec[0] + j*bvec[1] + k*bvec[2]
                kpoints.append(vec)
    wigner_seitz = Voronoi(np.array(kpoints))
    faces = []
    for idict in wigner_seitz.ridge_dict:
        if idict[0] == 13 or idict[1] == 13:
            faces.append(wigner_seitz.ridge_dict[idict])
    verts = wigner_seitz.vertices
    poly = []
    for ix in range(len(faces)):
        temp = []
        for iy in range(len(faces[ix])):
            temp.append(verts[faces[ix][iy]])
        poly.append(np.array(temp))
    BZ = []
    for iface in np.array(poly):
        iface = np.pad(iface, ((0, 1), (0, 0)), 'wrap')
        x, y, z = iface[:, 0], iface[:, 1], iface[:, 2]
        plane = go.Scatter3d(x=x, y=y, z=z, mode='lines',
                             line=dict(color='black', width=4))
        BZ.append(plane)
    annotations = []
    i = 1
    for axis in bvec:
        BZ.append(go.Scatter3d(
            x=[0, axis[0]],
            y=[0, axis[1]],
            z=[0, axis[2]],
            mode='lines',
            marker=dict(
                color='rgb(100,100,200)',
                size=5,
                opacity=0.8
            )
        ))
        annotations.append(
            dict(
                showarrow=False, x=axis[0], y=axis[1], z=axis[2], text=f"b{i}"
            ))
        i += 1
    BZ_layout = go.Layout(title=BZ_title, scene=dict(
        xaxis_showgrid=False,
        xaxis_showline=False,
        xaxis_showticklabels=False,
        xaxis_zeroline=False,
        xaxis_autorange=True,
        xaxis_showbackground=False,
        xaxis_title=' ',
        yaxis_showgrid=False,
        yaxis_showline=False,
        yaxis_showticklabels=False,
        yaxis_zeroline=False,
        yaxis_autorange=True,
        yaxis_title=' ',
        yaxis_showbackground=False,
        zaxis_showgrid=False,
        zaxis_showline=False,
        zaxis_showticklabels=False,
        zaxis_zeroline=False,
        zaxis_autorange=True,
        zaxis_showbackground=False,
        zaxis_title=' ',
        annotations=annotations),
        showlegend=False)
    BZ_figure = go.Figure(data=BZ, layout=BZ_layout)
    return BZ_figure


if __name__ == '__main__':
    import sys
    import os
    from pymatgen import Structure
    env_file = sys.argv[1]
    dir_name = os.path.dirname(env_file)
    env =
    structure_input_figure = draw_structure(
        dir_name+"_input", env["raw_crystal"], env["raw_spin"])
    pltly.plot(structure_input_figure, filename=dir_name+"_input")
    structure_output_figure = draw_structure(
        dir_name+"_output", env["crystal"], env["spin"])
    pltly.plot(structure_output_figure, filename=dir_name+"_output")
    BZ_figure = draw_BZ(dir_name+"_BZ", env["bvec"], env["klabel"])
    pltly.plot(BZ_figure,
               filename=dir_name+"_BZ", include_mathjax='cdn')
