import numpy as np
from scipy.spatial import Voronoi
import plotly.offline as pltly
import plotly.graph_objs as go
import re
from mendeleev import element
from pymatgen import Structure, Lattice
from pymatgen.core.periodic_table import Element
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


def draw_crystal(crystal_title, structure):
    # TODO : spin structure
    data = []
    annotations = []
    # 球
    dist = np.linspace(0, 1, 12)
    theta = 2*np.pi*dist
    phi = np.arccos(2*dist-1)
    rx = np.outer(np.cos(theta), np.sin(phi))
    ry = np.outer(np.sin(theta), np.sin(phi))
    rz = np.outer(np.ones(12), np.cos(phi))
    ''' TODO : bonds
    bonds = {}
    for i in range(len(structure)):
    nears = BrunnerNN_real().get_nn_info(structure,i)
    bonds[structure[i]] = [near["site"] for near in nears if all([0<= n <= 1 for n in list(near["image"])]) & (not bonds.get(near["site"]))]
    '''
    # make supercell to get boundary atoms
    structure_copy = structure.copy()
    structure_copy.make_supercell(
        [[2, 0, 0], [0, 2, 0], [0, 0, 2]], to_unit_cell=False)
    sites = [site_copy for site_copy in structure_copy if all(
        [-0.00001 <= frac <= 0.50001 for frac in site_copy.frac_coords])]
    # node
    for site in set(sites):
        coord = site.coords
        el = site.species_string
        r = float(Element(el).atomic_radius)
        atomTrace = go.Mesh3d(
            x=(coord[0]+r*rx).flatten(),
            y=(coord[1]+r*ry).flatten(),
            z=(coord[2]+r*rz).flatten(),
            alphahull=0,
            color=element(el).jmol_color,
            hoverinfo='none'
        )
        data.append(atomTrace)
    # edge

    # unitcell

    def get_edgeTrace(edge):
        x, y, z = zip(*edge)
        edgeTrace = go.Scatter3d(
            x=x,
            y=y,
            z=z,
            mode='lines',
            line=dict(color="#666666"),
            showlegend=False,
            hoverinfo='none'
        )
        return edgeTrace

    v0 = [0, 0, 0]
    v1, v2, v3 = structure.lattice.matrix
    v4, v5, v6, v7 = v2+v3, v1+v3, v1+v2, v1+v2+v3
    edges = [[v0, v1], [v0, v2], [v0, v3], [v1, v5], [v1, v6], [v2, v4],
             [v2, v6], [v3, v4], [v3, v5], [v4, v7], [v5, v7], [v6, v7]]
    for edge in edges:
        data.append(get_edgeTrace(edge))
    # axis
    xyz = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    xyz_name = ["x", "y", "z"]
    scale = np.max(np.array([site.coords for site in sites]), axis=0)
    for i, axis in enumerate(xyz):
        data.append(go.Scatter3d(
            x=[0, axis[0]*scale[0]*2],
            y=[0, axis[1]*scale[1]*2],
            z=[0, axis[2]*scale[2]*2],
            mode='lines',
            line=dict(
                color='rgb(100,100,200)',
                width=5
            ),
            hoverinfo='none'
        ))
        annotations.append(
            dict(
                showarrow=False, x=axis[0]*scale[0]*2, y=axis[1]*scale[1]*2, z=axis[2]*scale[2]*2, text=f"{xyz_name[i]}"
            ))
    abc = structure.lattice.matrix / \
        np.linalg.norm(structure.lattice.matrix, axis=1)
    abc_name = ["a", "b", "c"]
    displace = scale[0]
    for i, axis in enumerate(abc):
        data.append(go.Scatter3d(
            x=[-displace, axis[0]-displace],
            y=[0, axis[1]],
            z=[0, axis[2]],
            mode='lines',
            line=dict(
                color='rgb(100,100,200)',
                width=5
            ),
            hoverinfo='none'
        ))
        annotations.append(
            dict(
                showarrow=False, x=axis[0]-displace, y=axis[1], z=axis[2], text=f"{abc_name[i]}"
            ))
    crystal_figure = go.Figure(data=data, layout=layout(annotations))
    return crystal_figure


def draw_BZ(BZ_title, bvec, klabel):
    # TODO : special_point_annotation
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
    BZ_figure = go.Figure(data=BZ, layout=layout(annotations))
    return BZ_figure


def layout(annotations):
    layout = go.Layout(
        scene=dict(
            xaxis=dict(
                showspikes=False,
                showgrid=False,
                zeroline=False,
                ticks='',
                title='',
                showticklabels=False,
                showbackground=False
            ),
            yaxis=dict(
                showspikes=False,
                showgrid=False,
                zeroline=False,
                ticks='',
                title='',
                showticklabels=False,
                showbackground=False
            ),
            zaxis=dict(
                showspikes=False,
                showgrid=False,
                zeroline=False,
                ticks='',
                title='',
                showticklabels=False,
                showbackground=False
            ),
            annotations=annotations
        ),
        showlegend=False
    )
    return layout


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
