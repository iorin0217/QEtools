
# %%
import xml.etree.ElementTree as ET
import numpy as np
import plotly.offline as pltly
import plotly.graph_objs as go
import re
# %%
fermi_energy = 11.1498
orb_p = [r"$z$", r"$x$", r"$y$"]
orb_d = [r"$z^2$", r"$xz$", r"$yz$", r"$x^2-y^2$", r"$xy$"]
orb_f = [r"$z^3$", r"$xz^2$", r"$yz^2$", r"$zx^2-zy^2$",
         r"$xyz$", r"$x^3-3xy^2$", r"$3yx^2-y^3$"]
# TODO : jlmと軌道の対応
with open("/home/CMD35/cmd35stud07/experiments/SRO/pwscf.pdos_atm#3(Ru)_wfc#5(d_j2.5)", "r") as f:
    tmp = f.readlines()
    atm_3_wfc_5 = np.array(
        [np.array(list(map(float, i.rstrip().split()))) for i in tmp[1:]])
atm_3_wfc_5[0]
# %%
dosxaxis = go.XAxis(
    title="Density of states",
    showgrid=True,
    showline=True,
    range=[.01, 3],
    mirror="ticks",
    ticks="inside",
    linewidth=2,
    tickwidth=2
)
dosyaxis = go.YAxis(
    title=r"$E - E_f \quad / \quad \\text{eV}$",
    showgrid=True,
    showline=True,
    ticks="inside",
    mirror='ticks',
    linewidth=2,
    tickwidth=2,
    zerolinewidth=2
)
dos_layout = go.Layout(
    title="Density of states of Silicon",
    xaxis=dosxaxis,
    yaxis=dosyaxis
)
data = go.Scatter(x=atm_3_wfc_5.T[1], y=atm_3_wfc_5.T[0]-fermi_energy,
                  mode="lines", name="Rud", line=go.Line(color="#444444"), fill="tozeroy")
bandfig = go.Figure(data=go.Data([data]), layout=dos_layout)
pltly.iplot(bandfig)
# include_mathjax='cdn')

# %%
# pdos.out / fat.out : orbital
# band.out.gnu : [linear k , E] * nband
# atomic_proj.xml : [k, orb, a+bj * band]
RytoeV = 13.605698066
projxml = ET.parse(
    "/home/CMD35/cmd35stud07/experiments/SRO/pwscf.save/atomic_proj.xml").getroot()
nkpoints = int(projxml.findall(
    "HEADER/NUMBER_OF_K-POINTS")[0].text.strip())
# Read the number of BANDS
nbands = int(projxml.find("HEADER/NUMBER_OF_BANDS").text)
# get number of projections
nproj = int(projxml.find("HEADER/NUMBER_OF_ATOMIC_WFC").text)
# get weights of kpoints projections
weights = list(
    map(float, projxml.find("WEIGHT_OF_K-POINTS").text.split()))
# get kpoints
kpoints_lines = projxml.find("K-POINTS").text.strip().split('\n')
kpoints_float = [list(map(float, kline.split())) for kline in kpoints_lines]
kpoints = np.array(kpoints_float)
# %%
proj = np.zeros([nkpoints, nproj, nbands], dtype=complex)
for ik in range(nkpoints):
    for ip in range(nproj):
        projlist = projxml.find("PROJECTIONS/K-POINT.%d/ATMWFC.%d" %
                                (ik+1, ip+1)).text.splitlines()[1:-1]
        proj[ik, ip] = [(lambda x, y: complex(float(x), float(y)))(
            *c.split(',')) for c in projlist]
np.array(proj)
# %%
eigen = []
for ik in range(nkpoints):
    eigen.append(list(map(float, projxml.find(
        "EIGENVALUES/K-POINT.%d/EIG" % (ik+1)).text.split())))
np.array(eigen)*RytoeV
# %%
kpoints = []


def band_plot(emin, emax, band_gnu_file, fermi_energy, band_gp_file):
    k_list, bands = band_set(band_gnu_file, fermi_energy)
    band_traces = []  # spaghetti
    for band in bands:
        band_traces.append(
            go.Scatter(
                x=k_list,
                y=band,
                mode="lines",
                line=go.Line(color="#ff2b2b"),
                showlegend=False
            )
        )
    k_points = get_k_points(band_gp_file)
    vlines = []  # auxiliary lines
    for corrd in k_points:
        vlines.append(
            go.Scatter(
                x=[corrd, corrd],
                y=[np.min(bands), np.max(bands)],
                mode="lines",
                line=go.Line(color="#aaaaaa", width=1),
                showlegend=False
            )
        )
    annotations = []
    for corrd, K in k_points.items():
        if K == "GAMMA":
            K = r"\Gamma"
        K = r"$"+K+"$"
        annotations.append(
            go.Annotation(
                x=corrd, y=-0.01,
                xref="x1", yref="paper",
                text=K,
                xanchor="center", yanchor="top",
                showarrow=False
            )
        )
    bandxaxis = go.XAxis(
        title="k-points",
        range=[0, k_list[-1]],
        showgrid=False,
        showline=True,
        ticks="",
        showticklabels=False,
        mirror=True,
        linewidth=2
    )
    bandyaxis = go.YAxis(
        title="$E - E_f \quad / \quad \\text{eV}$",
        range=[emin, emax],
        showgrid=True,
        showline=True,
        zeroline=True,
        mirror="ticks",
        ticks="inside",
        linewidth=2,
        tickwidth=2,
        zerolinewidth=2
    )
    bandlayout = go.Layout(
        title="Bands diagram of SrP8",
        xaxis=bandxaxis,
        yaxis=bandyaxis,
        annotations=go.Annotations(annotations)
    )
    bandfig = go.Figure(data=band_traces + vlines, layout=bandlayout)
    pltly.plot(bandfig, filename="Bands_SrP8", include_mathjax='cdn')
