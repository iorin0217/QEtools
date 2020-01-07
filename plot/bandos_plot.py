import numpy as np
import plotly.offline as pltly
import plotly.graph_objs as go
import re

'''
python bandos_plot.py /path/to/calc/dir
This code automatically search the result files:
- structure.hdf5
- fermi.hdf5
- 
- prefix.projwfc.pdos_atm#num(atom)_wfc#num(orb)
- 
This code plots:
- band
- BZ of it
'''

def bandos_plot():
    # use subplot
    band_set
    dos_set


def band_set(band_gnu_file, fermi_energy):
    '''eV単位にしてFermi energyを0点にとる
    '''
    band_data = np.loadtxt(band_gnu_file)
    k_list = np.unique(band_data.T[0])  # k axis
    num_band = int(len(band_data)/len(k_list))
    # energies at k axis
    bands = np.array(np.split(band_data.T[1], num_band)) - float(fermi_energy)
    return k_list, bands


def get_k_points(band_gp_file):
    k_points = {}  # k point name and corrd
    with open(band_gp_file, "r") as f:
        band_gp = f.readlines()
    k_points_corrd = [float(s[:-1].split("=")[1])
                      for s in band_gp if re.match("x\d+ = ", s)]
    k_points_name = [s.split(" ")[0][1:-1]
                     for s in band_gp if re.match("\S+\sx\d+", s)]
    k_points = {corrd: K for (corrd, K) in zip(k_points_corrd, k_points_name)}
    return k_points


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


def dos_set(dos_gnu_file):
    return None


def dos_plot():
    with open('213p5_U.projwfc.pdos_atm#5(Pd)_wfc#3(d)','r') as dat:
    lines = dat.readlines()
    begin = 1
    end = len(lines)
    names=['total',r'$ z^2 $','xz','yz',r'$ x^2-y^2 $','xy']
    trace=[]
    for j in range(1,7):
        energy = []
        DOS = []
        for i in range(begin, end):
            e = float(lines[i].split()[0]) - 8.6327
            dos = float(lines[i].split()[j])
            energy.append(e)
            DOS.append(dos)

        trace1 = go.Scatter(
            x = energy,
            y = DOS,
            name=names[j-1]
        )
        trace.append(trace1)
    layout = go.Layout(
        xaxis=dict(title = 'energy[eV]'),
        yaxis=dict(title = 'DOS[E]'),
    )

    fig = go.Figure(data=trace, layout=layout)

    ofl.plot(fig,filename='p-DOS.html', auto_open=False, include_mathjax='cdn')

    return None
