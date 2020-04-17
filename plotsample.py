
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
kpoints_linear = np.array([0., 0.02586777, 0.05173553, 0.0776033, 0.10347106,
                           0.12933883, 0.15520659, 0.18107436, 0.20694212, 0.23280989,
                           0.25867765, 0.28454542, 0.31041318, 0.33628095, 0.36214871,
                           0.38801648, 0.41388424, 0.43975201, 0.46561977, 0.49148754,
                           0.5173553, 0.54322307, 0.56909083, 0.5949586, 0.62082636,
                           0.64669413, 0.67256189, 0.69842966, 0.72429742, 0.75016519,
                           0.77603295, 0.80190072, 0.82776848, 0.85363625, 0.87950401,
                           0.90537178, 0.93123954, 0.95710731, 0.98297507, 1.00884284,
                           1.0347106, 1.06057837, 1.08644613, 1.1123139, 1.13818166,
                           1.16861885, 1.19905605, 1.22949324, 1.25993043, 1.29036762,
                           1.32080482, 1.35124201, 1.3816792, 1.40764101, 1.43360281,
                           1.45956462, 1.48552642, 1.51148823, 1.53745003, 1.56341184,
                           1.58937365, 1.61533545, 1.64129726, 1.66725906, 1.69322087,
                           1.71918267, 1.74514448, 1.77110628, 1.79706809, 1.82302989,
                           1.8489917, 1.87495351, 1.90091531, 1.92687712, 1.95283892,
                           1.97880073, 2.00476253, 2.03072434, 2.05668614, 2.08264795,
                           2.10860976, 2.13457156, 2.16053337, 2.18649517, 2.21277157,
                           2.23904796, 2.26532436, 2.29160076, 2.31787715, 2.34415355,
                           2.37042995, 2.39670634, 2.42298274, 2.44925914, 2.47553553,
                           2.50181193, 2.52808833, 2.55436472, 2.58064112, 2.60691752,
                           2.63319392, 2.65947031, 2.68574671, 2.71202311, 2.7382995,
                           2.7645759, 2.7908523, 2.81712869, 2.84340509, 2.86968149,
                           2.89595788, 2.92223428, 2.94851068, 2.97478707, 3.00106347,
                           3.02733987, 3.05439515, 3.08145043, 3.10850571, 3.13556099,
                           3.16261628, 3.18967156, 3.21672684, 3.24378212, 3.27083741,
                           3.29789269, 3.32494797, 3.35200325, 3.37905854, 3.40611382,
                           3.4331691, 3.46022438, 3.48727966, 3.51433495, 3.54044729,
                           3.56655964, 3.59267198, 3.61878433, 3.64489667, 3.67100901,
                           3.69712136, 3.7232337, 3.74934605, 3.77545839, 3.80157074,
                           3.82768308, 3.85379543, 3.87990777, 3.90602012, 3.93213246,
                           3.9582448, 3.98435715, 4.01046949, 4.03658184, 4.06269418,
                           4.08880653, 4.11491887, 4.14103122, 4.16714356, 4.19325591,
                           4.21936825, 4.24548059, 4.24548059, 4.27131843, 4.29715626,
                           4.32299409, 4.34883192, 4.37466976, 4.40050759, 4.42634542,
                           4.45218325, 4.47802108, 4.50385892, 4.52969675, 4.55553458,
                           4.58137241, 4.60721024, 4.63304808, 4.65888591, 4.68472374,
                           4.71056157, 4.73639941, 4.76223724, 4.78807507, 4.8139129,
                           4.83975073, 4.86558857, 4.8914264, 4.91726423, 4.94310206,
                           4.96893989, 4.99477773, 5.02061556, 5.04645339, 5.07229122,
                           5.09812906, 5.12396689, 5.12396689, 5.15869541, 5.19342393,
                           5.22815246, 5.22815246, 5.25400236, 5.27985226, 5.30570216,
                           5.33155207, 5.35740197, 5.38325187, 5.40910177, 5.43495168,
                           5.46080158, 5.48665148, 5.51250138, 5.53835128, 5.56420119,
                           5.59005109, 5.61590099, 5.64175089, 5.6676008, 5.6934507,
                           5.7193006, 5.7451505, 5.77100041, 5.79685031, 5.82270021,
                           5.84855011, 5.87440001, 5.90024992, 5.92609982, 5.95194972,
                           5.97779962, 6.00364953, 6.02949943, 6.05534933, 6.08119923,
                           6.10704913, 6.13289904, 6.15874894, 6.18459884, 6.21044874,
                           6.23629865, 6.26214855])

kticks = {0.0: 'GAMMA',
          1.13818166031401: 'X',
          1.3816792013920716: 'P',
          2.186495171622272: 'N',
          3.027339865039933: 'GAMMA',
          3.514334947196056: 'M',
          4.245480594390664: 'S|S_0',
          5.123966887656456: 'GAMMA|X',
          5.228152457637798: 'R|G',
          6.262148547970466: 'M'}
# %%
band_traces = []
for band in (np.array(eigen)*RytoeV).T:
    band_traces.append(go.Scatter(
        x=kpoints_linear, y=band-fermi_energy, mode="lines", line=go.Line(color="#666666"), showlegend=False))
# %%
annotations = []
bandlayout = go.Layout(
    width=600,
    height=450,  # energy*75
    autosize=True,
    paper_bgcolor='rgba(0,0,0,0)',
    plot_bgcolor='rgba(0,0,0,0)',
    font={"family": "arial", "size": 12},
    showlegend=False,
    xaxis={"range": [0, 3.514334947196056], "showgrid": False, "showline": True, "linecolor": "#666666", "zeroline": False,
           "linewidth": 2, "mirror": "ticks", "ticks": "inside",  "tickwidth": 2, "tickvals": list(kticks.keys()), "ticktext": list(kticks.values())},
    yaxis={"title": r"$E - E_f \quad / \quad \\text{eV}$",
           "range": [-2.5, 2.5], "showgrid": False, "showline": True, "linecolor": "#666666", "zeroline": False, "linewidth": 2, "mirror": "ticks", "ticks": "outside",  "tickwidth": 2}
)
bandfig = go.Figure(data=band_traces, layout=bandlayout)
pltly.iplot(bandfig)
# %%


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
