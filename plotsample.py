# %%
import plotly.offline as pltly
import plotly.figure_factory as ff
import plotly.graph_objs as go
import re
import numpy as np
from scipy.spatial import Voronoi
from pymatgen import Structure, Lattice
from pymatgen.core.periodic_table import Element
# %%

# plot
pltly.iplot(go.Figure(data=data, layout=layout))


# %%
