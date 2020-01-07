import numpy as np
import plotly.offline as pltly
import plotly.graph_objs as go
import re

'''
python structure_plot.py /path/to/calc/dir
This code automatically search the input files:
- structure.hdf5
This code plots:
- crystal structure at the submission stage
- BZ of it
'''