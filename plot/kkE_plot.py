import numpy as np
import plotly.offline as pltly
import plotly.graph_objs as go
import re

'''
python kkE_plot.py /path/to/calc/dir
This code automatically search the result files:
- structure.hdf5
- fermi.hdf5
- dense.out
This code calculates & plots:
- 3D band structure at a 2D k-plane (default: k_z=0)
- spin texture of the 2D k-plane at a selected energy (default: E_F=0)
'''
# electron or hole の判別
# BZから折り返す