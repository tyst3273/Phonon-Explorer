
from MDE_tools import c_MDE_tools, c_timer
import numpy as np

MDE_file = '../merged_mde/LSNO25_Ei_120meV_5K.nxs'
MDE_tools = c_MDE_tools(MDE_file)

E_shift = 1.2

u = [ 1, 0, 0]
v = [ 0, 1, 0]
w = [ 0, 0, 1]


out_file = 'HK_plane_E5_5K.hdf5'
H_bins = [ -4,  0.05, 12]
K_bins = [ -10, 0.05, 6]
L_bins = [ -1.5, 1.5]
E_bins = [ 5+E_shift-1.5, 5+E_shift+1.5]
MDE_tools.bin_MDE(H_bins,K_bins,L_bins,E_bins,u,v,w)
MDE_tools.save_to_hdf5(out_file)


