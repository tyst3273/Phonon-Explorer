
from MDE_tools import c_MDE_tools, c_timer
import numpy as np

MDE_file = '../merged_mde/LSNO25_Ei_120meV_300K.nxs'

out_file = 'no_couple_6k0_300K.hdf5'

u = [ 1, 0, 0]
v = [ 0, 1, 0]
w = [ 0, 0, 1]

H_bins = [  5.925, 6.075]
K_bins = [  -6, 0.05, 2]
L_bins = [ -2.5, 2.5]
E_bins = [ -10, 0.25, 100]


MDE_tools = c_MDE_tools(MDE_file)
MDE_tools.bin_MDE(H_bins,K_bins,L_bins,E_bins,u,v,w)
MDE_tools.save_to_hdf5(out_file)



