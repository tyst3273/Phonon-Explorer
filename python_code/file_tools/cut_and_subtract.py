

from MDE_tools import c_MDE_tools
from file_utils import c_timer

# --------------------------------------------------------------------------------------------------

MDE_timer = c_timer('MDE',units='m')

T = 300

MDE_file_name = f'../merged_mde/LSNO25_Ei_120meV_300K.nxs'
out_file_name = f'LSNO25_300K_parallel_fine.hdf5'
    
u = [ 1, 0, 0]
v = [ 0, 1, 0]
w = [ 0, 0, 1]

H_bins = [  -5.025,   0.050,    15.025]
K_bins = [ -12.025,   0.050,     7.525]
L_bins = [ -12.500,   5.000,    12.500]
E_bins = [ -20.250,   0.500,   100.250]

num_chunks = [4,4,4]

MDE_tools = c_MDE_tools(MDE_file_name)
MDE_tools.bin_MDE_chunks(H_bins,K_bins,L_bins,E_bins,num_chunks,out_file_name,u,v,w)

MDE_timer.stop()

# --------------------------------------------------------------------------------------------------

BG_timer = c_timer('background')

bg_tools = c_background_tools(out_file_name)
bg_tools.make_smoothed_file(smoothing_fwhm=1.5,num_blocks=10)
bg_tools.calculate_background(num_Q_point_procs=16)
bg_tools.subtract_background(num_blocks=10)

BG_timer.stop()

# --------------------------------------------------------------------------------------------------









