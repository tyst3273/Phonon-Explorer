

# m_ means 'module' (i.e. python file with classes in it), 
# c_ means class. if modules arent in this dir, can import with '.' 
# connecting the path
from file_tools.m_MDE_tools import c_MDE_tools
from file_tools.m_file_utils import c_timer, crash

# --------------------------------------------------------------------------------------------------
# bin raw events file into histo datasets. merge histo datasets into single
# hdf5 file

MDE_timer = c_timer('MDE',units='m')

# raw event file
MDE_file_name = f'../merged_mde/LSNO25_Ei_120meV_300K.nxs'

# load the raw event dataset
MDE_tools = c_MDE_tools(MDE_file_name)

# output hdf5 file
out_file_name = f'LSNO25_300K_test.hdf5'

# binning projections
u = [ 1, 0, 0]
v = [ 0, 1, 0]
w = [ 0, 0, 1]

#H_bins = [  -5.025,   0.050,    15.025]
#K_bins = [ -12.025,   0.050,     7.525]
#L_bins = [ -12.500,   5.000,    12.500]
#E_bins = [ -20.250,   0.500,   100.250]

# binning args: [start, step, stop]
H_bins = [  -0.025,   0.050,    15.025]
K_bins = [  -0.025,   0.050,     7.525]
L_bins = [ -12.500,   5.000,    12.500]
E_bins = [ -20.250,   0.500,   100.250]

# how many chunks along H, K, L to split binning into
num_chunks = [2,2,2]

# bin the file and write hdf5
MDE_tools.bin_MDE_chunks(H_bins,K_bins,L_bins,E_bins,num_chunks,out_file_name,u,v,w)

MDE_timer.stop()

# calls keyboard interrupt exception; 'exit()' doesnt work with ipython
crash()

# --------------------------------------------------------------------------------------------------
# can use these to smooth raw data and then 'calculate' background with rocking scans

BG_timer = c_timer('background')

bg_tools = c_background_tools(out_file_name)

# splits Q-points in file into 'blocks' and Gaussian smooths using fwhm given as arg
# writes new file with suffix 'SMOOTHED'
bg_tools.make_smoothed_file(smoothing_fwhm=1.5,num_blocks=10)

# performs 'rocking scans' around every Q-point in the histo hdf5 file
# and writes background file with suffix 'BACKGROUND'
bg_tools.calculate_background(num_Q_point_procs=16)

# subtracts background from all raw Q-points and writes new file 
# with suffix 'BACKGROUND_SUBTRACTED'
bg_tools.subtract_background(num_blocks=10)

BG_timer.stop()

# --------------------------------------------------------------------------------------------------









