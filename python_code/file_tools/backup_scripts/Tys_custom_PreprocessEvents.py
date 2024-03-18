
"""
this file is part of the phonon explorer package!
author: tyler c. sterling & dmitry reznik.
email: ty.sterling@colorado.edu
affil: University of Colorado Boulder, Raman Spectroscopy and Neutron Scattering Lab
date: 03/08/2024
description:
    take binning and projection args from the user and call mantid MDNorm (on the backend)
    to "prebin" data for later analysis with phonon explorer. can be split into "chunks" to 
    save memory. see the description of the variables below.
"""

# m_ means 'module' (i.e. python file with classes in it), 
# c_ means class. if modules arent in this dir, can import with '.' 
# connecting the path

from file_tools.m_MDE_tools import c_MDE_tools
from file_tools.m_background_tools import c_background_tools
from file_tools.m_file_utils import c_timer, crash

### DEV
# reload the modules if running interactively while modifying modules
#from importlib import reload
#import file_tools.m_MDE_tools as m_MDE_tools
#reload(m_MDE_tools)
#c_MDE_tools = m_MDE_tools.c_MDE_tools
### DEV

# --------------------------------------------------------------------------------------------------

# raw MDE file from SNS experiment
MDE_file_name = f'../../../merged_mde/LSNO25_Ei_120meV_300K.nxs'

# binned sparse hdf5 output file
out_file_name = f'LSNO25_300K_test.hdf5'

# binning projections. same meaning as in NormMD
u = [ 1, 0, 0]
v = [ 0, 1, 0]
w = [ 0, 0, 1]

# binning args: [start, step, stop]. these are interpreted like the args to mantid MDNorm and to 
# Horace cut_sqw algorithms. the binning will actually start 'start' with spacing equal to the
# step size. i.e. the  bin centers will be [start+1*step/2, start+2*step/2, start+3*step/2, ...]. 
# e.g. for H_bins = [-0.05, 0.10, 1.05], the bin centers will be H = [ 0, 0.1, ..., 1] integrated 
# -0.05 to 0.05, 0.05 to 0.15, etc.
#H_bins = [  -5.050,   0.100,    15.050]
#K_bins = [ -12.050,   0.100,     7.550]
#L_bins = [ -10.500,   1.000,    10.500]
#E_bins = [ -20.500,   1.000,   100.500]
H_bins = [    3.950,    0.100,     8.050]
K_bins = [   -2.050,    0.100,     2.050]
L_bins = [   -5.000,    2.000,     5.000]
E_bins = [  -20.500,    1.000,   100.500]

# for the binning args above, bin data with this offset. e.g. for 0.0, it is just the binning args 
# above. for 0.1, the bins above are offset by 0.1 etc. it is unecessary to give 0.0 as an offset 
# as this is included automatically and, no matter what offsets are requested by the user, zero
# offset is always calculated (the default). as an example, assume H_offsets = [0.05] and  
# K_offsets = L_offsets = None. Then data will be binned on the grid specified by bin args above 
# and offset along H by 0.05 rlu. e.g. if H_bins = [-0.05, 0.10, 1.05], data in the file will have 
# bin centers H = [ 0, 0.1, ..., 1,  0.05, 0.15, ..., 0.95, 1.05]. keep in mind that the  bins will
# using this scheme.
H_offsets = [0.050]
K_offsets = [0.050]
L_offsets = [1.000]
#H_offsets = None; K_offsets = None; L_offsets = None ## the default

# split binning over these chunks
num_chunks = [2,2,1]

# --------------------------------------------------------------------------------------------------
# you don't have to change anything below here!

# just a timer. dont change this
MDE_timer = c_timer('MDE',units='m')

# load the raw event dataset. dont change this.
MDE_tools = c_MDE_tools(MDE_file_name)

# bin the file and write hdf5. dont change this.
MDE_tools.bin_MDE_with_offsets(H_bins,K_bins,L_bins,E_bins,H_offsets,K_offsets,L_offsets,
        num_chunks,out_file_name,u,v,w)

MDE_timer.stop()

# --------------------------------------------------------------------------------------------------
# can use these to smooth raw data and then 'calculate' background with rocking scans. 
# This is not used in current version of Phonon-Explorer but you can try it.

### THIS STUFF IS JUST FOR DEVELOPERS AND CAN BE IGNORED.

subtract_bg = False

if subtract_bg:

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









