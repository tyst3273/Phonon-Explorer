
"""
this file is part of the phonon explorer package!
author: tyler c. sterling & dmitry reznik.
email: ty.sterling@colorado.edu
affil: University of Colorado Boulder, Raman Spectroscopy and Neutron Scattering Lab
date: 05/16/2024
description:
    take binning and projection args from the user and call mantid MDNorm (on the backend)
    to "prebin" data for later analysis with phonon explorer. can be split into "chunks" to 
    save memory. see the description of the variables below.
"""

# m_ means 'module' (i.e. python file with classes in it), 
# c_ means class. if modules arent in this dir, can import with '.' 
# connecting the path

from file_tools.m_MDE_tools import bin_MDE

"""
### DEV
# reload the modules if running interactively while modifying modules
from importlib import reload
import file_tools.m_MDE_tools as m_MDE_tools
reload(m_MDE_tools)
bin_MDE = m_MDE_tools.bin_MDE
### DEV
"""

# --------------------------------------------------------------------------------------------------

# raw MDE file from SNS experiment
MDE_file_name = f'../../merged_mde/LSNO25_Ei_120meV_300K.nxs'

# binned sparse hdf5 output file
merged_file_name = f'LSNO25_300K_test.hdf5'

# binning projections. same meaning as in NormMD
u = [ 1, 0, 0]
v = [ 0, 1, 0]
w = [ 0, 0, 1]

# H_lo = lower bin center. 
# H_hi = *requested* upper bin center. if it's not commensurate with H_step, it will be tweaked!
# H_step = widths of the bins centered on the requested array of bin-centers.
#    
#       e.g. for H_lo = 0.0, H_hi = 1.0, H_step = 0.1, the bin centers will be H = [0.0, 0.1, 0.2,
#       ..., 1.0] which are integrated from -0.05 to 0.05, 0.05 to 0.15, etc.
#
#       e.g. for H_lo = 0.0, H_hi = 1.02, H_step = 0.1, you will get the same binning as above.
#       
#       K_*, L_*, and E_* have the same meaning. 
#
# H_bin, K_bin, and L_bin are number of steps *within* a bin to take. data are integrated 
# once for each step by shifting bin centers by that amount. but keeping bin widths the same:
# e.g. if H_step = 0.1 and H_bin = 2, first bin centers will be shifted by 0.0 (i.e. unshifted),
# then shifted by 0.05 rlu and integrated again. shifts along different directions are 'meshed', 
# i.e. H_bin = 2, K_bin = 2, L_bin = 1, H_step = 0.1, K_step = 0.1, L_step = 0.1 will result in 
# 4 integrations with shifts
#       1) [0.00, 0.00, 0.00]
#       2) [0.00, 0.05, 0.00]
#       3) [0.05, 0.00, 0.00]
#       4) [0.05, 0.05, 0.00]
#
# each shifted integration grid will be appended to the same file. e.g. for if H_bin = 2 and 
# H_lo = 0.0, H_hi = 1.0, and H_step = 0.1, the bin centers will be H = [0.00, 0.05, 0.10, 
# 0.15, ..., 0.95, 1.00, 1.05] intergrated from -0.05 to 0.05, 0.00 to 0.10, 0.05 to 0.15, etc.

H_lo = 5
H_hi = 7
H_step = 0.10
H_bin = 2

K_lo = -2
K_hi = 2
K_step = 0.10
K_bin = 2

L_lo = -4
L_hi = 4
L_step = 2
L_bin = 4

E_lo = -50
E_hi = 100
E_step = 1

# split binning over these chunks
num_chunks = [1,1,1]

# --------------------------------------------------------------------------------------------------
# you don't have to change anything below here!

# load the raw event dataset. dont change this.
bin_MDE(MDE_file_name,H_lo,H_hi,H_step,K_lo,K_hi,K_step,L_lo,L_hi,L_step,E_lo,E_hi,E_step,
        H_bin,K_bin,L_bin,merged_file_name,u,v,w,num_chunks)


