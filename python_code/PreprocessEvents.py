
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
MDE_file_name = f'../../../merged_mde/LSNO25_Ei_120meV_300K.nxs'

# binned sparse hdf5 output file
merged_file_name = f'LSNO25_300K_test.hdf5'

# binning projections. same meaning as in NormMD
u = [ 1, 0, 0]
v = [ 0, 1, 0]
w = [ 0, 0, 1]

# H_lo = lower bin center. H_hi = *requested* upper bin center. if it's not commensurate with 
# H_bin, it will be tweaked! H_bin = widths of the bins centered on the requested array of 
# bin-centers. e.g. for H_lo = 0.0, H_hi = 1.0, H_bin = 0.1, the bin centers will be H = [0.0, 
# 0.1, 0.2,..., 1.0] which are integrated from -0.05 to 0.05, 0.05 to 0.15, etc. e.g. for 
# H_lo = 0.0, H_hi = 1.02, H_bin = 0.1, you will get the same binning as above. K_* and L_* 
# have the same meaning. H_step, K_step, and L_step are "offsets" applied to shift the binning 
# grid centers by a fixed amount so that users can have the same binning with arbitrary bin 
# centers. e.g. if H_step = 0.03 and H_lo = 0.0, H_hi = 1.0, and H_bin = 0.1, the bin centers 
# will be H = [0.03, 0.13, ..., 1.03] intergrated from -0.02 to 0.08, 0.08 to 0.18 etc.

H_lo = 5
H_hi = 7
H_bin = 0.10
H_step = 0.0

K_lo = -2
K_hi = 2
K_bin = 0.10
K_step = 0.05

L_lo = -5
L_hi = 5
L_bin = 5.0
L_step = 0.0

E_lo = -20
E_hi = 100
E_bin = 1.0

# split binning over these chunks
num_chunks = [1,1,1]

# whether or not to append to file (or overwrite).
# e.g. for when you want to add multiple binning schemes to a single file 
append = True

# --------------------------------------------------------------------------------------------------
# you don't have to change anything below here!

# load the raw event dataset. dont change this.
bin_MDE(MDE_file_name,H_lo,H_hi,H_bin,K_lo,K_hi,K_bin,L_lo,L_hi,L_bin,E_lo,E_hi,E_bin,
        H_step,K_step,L_step,merged_file_name,u,v,w,num_chunks,append)


