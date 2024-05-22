
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
MDE_file_name = f'/SNS/ARCS/IPTS-26347/shared/tys_stuff/merged_mde/LSNO25_Ei_120meV_300K.nxs'

# binned sparse hdf5 output file
merged_file_name = f'LSNO25_300K_test.hdf5'

# binning projections. same meaning as in NormMD
u = [ 1, 0, 0]
v = [ 0, 1, 0]
w = [ 0, 0, 1]

# ...

H_lo = 4
H_hi = 6
H_step = 0.05
H_bin = 2

K_lo = -2
K_hi = 2
K_step = 0.05
K_bin = 2

L_lo = -4
L_hi = 4
L_step = 1
L_bin = 2

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


