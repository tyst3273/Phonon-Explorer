
"""
author: Tyler C. Sterling
Email: ty.sterling@colorado.edu
Affil: University of Colorado Boulder, Raman Spectroscopy and Neutron Scattering Lab

Date: 03/05/2024
Description:
    - parse files NXSPE files from JPARC, merge into MDE dataset, use MDNorm to bin, 
        then write the histogrammed data to .hdf5 file that can be used by phonon 
        explorer. 
    - all the code to talk to the JPARC files was written by Andrei. Tyler wrote the 
        wrapper that loops over Q-pt grids and he wrote the code that writes histo 
        data to file.

    - note: to get mantid to correctly accept files from JPARC, you have to use a 
        modified "Facilities.xml" file. there is one provied in 'file_tools' dir. 
        allegedly we can just copy this to our $HOME/.mantid/instrument dir, but that 
        didnt work recently (it used to). I had to edit $HOME/.mantid/Mantid.user.properties 
        file and add instrumentDefinition.directory=/SNS/users/YOUR-USER-NAME/.mantid/instrument
        where you should replace YOUR-USER-NAME by your user name on the SNS server. 
        for other computers, the path should ultimately it should point to 
        $HOME/.mantid/instrument
"""

import glob
from file_tools.m_PreprocessNXSPE import bin_NXSPE_with_offsets
from file_tools.m_file_utils import c_timer

### DEV
# reload the modules if running interactively while modifying modules
#from importlib import reload
#import file_tools.m_PreprocessNXSPE as m_PreprocessNXSPE
#reload(m_PreprocessNXSPE)
#bin_NXSPE_with_offsets = m_PreprocessNXSPE.bin_NXSPE_with_offsets
### DEV

# --------------------------------------------------------------------------------------------------

# binning args: [start, step, stop]. these are interpreted like the args to mantid MDNorm and to 
# Horace cut_sqw algorithms. the binning will actually start 'start' with spacing equal to the
# step size. i.e. the  bin centers will be [start+1*step/2, start+2*step/2, start+3*step/2, ...]. 
# e.g. for H_bins = [-0.05, 0.10, 1.05], the bin centers will be H = [ 0, 0.1, ..., 1] integrated 
# -0.05 to 0.05, 0.05 to 0.15, etc.
H_bins = [  0.90,  0.20,  3.10]
K_bins = [ -6.10,  0.20,  6.10]
L_bins = [ -4.10,  0.20,  4.10]

# for the binning args above, bin data with this offset. e.g. for 0.0, it is just the binning args 
# above. for 0.1, the bins above are offset by 0.1 etc. it is unecessary to give 0.0 as an offset 
# as this is included automatically and, no matter what offsets are requested by the user, zero
# offset is always calculated (the default). as an example, assume H_offsets = [0.05] and  
# K_offsets = L_offsets = None. Then data will be binned on the grid specified by bin args above 
# and offset along H by 0.05 rlu. e.g. if H_bins = [-0.05, 0.10, 1.05], data in the file will have 
# bin centers H = [ 0, 0.1, ..., 1,  0.05, 0.15, ..., 0.95, 1.05]. keep in mind that the  bins will
# using this scheme.
H_offsets = [  0.10]
K_offsets = [  0.00] 
L_offsets = [  0.00]

# event files (i.e. neutron => detector) already binned in energy
event_files = sorted(
    glob.glob('/SNS/ARCS/IPTS-26347/shared/jpark/nxspe120Ei_010K_Ebin1_RCmasked/S0*.nxspe'))

# file that lists the goniometer angles for each NXPSE file
goniometer_file = '/SNS/ARCS/IPTS-26347/shared/jpark/120meV300Hz010K_run_list.txt'
hdf5_output_file = 'nxspe120Ei_010K_Ebin1ZBk.hdf5'

# this should be aligned UV matrix
UB_params={'a':3.93,'b':3.93,'c':3.93,'alpha':90,'beta':90,'gamma':90,
    'u':'1,0,-0.1','v':'0,1,-0.1'}

# --------------------------------------------------------------------------------------------------
# you don't have to change anything below here!

_timer = c_timer('bin_with_offsets',units='m')

# this goes and does the stuff
bin_NXSPE_with_offsets(event_files,goniometer_file,hdf5_output_file,UB_params,
    H_bins,K_bins,L_bins,H_offsets,K_offsets,L_offsets)

_timer.stop()







