
"""
Author: Andrei T. Savici
Email: saviciat@ornl.gov

Updated by: Tyler C. Sterling
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

from importlib import reload
import file_tools.m_PreprocessNXSPE as m_PreprocessNXSPE
reload(m_PreprocessNXSPE)
bin_NXSPE_with_offsets = m_PreprocessNXSPE.bin_NXSPE_with_offsets

# --------------------------------------------------------------------------------------------------

# [start, step, stop]. these are interpreted like the args to mantid MDNorm and to 
# Horace cut_sqw: the binning will actually start at 'start' with spacing equal to 
# the step size. i.e. the bin centers will be [start+1*step/2, start+2*step/2,
# start+3*step/2, ...]

H_bins = [  0.90,  0.20,  3.10]
K_bins = [ -6.10,  0.20,  6.10]
L_bins = [ -4.10,  0.20,  4.10]

# [offset_1, offset_2, ..., offset_N]. loop over the offets in *_offsets and bin using 
# the args in *_bin shifted by the offets. e.g. for H_offset = [0.0, 0.01] and H_bins = 
# [-0.05, 0.1, 0.55], we will generate data on a grid intergrated around H=0.0+-0.05, 
# H=0.1+-0.05, ..., H=0.5+-0.05 and similarly H=0.01+-0.05, H=0.11+-0.05, ..., H=0.51+-0.05 
# etc. similarly for K_offets, L_offsets. 

H_offsets = [  0.10]
K_offsets = [  0.00] 
L_offsets = [  0.00]

# event files (i.e. neutron => detector) binned in energy
event_files = sorted(
    glob.glob('/SNS/ARCS/IPTS-26347/shared/jpark/nxspe120Ei_010K_Ebin1_RCmasked/S0*.nxspe'))

# file that lists the goniometer angles for each NXPSE file
goniometer_file = '/SNS/ARCS/IPTS-26347/shared/jpark/120meV300Hz010K_run_list.txt'
hdf5_output_file = 'nxspe120Ei_010K_Ebin1ZBk.hdf5'

# this should be aligned UV matrix
UB_params={'a':3.93,'b':3.93,'c':3.93,'alpha':90,'beta':90,'gamma':90,
    'u':'1,0,-0.1','v':'0,1,-0.1'}

_timer = c_timer('bin_with_offsets',units='m')

# this goes and does the stuff
bin_NXSPE_with_offsets(event_files,goniometer_file,hdf5_output_file,UB_params,
    H_bins,K_bins,L_bins,H_offsets,K_offsets,L_offsets)

_timer.stop()







