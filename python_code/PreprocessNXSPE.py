
"""
author: Tyler C. Sterling
Email: ty.sterling@colorado.edu
Affil: University of Colorado Boulder, Raman Spectroscopy and Neutron Scattering Lab

Date: 05/16/2024
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
        file and add instrumentDefinition.directory=$HOME/.mantid/instrument where you should 
        replace $HOME by the actual path. this might break some stuff, however, as mantid 
        needs to find the correct *.xml configuration files for each instrument. i found these
        under /opt/mantid*/instrument where the * means the correct version. i copied them 
        to $HOME/.mantid/instrument and all is well! :)
    - update: this broke some other mantid functionality later and i had to undo it, but 
        this script is still work with the unmodified .xml file! 
"""

import glob
from file_tools.m_NXSPE_tools import bin_NXSPE
from file_tools.m_file_utils import c_timer

"""
### DEV
# reload the modules if running interactively while modifying modules
from importlib import reload
import file_tools.m_NXSPE_tools as m_NXSPE_tools
reload(m_NXSPE_tools)
bin_NXSPE = m_NXSPE_tools.bin_NXSPE
### DEV
"""

# --------------------------------------------------------------------------------------------------

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

H_lo = 1
H_hi = 3
H_step = 0.20
H_bin = 2

K_lo = -6
K_hi = 6
K_step = 0.20
K_bin = 2

L_lo = -4
L_hi = 4
L_step = 0.20
L_bin = 2

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

# this goes and does the stuff
bin_NXSPE(event_files,goniometer_file,hdf5_output_file,UB_params,
          H_lo,H_hi,H_step,K_lo,K_hi,K_step,L_lo,L_hi,L_step,H_bin,K_bin,L_bin)








