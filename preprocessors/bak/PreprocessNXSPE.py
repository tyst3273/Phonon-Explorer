
"""
author: Tyler C. Sterling
Email: ty.sterling@colorado.edu
Affil: University of Colorado Boulder, Raman Spectroscopy and Neutron Scattering Lab

Date: 08/22/2024
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

# ...

H_lo = 0
H_hi = 6
H_step = 0.05
H_bin = 1

K_lo = -3
K_hi = 3
K_step = 0.05
K_bin = 1

L_lo = -3
L_hi = 3
L_step = 0.05
L_bin = 1

# event files (i.e. neutron => detector) already binned in energy
event_files = sorted(
    glob.glob('/SNS/ARCS/IPTS-26347/shared/jpark/nxspe054Ei_335K_Ebin0p5_RCmasked/S0*.nxspe'))

# file that lists the goniometer angles for each NXPSE file
goniometer_file = '/SNS/ARCS/IPTS-26347/shared/jpark/120meV300Hz335K_run_list.txt'
hdf5_output_file = 'nxspe054Ei_335K_Ebin0p5Hres2_symm.hdf5'

# this should be aligned UV matrix
UB_params={'a':3.92,'b':3.92,'c':3.92,'alpha':90,'beta':90,'gamma':90,
    'u':'1,0,-0.0935','v':'0,1,-0.0987'}

# symmetry args for MDNorm
# -- dunno what you actually need, this is an example with identity, reflections, and inversion
# -- you can just set to None if you dont want to use symmetry
SymmetryOperations = 'x,y,z;-x,y,z;x,-y,z;x,y,-z;-x,-y,-z'
#SymmetryOperations = None

# --------------------------------------------------------------------------------------------------
# you don't have to change anything below here!

# this goes and does the stuff
bin_NXSPE(event_files,goniometer_file,hdf5_output_file,UB_params,
          H_lo,H_hi,H_step,K_lo,K_hi,K_step,L_lo,L_hi,L_step,H_bin,K_bin,L_bin,
          SymmetryOperations=SymmetryOperations)








