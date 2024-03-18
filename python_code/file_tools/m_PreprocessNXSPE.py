
"""
PreprocessNXSPE routine:
author: Andrei T. Savici
email: saviciat@ornl.gov

other stuff:
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

import numpy as np
from mantid.simpleapi import *

from file_tools.m_save_MDE_to_hdf5 import save_MDE_to_hdf5
from file_tools.m_file_utils import crash, c_timer

"""
### DEV
# reload the modules if running interactively while modifying modules
from importlib import reload
import file_tools.m_save_MDE_to_hdf5 as m_save_MDE_to_hdf5
reload(m_save_MDE_to_hdf5)
save_MDE_to_hdf5 = m_save_MDE_to_hdf5.save_MDE_to_hdf5
### DEV
"""

# --------------------------------------------------------------------------------------------------

def PreprocessNXSPE(H_bins,K_bins,L_bins,event_files,goniometer_file,UB_params):
    """
    wrap the code to load and preprocess NXSPE files from JPARC measurement.

    the code below was written by Andrei T. Savici at ORNL (saviciat@ornl.gov). I have no idea
    what it does and don't know how to modify it. I am simply wrapping it in a function to 
    programmatically pass binning args in and get the data out.
    """

    msg = '\n*** WARNING ***\n'
    msg += 'PreprocessNXSPE is currently configured for the manganite experiment\n'
    msg += 'at JPARC! It is not generic and must be redone for other experiments!\n'
    print(msg)

    d = np.loadtxt(goniometer_file,skiprows=1,dtype=str)
    runs = d[:,0]
    gon = d[:,11].astype(float)-82.388
    gon_dict = dict()
    for r,g in zip(runs, gon):
        gon_dict[r] = g

    event_files = ','.join(event_files)
    wg = Load(event_files)
    for wsn in wg.getNames():
        ang = gon_dict[wsn.split('S0')[1]]
        SetGoniometer(Workspace=wsn, Axis0 = f'{ang},0,1,0,1')

    axis_deltaE = wg[0].readX(0)
    emin = axis_deltaE[0]
    emax = axis_deltaE[-1]
    de_step = axis_deltaE[1]-emin

    ConvertFromDistribution(wg)
    wge = ConvertToEventWorkspace(wg)
    wge = CropWorkspaceForMDNorm(wge,XMin=emin,XMax=emax)

    SetUB(wge,**UB_params)

    AddSampleLog(Workspace='wge', LogName='gd_prtn_chrg', LogText='1.0', LogType='Number',
        NumberType='Double')

    mdparts = ConvertToMD(InputWorkspace=wge, QDimensions='Q3D', dEAnalysisMode='Direct',
        Q3DFrames='Q_sample')

    mde = MergeMD(mdparts)
    #SaveMD(InputWorkspace="mde",Filename="/SNS/ARCS/IPTS-26347/shared/jpark/MDEventFile.nxs")

    #UB_params = {'a':3.93,'b':3.93,'c':3.93,'alpha':90,'beta':90,'gamma':90,
    #    'u':'1,0,-0.1','v':'0,1,-0.1'}
    #SetUB(mde,**UB_params)

    #print(axis_deltaE)
    # test MDNorm - elastic slice
    MDNorm(InputWorkspace=mde,
           QDimension0='1,0,0',
           QDimension1='0,1,0',
           QDimension2='0,0,1',
           Dimension0Name='QDimension0',
           Dimension0Binning="{},{},{}".format(H_bins[0],H_bins[1],H_bins[2]),
           Dimension1Name='QDimension1',
           Dimension1Binning="{},{},{}".format(K_bins[0],K_bins[1],K_bins[2]),
           Dimension2Name='QDimension2',
           Dimension2Binning="{},{},{}".format(L_bins[0],L_bins[1],L_bins[2]),
           Dimension3Name='DeltaE',
           Dimension3Binning="{}".format(de_step),
    #       SymmetryOperations='x,y,z;-x,-y,z;x,y,-z;-x,-y,-z;y,x,z;y,x,-z',
           OutputWorkspace='o',
           OutputDataWorkspace='d',
           OutputNormalizationWorkspace='n')

    #SaveMD(InputWorkspace="o",Filename="/SNS/ARCS/IPTS-26347/shared/jpark/MDHistoFile.nxs")

    MD_workspace = mtd['o']

    # to save data from MDNorm, we just need the output MD_workspace
    return MD_workspace

    # ----------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------

def _get_bins(lo,hi,bin,offset=0.0):
    """
    reformat binning args for MDNorm
    """
    lo = float(lo); hi = float(hi); bin = float(bin); offset = float(offset)
    _mod = np.modf((hi-lo)/bin)
    num_bins = _mod[1].astype(int)
    decimal = _mod[0].round(6)
    if decimal != 0.0:
        msg = '\n*** WARNING ***\n'
        msg += f'upper bin center {hi: 6.3f} is not commensurate with\n'
        msg += f'lower bin center {lo: 6.3f} and bin size {bin:6.3f}\n'
        hi = lo+bin*(num_bins+1)
        msg += f'tweaking upper bin center to {hi: 6.3f}\n'
        print(msg)
    lo = np.round(lo-bin/2+offset,6)
    hi = np.round(hi+bin/2+offset,6)
    bin = np.round(bin,6)
    return [lo,bin,hi]

# --------------------------------------------------------------------------------------------------

def bin_NXSPE(event_files,goniometer_file,output_file,UB_params,H_lo,H_hi,H_bin,K_lo,K_hi,K_bin,
            L_lo,L_hi,L_bin,H_step=0.0,K_step=0.0,L_step=0.0,append=True):
    """
    wrapper to translate variable names for interface with PreprocessNXSPE. note, E-binning is 
    fixed by the binning already present in the *.nxspe files

    H_lo = lower bin center. 
    H_hi = *requested* upper bin center. if it's not commensurate with H_bin, it will be tweaked!
    H_bin = widths of the bins centered on the requested array of bin-centers.
    
        e.g. for H_lo = 0.0, H_hi = 1.0, H_bin = 0.1, the bin centers will be H = [0.0, 0.1, 0.2,
        ..., 1.0] which are integrated from -0.05 to 0.05, 0.05 to 0.15, etc.

        e.g. for H_lo = 0.0, H_hi = 1.02, H_bin = 0.1, you will get the same binning as above.

    K_* and L_* have the same meaning. 

    H_step, K_step, and L_step are "offsets" applied to shift the binning grid centers by a fixed
    amount so that users can have the same binning with arbitrary bin centers.
        
        e.g. if H_step = 0.03 and H_lo = 0.0, H_hi = 1.0, and H_bin = 0.1, the bin centers will
        be H = [0.03, 0.13, ..., 1.03] intergrated from -0.02 to 0.08, 0.08 to 0.18 etc.

    output_file is the name where the output histogram data will be written.

    append deterines whether or not the file is overwritten

    see the descriptions in the $ROOT/python_code/PreprocessNXSPE.py script for explanation of
    the other args to this function
    """

    timer = c_timer('bin_NXSPE',units='m')

    # get binning args for c_MDE_tools.bin_MDE_chunks
    H_bins = _get_bins(H_lo,H_hi,H_bin,H_step)
    K_bins = _get_bins(K_lo,K_hi,K_bin,K_step)
    L_bins = _get_bins(L_lo,L_hi,L_bin,L_step)
    
    # bin and save data with on offset grids
    MD_workspace = PreprocessNXSPE(H_bins,K_bins,L_bins,event_files,goniometer_file,UB_params)

    # append to hdf5 file
    save_MDE_to_hdf5(MD_workspace,output_file,append)

    timer.stop()

# --------------------------------------------------------------------------------------------------
