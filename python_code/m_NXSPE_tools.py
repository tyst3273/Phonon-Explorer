
"""
PreprocessNXSPE function:
author: Andrei T. Savici
email: saviciat@ornl.gov

Everything else:
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

    """
    # !!! DEV !!!
    # set _stride != 0 to use only every _stride^th files to test run
    _stride = 0
    if _stride:
        event_files = event_files[::_stride]  
    # !!! DEV !!!
    """

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

# --------------------------------------------------------------------------------------------------

def _setup_shifts(H_step,K_step,L_step,H_bin,K_bin,L_bin):
    """
    setup shifts to loop over for binning grids with different shifts
    """

    H_bin = int(H_bin)
    if H_bin == 0: H_bin = 1
    K_bin = int(K_bin)
    if K_bin == 0: K_bin = 1
    L_bin = int(L_bin)
    if L_bin == 0: L_bin = 1
    print('H_bin:',H_bin)
    print('K_bin:',K_bin)
    print('L_bin:',L_bin,'\n')

    _H_step = float(H_step)/H_bin
    _H_shifts = np.round(np.arange(H_bin)*_H_step,4)#[1:]
    _K_step = float(K_step)/K_bin
    _K_shifts = np.round(np.arange(K_bin)*_K_step,4)#[1:]
    _L_step = float(L_step)/L_bin
    _L_shifts = np.round(np.arange(L_bin)*_L_step,4)#[1:]

    print('H_shifts:',_H_shifts)
    print('K_shifts:',_K_shifts)
    print('L_shifts:',_L_shifts)

    _H_shifts, _K_shifts, _L_shifts = np.meshgrid(_H_shifts,_K_shifts,_L_shifts,indexing='ij')
    _H_shifts = _H_shifts.flatten()
    _K_shifts = _K_shifts.flatten()
    _L_shifts = _L_shifts.flatten()
    shifts = np.c_[_H_shifts,_K_shifts,_L_shifts]

    ind = np.argwhere((shifts == 0.0).all(axis=1))
    shifts = np.delete(shifts,ind,axis=0)

    return np.atleast_2d(shifts)

# --------------------------------------------------------------------------------------------------

def _get_bins(lo,hi,bin,offset=0.0):
    """
    get binning args for MDNorm
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

def bin_NXSPE(event_files,goniometer_file,output_file,UB_params,H_lo,H_hi,H_step,K_lo,K_hi,K_step,
            L_lo,L_hi,L_step,H_bin=1,K_bin=1,L_bin=1):

    """
    wrapper to translate variable names for interface with PreprocessNXSPE. note, E-binning is 
    fixed by the binning already present in the *.nxspe files

    H_lo = lower bin center. 
    H_hi = *requested* upper bin center. if it's not commensurate with H_step, it will be tweaked!
    H_step = widths of the bins centered on the requested array of bin-centers.
    
        e.g. for H_lo = 0.0, H_hi = 1.0, H_step = 0.1, the bin centers will be H = [0.0, 0.1, 0.2,
        ..., 1.0] which are integrated from -0.05 to 0.05, 0.05 to 0.15, etc.

        e.g. for H_lo = 0.0, H_hi = 1.02, H_step = 0.1, you will get the same binning as above.

    K_* and L_* have the same meaning. 

    H_bin, K_bin, and L_bin are number of steps *within* a bin to take. data are integrated 
    once for each step by shifting bin centers by that amount. but keeping bin widths the same:
    e.g. if H_step = 0.1 and H_bin = 2, first bin centers will be shifted by 0.0 (i.e. unshifted),
    then shifted by 0.05 rlu and integrated again. shifts along different directions are 'meshed', 
    i.e. H_bin = 2, K_bin = 2, L_bin = 1, H_step = 0.1, K_step = 0.1, L_step = 0.1 will result in 
    4 integrations with shifts
        1) [0.00, 0.00, 0.00]
        2) [0.00, 0.05, 0.00]
        3) [0.05, 0.00, 0.00]
        4) [0.05, 0.05, 0.00]
        
    each shifted integration grid will be appended to the same file. e.g. for if H_bin = 2 and 
    H_lo = 0.0, H_hi = 1.0, and H_step = 0.1, the bin centers will be H = [0.00, 0.05, 0.10, 
    0.15, ..., 0.95, 1.00, 1.05] intergrated from -0.05 to 0.05, 0.00 to 0.10, 0.05 to 0.15, etc.

    output_file is the name where the output histogram data will be written.

    see the descriptions in the $ROOT/python_code/PreprocessNXSPE.py script for explanation of
    the other args to this function
    """

    timer = c_timer('bin_NXSPE',units='m')

    # convert binning args to ones used by c_MDE_tools
    H_bins = _get_bins(H_lo,H_hi,H_step)
    K_bins = _get_bins(K_lo,K_hi,K_step)
    L_bins = _get_bins(L_lo,L_hi,L_step)
    print('H_bins:',H_bins)
    print('K_bins:',K_bins)
    print('L_bins:',L_bins,'\n')

    # convert shift args to useful form
    shifts = _setup_shifts(H_step,K_step,L_step,H_bin,K_bin,L_bin)
    if shifts.size == 0:
        loop_over_shifts = False
    else:
        loop_over_shifts = True
    print('\nbin shifts:')
    print(shifts,'\n')

    # bin and save data with no offset
    MD_workspace = PreprocessNXSPE(H_bins,K_bins,L_bins,event_files,goniometer_file,UB_params)
    save_MDE_to_hdf5(MD_workspace,output_file)

    # only loop over offsets if atleast one is defined
    if loop_over_shifts:

        num_shifts = shifts.shape[0]

        print('\n----------------------------------------------------------------\n')
        print('num_shifts:',num_shifts)
        print('shifts:')
        print(shifts)

        msg = '\nlooping over shifts'
        print(msg)

        _H_shift = np.array([0,0,0],dtype=float)
        _K_shift = np.array([0,0,0],dtype=float)
        _L_shift = np.array([0,0,0],dtype=float)

        for ii in range(num_shifts):

            _H, _K, _L = shifts[ii,:]
            print(f'\nshifts[{ii+1}]:',_H,_K,_L,'\n')

            _H_shift[0] = _H; _H_shift[2] = _H
            _K_shift[0] = _K; _K_shift[2] = _K
            _L_shift[0] = _L; _L_shift[2] = _L

            # go and bin the data -- these MUST be appended.
            MD_workspace = PreprocessNXSPE(H_bins+_H_shift,K_bins+_K_shift,L_bins+_L_shift,
                                           event_files,goniometer_file,UB_params)
            save_MDE_to_hdf5(MD_workspace,output_file,append=True)

    timer.stop()

# --------------------------------------------------------------------------------------------------    


