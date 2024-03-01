################################
#
# Copy Facilities.xml to your 
# .mantid/instrument folder
#
################################


#from file_tools.m_save_MDE_to_hdf5 import save_MDE_to_hdf5
#from mantid.simpleapi import *

import glob
import numpy as np
import math


#EDIT HERE

stepH=0.1
stepK=0.2
stepL=0.2
HMult=2 #number of steps in each bin (This way bins can overlap)
KMult=2
LMult=2
HStart=1
HEnd=2
KStart=1
KEnd=2
LStart=1
LEnd=2
###############################

NHSteps=math.ceil((HEnd-HStart)/stepH)
CorrectHEnd=HStart+NHSteps*stepH

KStart=-6
KEnd=6
LStart=-4
LEnd=4

#HBin=zeros[Hmult][3]
HBin = np.zeros((HMult,3))
for i in range(HMult):
    HBin[i,:]=[HStart-stepH/2+i*stepH/HMult,stepH,CorrectHEnd+stepH/2+i*stepH/HMult]

print(HBin)

HBin=[1-stepH/2,stepH,2+stepH/2]
KBin=[-6.1-stepK/2,stepK,6.1+stepK/2]
LBin=[-4-stepL/2,stepL,4+stepL/2]
files = sorted(glob.glob('/SNS/ARCS/IPTS-26347/shared/jpark/nxspe120Ei_010K_Ebin1_RCmasked/S0*.nxspe'))

#This is the file that lists the goniometer angler for each NXPCE file
fn = '/SNS/ARCS/IPTS-26347/shared/jpark/120meV300Hz010K_run_list.txt'
hdf5_file_name = 'nxspe120Ei_010K_Ebin1ZBk.hdf5'

#This should be aligned UV matrix
UBpars={'a':3.93,'b':3.93,'c':3.93,'alpha':90,'beta':90,'gamma':90,'u':'1,0,-0.1','v':'0,1,-0.1'}

########################################


d=np.loadtxt(fn,skiprows=1,dtype=str)
runs=d[:,0]
gon=d[:,11].astype(float)-82.388
gon_dict=dict()
for r,g in zip(runs, gon):
    gon_dict[r]=g

wg=Load(','.join(files))
for wsn in wg.getNames():
    ang=gon_dict[wsn.split('S0')[1]]
    SetGoniometer(Workspace=wsn, Axis0 = f'{ang},0,1,0,1')
axis_deltaE=wg[0].readX(0)
emin=axis_deltaE[0]
emax=axis_deltaE[-1]
de_step=axis_deltaE[1]-emin
ConvertFromDistribution(wg)
wge=ConvertToEventWorkspace(wg)
wge=CropWorkspaceForMDNorm(wge,XMin=emin,XMax=emax)
SetUB(wge,**UBpars)
AddSampleLog(Workspace='wge', LogName='gd_prtn_chrg', LogText='1.0', LogType='Number', NumberType='Double')
mdparts=ConvertToMD(InputWorkspace=wge, QDimensions='Q3D', dEAnalysisMode='Direct', Q3DFrames='Q_sample')
mde=MergeMD(mdparts)
#SaveMD(InputWorkspace="mde",Filename="/SNS/ARCS/IPTS-26347/shared/jpark/MDEventFile.nxs")

#UBpars={'a':3.93,'b':3.93,'c':3.93,'alpha':90,'beta':90,'gamma':90,'u':'1,0,-0.1','v':'0,1,-0.1'}
#SetUB(mde,**UBpars)

#print(axis_deltaE)
# test MDNorm - elastic slice
MDNorm(InputWorkspace=mde,
       QDimension0='1,0,0',
       QDimension1='0,1,0',
       QDimension2='0,0,1',
       Dimension0Name='QDimension0',
       Dimension0Binning="{},{},{}".format(HBin[0],HBin[1],HBin[2]),
       Dimension1Name='QDimension1',
       Dimension1Binning="{},{},{}".format(KBin[0],KBin[1],KBin[2]),
       Dimension2Name='QDimension2',
       Dimension2Binning="{},{},{}".format(LBin[0],LBin[1],LBin[2]),
       Dimension3Name='DeltaE',
       Dimension3Binning="{}".format(de_step),
#       SymmetryOperations='x,y,z;-x,-y,z;x,y,-z;-x,-y,-z;y,x,z;y,x,-z',
       OutputWorkspace='o',
       OutputDataWorkspace='d',
       OutputNormalizationWorkspace='n')

#SaveMD(InputWorkspace="o",Filename="/SNS/ARCS/IPTS-26347/shared/jpark/MDHistoFile.nxs")

MD_workspace = mtd['o']
save_MDE_to_hdf5(MD_workspace,hdf5_file_name)





