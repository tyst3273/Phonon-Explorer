from TextFile import *
import numpy
from numpy import *
import math
import os
import random
import io
import matlab.engine
from RSE_Constants import *

#This class deals with generation of wavevectors used for background determination
class BackgroundQ:
    def __init__(self,h,k,l,params,index): #h,k,l for the dataset for which background is being calculated
        #self.Qslash gives the values of the wavevector (in r.l.u) used for background determination
        #index enumerates how many time it has been called already. It allows the module to return
        #different values each time.
        self.params=params
        self.H=h
        self.K=k
        self.L=l
        self.index=index
        self.mult=1
        self.flag=0
        #CalcQslash takes input in reciprocal angstoms, result returned in r.l.u
        self.Qslash=self.CalcQslash()
#        self.fileName=str("H%5.2f K%5.2f L%5.2f" % (self.Qslash[0],self.Qslash[1],self.Qslash[2]))
#  REPLACE this line with the following that is commented out
        if self.flag==0:
            self.fileName=str(RSE_Constants.FILENAME_FORMAT % (self.Qslash[0],self.Qslash[1],self.Qslash[2]))    

    def Qabs(self):
        Qx=self.H*2*math.pi/self.params.a
        Qz=self.K*2*math.pi/self.params.b
        Qy=self.L*2*math.pi/self.params.c
        A=sqrt(Qx**2+Qy**2+Qz**2)
        return A

    def CalcQslash(self):

#a-lattice param along proj.u, b- along proj.v

        if self.index>4:  #This statement controls when to stop generating Q's. In this 
            self.flag=1    #example it will stop after 4 times (3+1).
            Qslash=[0,0,0]
            return

        hs=[-6,-5,-6.5,-5,-6]
        ks=[-5,-4,-5.5,-5,-4]
            
        hslash=hs[self.index]
        kslash=ks[self.index]
        self.mult=(self.H**2+self.K**2)/(hslash**2+kslash**2)


        lslash=0
#        lslash=((self.Qabs()**2-(hslash*2*math.pi/self.params.a)**2-(kslash*2*math.pi/self.params.a)**2)**0.5)*self.params.c/(2*math.pi)

        Qslash=[hslash,kslash,lslash]

            
#        print(self.H,hslash,self.mult)
        
#        Qslash=[-4,0,7.5]
#        Qslash=[-5,0,1.5]
#        if self.H<-5:
#            Qslash=[-5,0,5.8]
#            Qslash=[-5.5,0,1.5]
#        Qs_s0=Qu**2+Qv**2+Qp**2
#        Qs_s=Qu1slash**2+Qv1slash**2+Qp1slash**2
#        print(math.sqrt(Qs_s))
#        print(math.sqrt(Qs_s0))
        return Qslash
