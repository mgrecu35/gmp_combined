
from numpy import *
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.colors as col
from matplotlib.colors import BoundaryNorm
#from pysimWrapKu import *
from pyHB import *
from readData import *
#import SDSU as sdsu
#sdsu.setparams()
#sdsu.simradarlut(7)
nmfreq=8
ngates=76
nmemb1=30


nr1=0
nr2=49
n1=900
n2=1200
initt()
f2AKu='DPRnew/GPMCOR_DPR_1802210352_0525_022628_L2S_DD2_05A.h5'
zFactKu,zFactKa,zFactKaHS,\
    nodeK,flag,precipType,sfcType,zeroDeg,binClut,\
    lon1,lat1,zKu,hZero,piaKu,binBBBottom,binBBTop,\
    binBBPeak,piaHB,sfcRain,binSf,reliab,\
    locZAngle,s0Ku,eps,dm=readDPR(f2AKu,n1,n2,nr1,nr2)
            
z13obs=log10(0.5*(10.**(0.1*zFactKu[i,j,0::2])+\
                  10**(0.1*zFactKu[i,j,1::2]))+1e-9)*10.
z13ku=log10(0.5*(10.**(0.1*zKu[i,j,0::2])+\
                 10**(0.1*zKu[i,j,1::2]))+1e-9)*10.
z35ka=log10(0.5*(10.**(0.1*zFactKa[i,j,0::2])+\
                 10**(0.1*zFactKa[i,j,1::2]))+1e-9)*10.

imu,hh,hfreez,pia13srt,log10dnw,\
    node,isurf,dr,itype=\
                         setProf(binSf,binClut,precipType,\
                                 zeroDeg,binBBPeak,piaKu,hZero,i,j,z13obs)
