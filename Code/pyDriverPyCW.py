
from numpy import *
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.colors as col
from matplotlib.colors import BoundaryNorm
from pysimWrapKu import *
from pyHB import *
from readData import *
import SDSU as sdsu
sdsu.setparams()
sdsu.simradarlut(7)
nmfreq=8
ngates=76
nmemb1=30


nr1=0
nr2=49
n1=5700
n2=6000

f2AKu='/gpmdata/2014/03/18/radar/2A.GPM.Ku.V7-20170308.20140318-S205917-E223142.000299.V05A.HDF5'
f2ADPR='/gpmdata/2014/03/18/radar/2A.GPM.DPR.V7-20170308.20140318-S205917-E223142.000299.V05A.HDF5'
f2KuENV='/gpmdata/2014/03/18/radar/2A-ENV.GPM.Ku.V7-20170308.20140318-S205917-E223142.000299.V05A.HDF5'
from netCDF4 import Dataset
#f=Dataset(f2AKu,'r')
#stop
files=open('fileList','r').readlines()
nprof=0
dBase=[]
from setProf import *
nmfreqm=8
import pyHB as pyHB
initt()
piaS=[]
piaS2=[]
kzHP=array([ 1.02332128, -5.4866022 ])
import pickle
ifi=-1
from retrGaussN import *
for f in files[:]:
    ifi+=1
    print ifi
    cont=open(f[:-1],'r').readlines()
    for l1 in cont:
        l1s=l1.split('=')
        if l1s[0]=='f2AKu':
            f2AKu=l1s[1][:-1]
        if l1s[0]=='f2ADPR':
            f2ADPR=l1s[1][:-1]
        if l1s[0]=='f2KuEnv':
            f2KuEnv=l1s[1][:-1]
        if l1s[0]=='fCMB':
            f2BCMB=l1s[1][:-1]
#
    print f2AKu,f2ADPR
        
    i1=f2ADPR.find('2A.GPM.DPR.V7')
    s1='2B.GPM.DPRGMI.CORRA2016'
    s2=f2ADPR[i1:i1+22]
    f2BCMB=f2ADPR.replace(s2,s1)
    print f2BCMB
        #stop
    n1=00
    n2=9000
    #n1=5700
    #n2=6000
    nr1=12
    nr2=37
    pressEnv,tEnv,qvEnv,sfcTempEnv,sknTempEnv=readENV(f2KuENV,n1,n2,nr1,nr2)
    zFactKu,nodeK,flag,precipType,sfcType,zeroDeg,binClut,\
        lon1,lat1,zKu,hZero,piaKu,binBBBottom,binBBTop,\
        binBBPeak,piaHB,sfcRain,binSf,reliab,\
        locZAngle,s0Ku,eps,dm=readKuR(f2AKu,n1,n2,nr1,nr2)
    zFactKa, piaKa,piaKuMS,reliabFact,s0, zdpr,zdpr2,sfcRain2= readDPR2(f2ADPR,n1,n2,nr1,nr2)
    du1,du2,zMS,nubf,sfcRainNS,sfcRainMS=read2BCMB(f2BCMB,n1,n2,12,37)
        #stop
    a=nonzero(sfcType>0)
    b=nonzero((precipType[a]/1e7).astype(int)==2)
    c=nonzero(abs(reliab[a][b]-1.5)<2.6)
    print len(c[0])
    piaMS=du1[:,:,1]
    piaNS=du1[:,:,0]
    for i,j in zip(a[0][b][c],a[1][b][c]):
        zL=list(zFactKu[i,j,:])
        zL.append(zeroDeg[i,j])
        zL.append(piaKu[i,j])
        zL.append(binSf[i,j])
        zL.append(piaHB[i,j])
        if binClut[i,j]>=167 and  binSf[i,j]>172:
            dBase.append(zL)
            
        z13obs=log10(0.5*(10.**(0.1*zFactKu[i,j,0::2])+\
                              10**(0.1*zFactKu[i,j,1::2]))+1e-9)*10.
        z13ku=log10(0.5*(10.**(0.1*zKu[i,j,0::2])+\
                             10**(0.1*zKu[i,j,1::2]))+1e-9)*10.
        z35ka=log10(0.5*(10.**(0.1*zFactKa[i,j,0::2])+\
                             10**(0.1*zFactKa[i,j,1::2]))+1e-9)*10.
        
        imu,hh,hfreez,pia13srt,log10dnw,node,isurf,dr,itype=\
            setProf(binSf,binClut,precipType,zeroDeg,binBBPeak,piaKu,hZero,i,j,z13obs)
        node[node>88]=88
        node[node<1]=1
        imembc=1
        a1=nonzero(z13obs[node[0]:node[2]]>45)
        z13obsNoHail=z13obs.copy()
        zmix=10*log10(10**(z13obs[node[2]-20:node[2]-4]*0.1).mean())
        log10dnw=log10dnw*0.-0.5

        if(zmix<40 and (reliabFact[i,j]==1 or reliabFact[i,j]==2) and \
               not(reliab[i,j]==1 or reliab[i,j]==2) and \
               zMS[i,j]<5):# and i==234 and j==13):
            R=2.
            B=2.
            dpia=piaKa[i,j]-piaKuMS[i,j]
            z13,z35,pia13,pia35,z35mod,\
                lwc,kext,salb,asym,rrate,dm,\
                log10dnw,zeta,pia130,pia350,rrate0=gaussN(z13obs,log10dnw,dr,node,\
                                         isurf,imu,nmfreqm,hh,\
                                         itype,hfreez,pia13srt,\
                                         imembc,pyHB,dpia,R,B)
            im=1
            piaS2.append([pia35-pia13,pia350-pia130,\
                              piaKa[i,j]-piaKuMS[i,j],pia13,zMS[i,j],\
                              rrate[node[4]-1],rrate0[node[4]-1],\
                              sfcRainNS[i,j],sfcRainMS[i,j],\
                              sfcRain[i,j],sfcRain2[i,j],\
                              log10dnw[node[4]-1],dm[node[4]-1]])
            #print zmix,pia132
            if pia13>1000:
                    #print pia13, pia132
                if reliab[i,j]==1 or reliab[i,j]==2:
                    piaS.append([pia13,pia132,piaKu[i,j],zMS[i,j],nubf[i,j],piaKuMS[i,j],\
                                     piaKa[i,j],piaMS[i,j],piaNS[i,j],reliabFact[i,j],z13obs[node[4]-1],sfcRain[i,j],sfcRain2[i,j],\
                                     sfcRainNS[i,j],sfcRainMS[i,j],rrate[node[4]-1]])
                else:
                    piaS2.append([pia13,pia132,piaKu[i,j],zMS[i,j],nubf[i,j],piaKuMS[i,j],\
                                      piaKa[i,j],piaMS[i,j],piaNS[i,j],reliabFact[i,j],z13obs[node[4]-1],sfcRain[i,j],sfcRain2[i,j],\
                                      sfcRainNS[i,j],sfcRainMS[i,j],rrate[node[4]-1]])
                    
    print len(piaS2)
    nprof+=len(c[0])
    pickle.dump([piaS2],open('piaSCW3.pklz','wb'))
    #print corrcoef(array(piaS).T)[-5:,-5:]
    #print corrcoef(array(piaS2).T)[-5:,-5:]
    print corrcoef(array(piaS2).T)[:5,:5]
    print corrcoef(array(piaS2).T)[5:11,5:11]
#
#zFactKa2, piaKa2,relibFact2,s02= read2AKa(f2AKa,n1,n2)



#    print corrcoef(array(piaS2).T)[-7:-3,-7:-3]
#print corrcoef(array(piaS2).T)[-7:-3,-7:-3]
#
