#f2py -c -m readEnkf readEnkf.f90 readenkfc.o
from readEnkf import *
from numpy import *
orbnumber=1456
f=open('enkFileList','r')
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.colors as col

lines=f.readlines()
iprofs=0
r1=[]
r2=[]
piaSRTKu=[]
piaHB=[]
piadHB=[]
piaKuDF=[]
piadDF=[]
rrateAv=zeros((88),float)
rrateL=[]
zKuEnsL=[]
for l in lines[:]:
    l=l.split('File')
    oNumb=int(l[1])
    #oNumb=2665
    openfile(oNumb)
    for i in range(1000):
        nmemb1,raintype,node5,nzka,wfractm,freezh = read1()
    #print nzka
        if nmemb1==0:
            closefile()
            break

        if(nzka>0):
            zkuObs,izka,yEns,yObs,\
                zkuEns,log10dnw,rrate,\
                sfcRainEns,dsrtPIAku,dsrtPIAka,\
                pia13mod,pia35mod,\
                srtpiaku,srtrelpiaku,\
                dsrtrelpia,tbobs,tbsim,i1f,i2f = read2(nmemb1,node5,nzka)
        else:
            zkuObs,\
                zkuEns,log10dnw,rrate,\
                sfcRainEns,dsrtPIAku,dsrtPIAka,\
                pia13mod,pia35mod,\
                srtpiaku,srtrelpiaku,\
                dsrtrelpia,tbobs,tbsim,i1f,i2f = read3(nmemb1,node5,nzka)
        for k in range(45):
            rrate[node5[4]:,k]=rrate[node5[4],k]
            rrate[rrate<0]=0
#            print rrate.min()
        if freezh>4.2 and raintype==200 and node5[0]<node5[1]:
            rrateL.append(rrate.mean(axis=1))
            zKuEnsL.append(zkuObs.copy())
            iprofs=iprofs+1
            covXY=cov(sfcRainEns,yEns)
            yMean=yEns.mean(axis=1)
            n=yObs.shape[0]
            covyy=covXY[1:,1:]+eye(n)*4
            covxy=covXY[:1,1:]
            kgain=dot(covxy,linalg.pinv(covyy))
            dx=dot(kgain,yObs-yMean)
            r1.append(sfcRainEns.mean())
            r21=sfcRainEns.mean()+dx[0]
            if r21>sfcRainEns.max():
                r21=sfcRainEns.max()
            if r21<sfcRainEns.min():
                r21=sfcRainEns.min()
            r2.append(r21)
            rrateAv=rrateAv+rrate.mean(axis=1)
            if srtrelpiaku>3 and dsrtrelpia>3 and (i1f==1) and (i2f==1):
                #print i1f
                piaSRTKu.append(srtpiaku)
                piaKuDF.append(dsrtPIAku)
                piadDF.append(dsrtPIAka-dsrtPIAku)
                piaHB.append(pia13mod.mean())
                piadHB.append(pia35mod.mean()-pia13mod.mean())
        #print tbobs[9:]
        #print tbsim[9:,0]
    #print sfcRainEns
            
    
print corrcoef(piadHB,piadDF)
plt.plot(rrateAv/iprofs,(arange(88)*0.25)[::-1])
#plt.plot(rrateL[1],(arange(88)*0.25)[::-1])
#plt.plot(rrateL[150],(arange(88)*0.25)[::-1])
#plt.plot(zKuEnsL[150],(arange(88)*0.25)[::-1])
plt.ylim(0,15)

