fComb='/gpmdata/2015/06/01/radar/2B.GPM.DPRGMI.CORRA2016.20150601-S065913-E083144.007136.V05A.HDF5'

#fComb='/gpmdata/2015/06/01/radar/2B.GPM.DPRGMI.CORRA2016.20150601-S130925-E144157.007140.V05A.HDF5'
fComb='/gpmdata/2015/07/12/radar/2B.GPM.DPRGMI.CORRA2016.20150712-S084637-E101909.007775.V05A.HDF5'
fComb='/gpmdata/2015/07/12/radar/2B.GPM.DPRGMI.CORRA2016.20150712-S071405-E084636.007774.V05A.HDF5'
fComb='/gpmdata/2015/12/12/radar/2B.GPM.DPRGMI.CORRA2016.20151212-S112350-E125623.010157.V05A.HDF5'#lighy
#fComb='/gpmdata/2015/12/13/radar/2B.GPM.DPRGMI.CORRA2016.20151213-S103224-E120457.010172.V05A.HDF5'
f1CGMI='/gpmdata/2015/12/12/1C/1C.GPM.GMI.XCAL2016-C.20151212-S112350-E125623.010157.V05A.HDF5'
import cPickle as pickle
gIR='merg_2015071208_4km-pixel.nc4'
gIR='merg_2015121311_4km-pixel.nc4'
f2AKu='/gpmdata/2015/12/22/radar/2A.GPM.Ku.V7-20170308.20151222-S011552-E024824.010306.V05A.HDF5'
f2AKu='/itedata/ITE601/2015/12/22/radar/2A.GPM.Ku.V8-20180302.20151222-S011552-E024824.010306.ITE601.HDF5'
import cPickle as pickle
zMean=pickle.load(open('zMean.pklz','rb'))

def readcmbll(fname):
    f=h5.File(fname,'r')
    Lat=f['NS']['Latitude'][:,:]
    Lon=f['NS']['Longitude'][:,:]
    return Lat,Lon,f

from mpl_toolkits.basemap import Basemap, shiftgrid, cm
import glob
import datetime
import h5py as h5
from numpy import *
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.colors as col
from netCDF4 import Dataset

f1c=Dataset(f2AKu,'r')
ns=3161
fh=f1c.groups['NS'].groups['PRE']
#3150 shallow over mountain
#3030 over oceans
#6140 over Africa
#54-64 southern ocean
#4594-4596 west of ireland
x=array([arange(49) for i in range(176)]).T 
x2=array([arange(12,37) for i in range(176)]).T 
for ns in range(3025,-3034):
    zm=fh.variables['zFactorMeasured'][ns,:,:]
    binClutF=fh.variables['binClutterFreeBottom'][ns,:]
    bin0=f1c.groups['NS'].groups['VER'].variables['binZeroDeg'][ns,:]
    zmc=zm.copy()
    zAngle=f1c.groups['NS'].groups['PRE'].variables['localZenithAngle'][ns,:]
    h2d=array([(arange(176)*0.125*cos(zAngle[i]/180*pi))[::-1] \
                  for i in range(49)])
    
    for i in range(49):
        dc=10.**(0.1*zmc[i,150:176])
        dc[dc<0.1]=0.1
        zmc[i,150:176]=log10(dc)*10.
        zmc[i,binClutF[i]:]=-99
    zmc=ma.array(zmc,mask=zmc<-90)
    plt.pcolormesh(x,h2d,zmc[:,::],cmap='jet',vmin=12,vmax=45)
    plt.plot(arange(49),(176-bin0)*cos(zAngle/180.*pi)*0.125)
    plt.ylim(0,8)
    plt.xlabel('Ray',fontsize=13)
    plt.ylabel('Height (km)',fontsize=13)
    plt.title('Orbit 10306, Scan %4.4i \n 2AKu Measured Reflectivity'%ns,fontsize=14)
    ax=plt.colorbar()
    ax.set_label('dBZ')
    plt.savefig('orbit10306_Zm%i.png'%ns)
    plt.show()

f2AKu='/gpmdata/2015/12/14/radar/2A.GPM.Ku.V7-20170308.20151214-S080822-E094056.010186.V05A.HDF5'
scan=4747
f2AKu='/gpmdata/2015/12/14/radar/2A.GPM.Ku.V7-20170308.20151214-S094057-E111330.010187.V05A.HDF5'
f2AKu='/itedata/ITE601/2015/12/14/radar/2A.GPM.Ku.V8-20180302.20151214-S002531-E015805.010181.ITE601.HDF5'
f2AKu='/gpmdata/2017/04/01/radar/2A.GPM.DPR.V7-20170308.20170401-S051548-E064822.017558.V05A.HDF5' #Boston Snow
#f2AKu='/gpmdata/2017/10/29/radar/2A.GPM.Ku.V7-20170308.20171029-S003019-E020252.020837.V05A.HDF5' #Pensilvania Light Rain
#f2AKu='/gpmdata/2017/12/05/radar/2A.GPM.DPR.V7-20170308.20171205-S133506-E150740.021421.V05A.HDF5'
f2AKu='/gpmdata/2014/06/11/radar/2A.GPM.Ku.V7-20170308.20140611-S171129-E184401.001619.V05A.HDF5' 
f2AKu='/gpmdata/2014/07/09/radar/2A.GPM.Ku.V7-20170308.20140709-S113622-E130853.002051.V05A.HDF5'
f2AKu='/itedata/ITE601/2014/08/23/radar/2A.GPM.DPR.V8-20180305.20140823-S221559-E234832.002758.ITE601.HDF5'
f2AKu='DPRnew/GPMCOR_DPR_1802210352_0525_022628_L2S_DD2_05A.h5'
#f2AKu='DPRnew/GPMCOR_DPR_1709280352_0525_020357_L2S_DD2_05A.h5'#GPMCOR_DPR_1709280047_0220_020355_L2S_DD2_05A.h5'
f1c=Dataset(f2AKu,'r')
fh=f1c.groups['NS'].groups['PRE']

#for scan in range(4931,4971):
#for scan in range(4830,4850):# for 07/09
#for scan in range(4750,4760):# for 08/23

scan=422
for iscan in range(scan-2,-scan-1):
    LatHS=f1c.groups['HS'].variables['Latitude'][iscan,:]  
    LonHS=f1c.groups['HS'].variables['Longitude'][iscan,:]       
    plt.plot(LonHS[:12],LatHS[:12],'+')
    LatHS=f1c.groups['HS'].variables['Latitude'][iscan-1,:]  
    LonHS=f1c.groups['HS'].variables['Longitude'][iscan-1,:]  
    plt.plot(LonHS[12:],LatHS[12:],'+')
    LatNS=f1c.groups['NS'].variables['Latitude'][iscan,:]  
    LonNS=f1c.groups['NS'].variables['Longitude'][iscan,:]       
    #plt.plot(LonNS,LatNS)

#for scan in range(2170,2180):#new DPR
for scan in range(1000,1025):#new DPR
    plt.figure(figsize=(8,8))
    zm=fh.variables['zFactorMeasured'][scan,:,:]
    #zm=f1c.groups['NS'].groups['SLV'].variables['zFactorCorrected'][scan,:,:]
    binClutF=fh.variables['binClutterFreeBottom'][scan,:]
    Lat=f1c.groups['NS'].variables['Latitude'][scan,:]
    x=array([range(49) for i in range(176)]).T 
    x2=array([range(12,37) for i in range(176)]).T 

    binClutFMS=f1c.groups['MS'].groups['PRE'].variables['binClutterFreeBottom'][scan,:]
    binSf=fh.variables['binRealSurface'][scan,:]
    bin0=f1c.groups['NS'].groups['VER'].variables['binZeroDeg'][scan,:]
    binMS=f1c.groups['MS'].groups['VER'].variables['binZeroDeg'][scan,:]
    zmKa=f1c.groups['MS'].groups['PRE'].variables['zFactorMeasured'][scan,:,:]
    zmKaHS=f1c.groups['HS'].groups['PRE'].variables['zFactorMeasured'][scan-1,:,:]
    zmKaHS[:12,:]=zmKaHS[12:,:]
    zmKaHS[12:,:]=f1c.groups['HS'].groups['PRE'].variables['zFactorMeasured'][scan,:12,:][:,:]
    zAngle=f1c.groups['NS'].groups['PRE'].variables['localZenithAngle'][scan,:]
    LatHS=f1c.groups['HS'].variables['Latitude'][scan,:]
    h2d=array([(arange(176)*0.125*cos(zAngle[i]/180*pi))[::-1] \
                  for i in range(49)])
    h2d2=array([(arange(176)*0.125*cos(zAngle[i]/180*pi))[::-1] \
                   for i in range(25)])
    zmc=zm.copy()
    #x3=array([LatHS for i in range(176/2)]).T 
    for i in range(49):
        zmc[i,:]=zm[i,:]
        #zmc[i,binClutF[i]:]=-99
    ax=plt.subplot(311)
    zmc=ma.array(zmc,mask=zmc<12)
    Lat=range(49)
    plt.pcolormesh(x,h2d,zmc[:,::],cmap='jet',vmin=12,vmax=45)
    plt.plot(Lat,(175-bin0)*cos(zAngle/180.*pi)*0.125)
    plt.plot(Lat,(175-binSf)*cos(zAngle/180.*pi)*0.125)
    plt.plot(Lat,(175-binClutF)*cos(zAngle/180.*pi)*0.125)
    ax.axes.get_xaxis().set_visible(False)
    plt.ylim(0,10)
    plt.xlabel('Ray',fontsize=13)
    plt.title('Ku NS')
    plt.ylabel('Height (km)',fontsize=13)
    ax=plt.colorbar()
    ax.set_label('dBZ')
    zmc=zmKa.copy()
    zmc2=zmKaHS.copy()
    for i in range(25):
        zmc[i,:]=zmKa[i,:]
        #zmc[i,binClutFMS[i]:]=-99
    ax2=plt.subplot(312)
    zmc=ma.array(zmc,mask=zmc<12)

    plt.pcolormesh(x2,h2d[12:37,::],zmc[:,::],cmap='jet',vmin=12,vmax=45)
    plt.plot(Lat,(175-bin0)*cos(zAngle/180.*pi)*0.125)
    plt.plot(Lat,(175-binSf)*cos(zAngle/180.*pi)*0.125)
    plt.plot(Lat,(175-binClutF)*cos(zAngle/180.*pi)*0.125)
    plt.ylim(0,10)
    plt.xlabel('Ray',fontsize=13)
    plt.ylabel('Height (km)',fontsize=13)
    ax2.axes.get_xaxis().set_visible(False)
    plt.title('Ka MS')

    plt.xlim(0,48)
    ax=plt.colorbar()
    ax.set_label('dBZ')
    plt.subplot(313)
    lat2=[]
    zmL=[]
    ic=0
    for i in range(49):
        if (i-12)*(i-36)> 0:
            lat2.append(LatHS[ic])
            zmL.append(zmKaHS[ic,:])
            ic+=1
            print ic, i
        else:
            lat2.append(Lat[i])
            zmL.append([-99 for k in range(88)])
    
    lat2=array(lat2)
    lat2=arange(49)
    x3=array([lat2 for i in range(176/2)]).T 
    zmL=array(zmL)
    zmc2=ma.array(zmL,mask=zmL<12)

    plt.pcolormesh(x3,h2d[:,::2],zmc2[:,::],cmap='jet',vmin=12,vmax=45)
    plt.plot(Lat,(175-bin0)*cos(zAngle/180.*pi)*0.125)
    plt.plot(Lat,(175-binSf)*cos(zAngle/180.*pi)*0.125)
    plt.plot(Lat,(175-binClutF)*cos(zAngle/180.*pi)*0.125)
    plt.title('Ka HS')
    plt.ylim(0,10)
    plt.xlabel('Ray',fontsize=13)
    plt.ylabel('Height (km)',fontsize=13)
    plt.suptitle('Orbit 22628, Scan %4.4i'%scan,x=0.45,fontsize=14)
    ax=plt.colorbar()
    ax.set_label('dBZ')
    #stop
    plt.savefig('orbit22628_Zm%i.png'%scan)
    #plt.show()
