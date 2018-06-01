import h5py as h5

def readKuR(fnamein,n1,n2,nr1,nr2):
    f=h5.File(fnamein,'r')
    #print f['NS/PRE/zFactorMeasured'].shape
    zFactKu=f['NS/PRE/zFactorMeasured'][n1:n2,nr1:nr2,:]
    node=f['NS/DSD/binNode'][n1:n2,nr1:nr2,:]
    flag=f['/NS/PRE/flagPrecip'][n1:n2,nr1:nr2]
    precipType=f['/NS/CSF/typePrecip'][n1:n2,nr1:nr2]
    sfcType=f['/NS/PRE/landSurfaceType'][n1:n2,nr1:nr2]
    zeroDeg=f['/NS/VER/binZeroDeg'][n1:n2,nr1:nr2]
    binClut=f['/NS/PRE/binClutterFreeBottom'][n1:n2,nr1:nr2]
    binSf=f['/NS/PRE/binRealSurface'][n1:n2,nr1:nr2]
    binBBPeak=f['/NS/CSF/binBBPeak'][n1:n2,nr1:nr2]
    binBBBottom=f['/NS/CSF/binBBBottom'][n1:n2,nr1:nr2]
    binBBTop=f['/NS/CSF/binBBTop'][n1:n2,nr1:nr2]
    lon1=f['/NS/Longitude'][n1:n2,nr1:nr2]
    lat1=f['/NS/Latitude'][n1:n2,nr1:nr2]
    zKu=f['/NS/SLV/zFactorCorrected'][n1:n2,nr1:nr2,:]
    hZero=f['/NS/VER/heightZeroDeg'][n1:n2,nr1:nr2]
    piaKu=f['/NS/SRT/pathAtten'][n1:n2,nr1:nr2]
    piaHB=f['/NS/SLV/piaFinal'][n1:n2,nr1:nr2]
    sfcRain=f['NS/SLV/precipRateNearSurface'][n1:n2,nr1:nr2]
    eps=f['NS/SLV/epsilon'][n1:n2,nr1:nr2]
    reliabF=f['/NS/SRT/reliabFlag'][n1:n2,nr1:nr2]
    locZAngle=f['/NS/PRE/localZenithAngle'][n1:n2,nr1:nr2]
    s0=f['NS/PRE/sigmaZeroMeasured'][n1:n2,:]
    dm=f['NS/SLV/paramDSD'][n1:n2,nr1:nr2,:,:]
    f.close()
    return zFactKu,node,flag,precipType,sfcType,zeroDeg,binClut,\
        lon1,lat1,zKu,hZero,piaKu,binBBBottom,binBBTop,\
        binBBPeak,piaHB,sfcRain,binSf,reliabF,locZAngle,s0,eps,dm

def readDPR(fnamein,n1,n2,nr1,nr2):
    f=h5.File(fnamein,'r')
    #print f['NS/PRE/zFactorMeasured'].shape
    zFactKu=f['NS/PRE/zFactorMeasured'][n1:n2,nr1:nr2,:]
    zFactKa=f['MS/PRE/zFactorMeasured'][n1:n2,:,:]
    zFactKaHS=f['MS/PRE/zFactorMeasured'][n1:n2,:,:]
    node=f['NS/DSD/binNode'][n1:n2,nr1:nr2,:]
    flag=f['/NS/PRE/flagPrecip'][n1:n2,nr1:nr2]
    precipType=f['/NS/CSF/typePrecip'][n1:n2,nr1:nr2]
    sfcType=f['/NS/PRE/landSurfaceType'][n1:n2,nr1:nr2]
    zeroDeg=f['/NS/VER/binZeroDeg'][n1:n2,nr1:nr2]
    binClut=f['/NS/PRE/binClutterFreeBottom'][n1:n2,nr1:nr2]
    binSf=f['/NS/PRE/binRealSurface'][n1:n2,nr1:nr2]
    binBBPeak=f['/NS/CSF/binBBPeak'][n1:n2,nr1:nr2]
    binBBBottom=f['/NS/CSF/binBBBottom'][n1:n2,nr1:nr2]
    binBBTop=f['/NS/CSF/binBBTop'][n1:n2,nr1:nr2]
    lon1=f['/NS/Longitude'][n1:n2,nr1:nr2]
    lat1=f['/NS/Latitude'][n1:n2,nr1:nr2]
    zKu=f['/NS/SLV/zFactorCorrected'][n1:n2,nr1:nr2,:]
    hZero=f['/NS/VER/heightZeroDeg'][n1:n2,nr1:nr2]
    piaKu=f['/NS/SRT/pathAtten'][n1:n2,nr1:nr2]
    piaHB=f['/NS/SLV/piaFinal'][n1:n2,nr1:nr2]
    sfcRain=f['NS/SLV/precipRateNearSurface'][n1:n2,nr1:nr2]
    eps=f['NS/SLV/epsilon'][n1:n2,nr1:nr2]
    reliabF=f['/NS/SRT/reliabFlag'][n1:n2,nr1:nr2]
    locZAngle=f['/NS/PRE/localZenithAngle'][n1:n2,nr1:nr2]
    s0=f['NS/PRE/sigmaZeroMeasured'][n1:n2,:]
    dm=f['NS/SLV/paramDSD'][n1:n2,nr1:nr2,:,:]
    f.close()
    return zFactKu,zFactKa,zFactKaHS,node,flag,precipType,\
        sfcType,zeroDeg,binClut,\
        lon1,lat1,zKu,hZero,piaKu,binBBBottom,binBBTop,\
        binBBPeak,piaHB,sfcRain,binSf,reliabF,locZAngle,s0,eps,dm

def readKu(fnamein,n1,n2,nr1,nr2):
    f=h5.File(fnamein,'r')
    #print f['NS/PRE/zFactorMeasured'].shape
    zFactKu=f['NS/PRE/zFactorMeasured'][n1:n2,nr1:nr2,:]
    node=f['NS/DSD/binNode'][n1:n2,nr1:nr2,:]
    flag=f['/NS/PRE/flagPrecip'][n1:n2,nr1:nr2]
    precipType=f['/NS/CSF/typePrecip'][n1:n2,nr1:nr2]
    sfcType=f['/NS/PRE/landSurfaceType'][n1:n2,nr1:nr2]
    zeroDeg=f['/NS/VER/binZeroDeg'][n1:n2,nr1:nr2]
    binClut=f['/NS/PRE/binClutterFreeBottom'][n1:n2,nr1:nr2]
    binSf=f['/NS/PRE/binRealSurface'][n1:n2,nr1:nr2]
    binBBPeak=f['/NS/CSF/binBBPeak'][n1:n2,nr1:nr2]
    binBBBottom=f['/NS/CSF/binBBBottom'][n1:n2,nr1:nr2]
    binBBTop=f['/NS/CSF/binBBTop'][n1:n2,nr1:nr2]
    lon1=f['/NS/Longitude'][n1:n2,nr1:nr2]
    lat1=f['/NS/Latitude'][n1:n2,nr1:nr2]
    zKu=f['/NS/SLV/zFactorCorrected'][n1:n2,nr1:nr2,:]
    p3d=f['/NS/SLV/precipRate'][n1:n2,nr1:nr2,:]
    hZero=f['/NS/VER/heightZeroDeg'][n1:n2,nr1:nr2]
    piaKu=f['/NS/SRT/pathAtten'][n1:n2,nr1:nr2]
    piaHB=f['/NS/SLV/piaFinal'][n1:n2,nr1:nr2]
    sfcRain=f['NS/SLV/precipRateNearSurface'][n1:n2,nr1:nr2]
    eps=f['NS/SLV/epsilon'][n1:n2,nr1:nr2]
    reliabF=f['/NS/SRT/reliabFlag'][n1:n2,nr1:nr2]
    locZAngle=f['/NS/PRE/localZenithAngle'][n1:n2,nr1:nr2]
    fHeader=f.attrs['FileHeader']
    s0=f['NS/PRE/sigmaZeroMeasured'][n1:n2,:]
    f.close()
    return zFactKu,node,flag,precipType,sfcType,zeroDeg,binClut,\
        lon1,lat1,zKu,hZero,piaKu,binBBBottom,binBBTop,\
        binBBPeak,p3d,piaHB,sfcRain,fHeader,binSf,reliabF,locZAngle,s0,eps

def readDPR2(fnamein,n1,n2,nr1,nr2):
    f=h5.File(fnamein,'r')
    zFactKa=f['MS/PRE/zFactorMeasured'][n1:n2,:,:]
    zdpr=f['NS/SLV/zFactorCorrected'][n1:n2,nr1:nr2,:] 
    zdpr2=f['MS/SLV/zFactorCorrected'][n1:n2,:,:]
    piaKa=f['/MS/SRT/pathAtten'][n1:n2,:]
    piaKuMS=f['/NS/SRT/pathAtten'][n1:n2,12:37]
    relibFact=f['/MS/SRT/reliabFlag'][n1:n2,:]
    s0=f['MS/PRE/sigmaZeroMeasured'][n1:n2,:]
    sfcRain2=f['MS/SLV/precipRateNearSurface'][n1:n2,:]
    f.close()
    return zFactKa, piaKa,piaKuMS,relibFact,s0,zdpr,zdpr2,sfcRain2

def read2AKa(fnamein,n1,n2):
    f=h5.File(fnamein,'r')
    zFactKa=f['MS/PRE/zFactorMeasured'][n1:n2,:,:]
    piaKa=f['/MS/SRT/pathAtten'][n1:n2,:]
    relibFact=f['/MS/SRT/reliabFlag'][n1:n2,:]
    s0=f['MS/PRE/sigmaZeroMeasured'][n1:n2,:]
    f.close()
    return zFactKa, piaKa,relibFact,s0

def read2BCMB(fnamein,n1,n2,nr1,nr2):
    f=h5.File(fnamein,'r')
    piaMS=f['/MS/pia'][n1:n2,:,:]
    piaNS=f['/NS/pia'][n1:n2,nr1:nr2]
    sfcRainNS=f['/NS/surfPrecipTotRate'][n1:n2,nr1:nr2]        
    sfcRainMS=f['/MS/surfPrecipTotRate'][n1:n2,:]        
    zMS=f['MS/multiScatMaxContrib'][n1:n2,:]
    nubf=f['MS/nubfPIAfactor'][n1:n2,:]
    f.close()
    return piaMS,piaNS,zMS,nubf,sfcRainNS,sfcRainMS

def readENV(fname,n1,n2,nr1,nr2):
    f=h5.File(fname,'r')
    pressEnv=f['NS/VERENV/airPressure'][n1:n2,nr1:nr2]
    tEnv=f['NS/VERENV/airTemperature'][n1:n2,nr1:nr2,:]
    qvEnv=f['NS/VERENV/waterVapor'][n1:n2,nr1:nr2,:,:]
    sfcTempEnv=f['NS/VERENV/surfaceTemperature'][n1:n2,nr1:nr2]
    skTempEnv=f['NS/VERENV/skinTemperature'][n1:n2,nr1:nr2]
    return pressEnv,tEnv,qvEnv,sfcTempEnv,skTempEnv

