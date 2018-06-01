from numpy import *
def gaussN(z13obs,log10dnw,dr,node,\
               isurf,imu,nmfreqm,hh,\
               itype,hfreez,pia13srt,\
               imembc,pyHB,dpia,R,B):


    n=log10dnw.shape[0]
    dnw=interp((88),[0,node[1],node[3],87],[0.,0.,1.,1.])
    x=0.
    for it in range(4):
        z13,z35,pia13,pia35,z35mod,\
            lwc,kext,salb,asym,rrate,\
            dm,log10dnw,zeta = pyHB.fhb1r(z13obs,log10dnw+x*dnw,dr,node,\
                                              isurf,imu,nmfreqm,hh,\
                                              itype,hfreez,pia13srt,\
                                              imembc)
        if it==0:
            pia130=pia13
            pia350=pia35
            rrate0=rrate.copy()
        
        z131,z351,pia131,pia351,z35mod1,\
            lwc1,kext1,salb1,asym1,rrate1,dm1,\
            log10dnw,zeta =pyHB.fhb1r(z13obs,\
                                          log10dnw+(x+0.1)*dnw,dr,node,\
                                          isurf,imu,nmfreqm,hh,\
                                          itype,hfreez,pia13srt,\
                                          imembc)
        H=((pia351-pia131)-(pia35-pia13))/0.1
        A=(H/R*H+1/B)
        x=x+0.75/A*(H/R*(dpia-(pia35-pia13))-1/B*x)
        #print (dpia-(pia35-pia13)),pia131,pia13,\
        #    .25/A*(H/R*(dpia-(pia35-pia13))-1/B*x)
    #print x.25/A*(H/R*(dpia-(pia35-pia13))-1/B*x)
    if(x>4):
        x=4
    if(x<-4):
        x=-4
    z13,z35,pia13,pia35,z35mod,\
        lwc,kext,salb,asym,rrate,dm,\
        log10dnw,zeta =pyHB.fhb1r(z13obs,\
                                      log10dnw+x*dnw,dr,node,\
                                      isurf,imu,nmfreqm,hh,\
                                      itype,hfreez,pia13srt,\
                                      imembc)
    log10dnw+=x*dnw
    #print
    return z13,z35,pia13,pia35,z35mod,\
        lwc,kext,salb,asym,rrate,dm,\
        log10dnw,zeta,pia130,pia350,rrate0
