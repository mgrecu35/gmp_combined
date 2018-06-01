
! this is the Hitschfeld Bordan solution
! the variables in the command line are defined in 
! subroutine fModelFortran
! the formulae used in the solution are those in Grecu et al. 2011
 


subroutine fhb11(z13,z35,z13obs,&
     pia13,pia35,ic,jc,z35mod,lwc,log10dNw,dr,node,isurf,imu,&
     ngates,nmfreqm,hh,itype,kext,salb,asym,rrate,d0,hfreez,&
     pia13srt,imembC)
  !use cldclass
  use nbinMod
!  use geophysEns
  use ran_mod
  use gEnv  !! get the reliabFlag from here
  implicit none
  integer :: isurf
  integer :: ngates,node(5),imu(ngates), npart, nmfreqm
  real :: z13obs(ngates), z35mod(ngates), hh(ngates), dr
  real :: d0(ngates), rrate(ngates)
  real :: z13(ngates), z35(ngates),log10dNw(ngates), lwc(ngates)
  real :: pia13, pia35,zeta(ngates)
  real ::  alpha, beta, q, dzeta, dzetaOld, dn, hfreez
  integer :: i, j, k, i0, imembC, iEx
!begin  WSO declare ipix, jpix
  integer :: ipix, jpix



  real :: kext(ngates,nmfreqm), salb(ngates,nmfreqm),&
       asym(ngates,nmfreqm), z35ret(ngates)
  real :: dpia13, dpia35
  real :: pia13srt,pia35srt, ppia
  integer :: ic, jc, nodeA, itype
  real :: errpia35, dpia13n, pia13m, pia35u, dn2, rms, pia13v
!  real :: cldw(nlayer), rh(nlayer), wv_extMemb(nlayer), cld_extMemb(nlayer)
  real :: tavg
!  real :: rhPCij(nRhEofs), cldwPCij(nCldwEofs), 
  real :: rv(50)
  integer :: it, isinf, iNoAd
  real :: dz1,dz2
  rv=(/8.33649514e-04,   2.02394697e-02,   8.41044138e-01,&
       3.05986222e-01,   8.82082063e-01,   7.70981041e-01,&
       9.72449848e-01,   5.80631496e-01,   2.67921518e-01,&
       5.96091837e-01,   5.90155945e-01,   1.02158191e-02,&
       1.11711326e-01,   7.95037973e-01,   3.59357039e-01,&
       3.95581051e-01,   9.73541653e-01,   9.72111604e-01,&
       1.58652395e-01,   4.10102994e-01,   7.26038953e-01,&
       9.71678367e-01,   5.58166615e-01,   6.00756162e-01,&
       2.88898799e-01,   2.53188393e-01,   9.88474634e-02,&
       7.63488856e-01,   9.63377447e-01,   1.75159870e-01,&
       4.44191557e-01,   5.20691784e-01,   8.24856853e-01,&
       4.73616222e-01,   2.15273524e-01,   6.14520925e-01,&
       1.19503728e-01,   8.87660626e-01,   3.50617762e-01,&
       5.23541528e-01,   5.17400886e-01,   2.22747951e-01,&
       1.74581884e-01,   2.80254207e-01,   3.72906988e-01,&
       1.93337639e-02,   4.46672298e-01,   4.20528658e-01,&
       6.03508045e-01,   1.98176820e-01/)
  kext=-99.
  salb=-99.
  asym=-99.
  
  
  beta=0.75
  q=0.2*log(10.)*beta
  lwc=0
  z35mod=-99.
  dn=1.
  iNoAd=0
!  print*, node
!  write(*,*) z13obs(node(1):node(5))
!  stop
!  call interpolPC(imembC+1, rhPCij, cldwPCij, cldw, rh)
  !if(imembc+1==1) stop
  it=0
10 continue
  z13=z13obs
  !print*, node
  j=0
  rms=1e5
  !if(itype==1 .and. imembc==1) then
  !   print*, z13obs(node(2):node(4))
  !   print*, z13obs(node(2)), z13obs(node(3)), z13obs(node(4))
!     print*, 'max=',
  if(node(4)<88 .and. node(3)<87) then
     if(maxval(z13obs(node(2):node(4)))== z13obs(node(3))) then
        if(itype==1) then
           itype=11
        endif
     else
        if(maxval(z13obs(node(2):node(4)))== z13obs(node(3)-1) .and. &
             itype==1) itype=11
        if(maxval(z13obs(node(2):node(4)))== z13obs(node(3)+1) .and. &
             itype==1) itype=11
     endif
     if(itype==11) then
        dz1=maxval(z13obs(node(2):node(4)))-z13obs(node(2))
        dz2=maxval(z13obs(node(2):node(4)))-z13obs(node(4))
     endif
  endif
  !if(itype==1) itype=11
  if(node(1)<node(3)) then
     do i=max(node(1)-5,1),node(1)
        z13(i)=z13obs(node(1))-(2.*rv(mod(imembC,50)+1)+1)*(node(1)-i)  !2*ran1()
        log10dNw(i)=log10dNw(node(1))
        call integrateanvHB(z13,z35,i,i,pia13,&     ! this subroutine returns
             pia35,z35mod,lwc(:),log10dNw(:), &    ! precipitation parameters 
             ngates,nmfreqm,node,dr,imu(:),&
             dpia13,dpia35,&  !and electromagnetic properties                  
             kext(:,:),salb(:,:),asym(:,:),rrate,d0)
        z13(i)=z13obs(i)
     enddo
  endif

  !endif
  if(abs(log10dnw(node(1))+0.4589295)<1e-3) then
  !print*, log10dnw(node(1)),abs(log10dnw(node(1))+0.4589295)
     if(ipix==26 .and. jpix==4409) then
        !print*, log10dnw(node(1):node(5))
     endif
  endif
  do while(rms>1e-5 .and. j <40)
     j=j+1
     pia13=0
     pia13v=0
     pia35=0
     rms=0
     do i=node(1),node(5)
        if(z13(node(1))>60) then
        endif
!  SFM  start  06/22/2014; for M.Grecu (unknown justification)
        if(z13obs(i)>10) then
           if(isnan(log10dNw(i))) then
              log10dNw(i)=-0.5
              iNoAd=1
           endif
!  SFM  start  06/22/2014
           if(itype==11) then
              call integratestHB(z13,z13obs,z35,i,i,pia13,&                ! this subroutine returns
                   pia35,z35mod,lwc(:),log10dNw(:), &               ! precipitation parameters 
                   ngates,nmfreqm,node,dr,imu(:),dpia13,dpia35,&    ! and electromagnetic properties
                   kext(:,:),salb(:,:),asym(:,:),rrate,d0,dz1,dz2)          ! as a function of z13(i)
           else                                                     ! and log10dNw(i) 
              call integratecvHB(z13,z35,i,i,pia13,&
                   pia35,z35mod,lwc(:),log10dNw(:), &
                   ngates,nmfreqm,node,dr,imu(:),dpia13,dpia35,&
                   kext(:,:),salb(:,:),asym(:,:),rrate,d0)
           endif
           pia13=pia13+dpia13*dr
           pia35=pia35+dpia35*dr
           if(isnan(d0(i)) .or. isnan(rrate(i))) then
              log10dNw(i)=-0.5
              kext(i,:)=0
              salb(i,:)=0
              asym(i,:)=0
              rrate(i)=-99
              d0(i)=0.
              dpia13=0
              dpia35=0
           endif
           rms=rms+ (z13obs(i)+pia13-z13(i))**2
           z13(i)=z13obs(i)+pia13
!  SFM  start  06/22/2014; for M.Grecu (unknown justification)
           if(z13(i)>60) z13(i)=60.
!  SFM  start  06/22/2014
           
           if(itype==11) then
              call integratestHB(z13,z13obs,z35,i,i,pia13,&                ! this subroutine returns
                   pia35,z35mod,lwc(:),log10dNw(:), &               ! precipitation parameters 
                   ngates,nmfreqm,node,dr,imu(:),dpia13,dpia35,&    ! and electromagnetic properties
                   kext(:,:),salb(:,:),asym(:,:),rrate,d0,dz1,dz2)          ! as a function of z13(i)
           else                                                     ! and log10dNw(i) 
              call integratecvHB(z13,z35,i,i,pia13,&
                   pia35,z35mod,lwc(:),log10dNw(:), &
                  ngates,nmfreqm,node,dr,imu(:),dpia13,dpia35,&
                   kext(:,:),salb(:,:),asym(:,:),rrate,d0)
           endif
           if(isnan(d0(i)) .or. isnan(rrate(i))) then
              print*, 'Exception in fhb1.f90'
              print*, z13(i), z13obs(i), iEx, ipix,jpix
              print*, log10dNw(i), pia13,i
              iEx=iEx+1
              log10dNw(i)=-0.5
              kext(i,:)=0
              salb(i,:)=0
              asym(i,:)=0
              rrate(i)=-99
              d0(i)=0.
              z13(i)=z13obs(i)
              pia13=0
              pia35=0
              dpia13=0
              dpia35=0
              log10dNw(node(3):node(5))=log10dNw(node(3):node(5))-0.1
              if(iEx==3) then
                 kext(i,:)=0
                 salb(i,:)=0
                 print*, pia13srt, reliabFlag
                 print*, itype
                 print*, node
                 print*, z13obs(node(1):node(5))
                 print*, log10dnw(node(1):node(5))
                 stop
              endif
              goto 10
           endif
           pia13=pia13+dpia13*dr
           pia35=pia35+dpia35*dr
           pia13v=pia13v+dpia13*2*dr
 
        else
           dpia13=0.
           dpia35=0.
           rrate(i)=0.
           lwc(i)=0.
           d0(i)=0
        endif
        !pia35=pia35+wv_extMemb(i0)*dr
        !pia35=pia35+cld_extMemb(i0)*dr

        !pia35=pia35+atm_extKa(i0,ic)*dr
        !pia35=pia35+cld_extKa(i0,jc)*dr

       
     enddo
   
     pia35=pia35+dpia35*(isurf-node(5))*2.*dr
     pia13=pia13+dpia13*(isurf-node(5))*2.*dr
     pia13v=pia13v+dpia13*(isurf-node(5))*2.*dr
     
     do i=node(1),node(5)
        if(rrate(i)>250 .and. iNoAd==0) then
           log10dNw(i)=log10dNw(i)+log10(250./rrate(i))
        endif
      
     enddo
     !write(*,*) pia13, 'ITER=',j
  enddo
  if(imembC/3*3==imembC .and.  (reliabFlag==1 .or. reliabFlag==2) &
      .and. itype==2) then
     if(pia13srt>6 .and. pia13<0.75*pia13srt .and. &
          it<3 .and. pia13>0.1 .and. iNoAd==0) then
        !do i=node(2),88
        !print*, node(2)
        do i=node(2),88
           log10dNw(i)=log10dNw(i)-log(pia13/pia13srt)
           if(log10dNw(i)<-3) log10dNw(i)=-3
           if(log10dNw(i)>3) log10dNw(i)=3
        enddo
        it=it+1
        goto 10
     endif
  endif

  do i=node(1),node(5)
     if(isnan(log10dNw(i))) then
        !print*, 'nan in fhb1',log10dNw(i), z13(i), it, pia13srt
        !print*,pia13
        !stop
     endif
     
  enddo

  

!  if(RMS>1e-5 ) write(*,*) 'RMS', rms, j
!  if(pia13>3) write(*,*) 'PIAS=',pia13, pia13v, node(5), isurf
!  write(*,*) isurf, node(5), pia13

!  if(rrate(node(5))<-9) then
!     write(*,*) rrate(node(5)), itype
!     stop
!  endif
  if(itype==11) itype=1
end subroutine fhb11


