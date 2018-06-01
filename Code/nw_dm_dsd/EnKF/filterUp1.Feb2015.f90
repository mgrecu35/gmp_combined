module writeEnkF
integer :: iwENKF
end module writeEnkF

!  SFM  begin  07/29/2014; for M.Grecu  eliminate NANs
subroutine filterUpNS(dPRData,dPRRet, Xens,Yens,Yobs,Xup,tb,&
     dprRain,sfcRain,nmemb1,ic,i,j,&
     nxu,nyu,wfractm,s0KuVar,s0KaVar,s0Cov,hFreqTb)
!  SFM  end    07/29/2014
  use f90DataTypes
  use f90Types
  use writeEnKF
  implicit none
  type (dPRDataType)     :: dPRData
  type (dPRRetType)     :: dPRRet
  real :: Xens(nxu,nmemb1),Yens(nyu,nmemb1),Yobs(nyu), Xup(nxu)
  integer :: nmemb1,ic
  integer :: nxu,nyu
  integer :: ink,nx, ny, i, j, k
  integer :: ipias(2)
  integer :: ipol(15), ifreq(15), iobs(15)
  real :: tb(49,300,15), dprrain(49,300)!, tbRgrid(9,49,9300), 
  real :: sfcRain(49,300),sfcRainStd(49,300)
  real :: s0KuVar(49,300), s0KaVar(49,300), s0Cov(49,300), hFreqTb(4)
!  SFM  begin  07/29/2014; for M.Grecu  eliminate NANs
  real :: ymean, ystd, wfractm
  integer :: imemb, ifreqS2(4), ipolS2(4), iprint
  real ::  stddev
!  SFM  end    07/29/2014

!!  integer*4 :: ob_flag(100)
!!  ob_flag = -9

  ink=0 

  Xens(1,1:1*nMemb1)=(dPRRet%sfcRainEns(i,j,1:1*nmemb1))
  nx=1 
  do k=dPRData%node(1,i,j)+1,dPRData%node(5,i,j)+1
!  SFM  begin  07/29/2014; for M.Grecu  eliminate NANs
!  SFM  begin  06/22/2014; for M.Grecu (unknonw justification)
!begin  WSO 9/11/4 change zka1c21 to zku1c21
     if( dPRData%zku1c21(k,i,j)>12 .and. dPRRet%rrate(1,k,i,j)>0) then
!end    WSO 9/11/14
!  SFM  end    06/22/2014
!  SFM  end    07/29/2014
        nx=nx+1
        Xens(nx,1:nMemb1)=dPRRet%rrate(1:nmemb1,k,i,j)
        nx=nx+1
        Xens(nx,1:nMemb1)=dPRRet%pwc(1:nmemb1,k,i,j)
        nx=nx+1
        Xens(nx,1:nMemb1)=dPRRet%d0(1:nmemb1,k,i,j)
        nx=nx+1
        Xens(nx,1:nMemb1)=dPRRet%log10dNw(1:nmemb1,k,i,j)
        nx=nx+1
        Xens(nx,1:nMemb1)=dPRRet%z13c(1:nmemb1,k,i,j)
     endif
  enddo
  !Start SJM  12/4/2014
  nx=nx+1
  Xens(nx,1:nMemb1)=dPRRet%pia13mod(i,j,1:nmemb1)
  if(wfractm .gt. 99) then
    nx=nx+1
    Xens(nx,1:nMemb1)=dPRRet%simSigmaZeroKu(i,j,1:nmemb1)
    nx=nx+1
    Xens(nx,1:nMemb1)=dPRRet%sfcWindEns(i,j,1:nmemb1)
  endif
  !end SJM 12/4/2014
  do k=1,10
     nx=nx+1
     Xens(nx,1:nmemb1)=dPRRet%cldwcoeff(i,j,k,1:nmemb1)
  enddo
  !print*, 'nx1=',nx
  do k=1,8
     nx=nx+1
     Xens(nx,1:nmemb1)=dPRRet%tb(i,j,1,k,1:1*nmemb1)
     nx=nx+1
     Xens(nx,1:nmemb1)=dPRRet%tb(i,j,2,k,1:1*nmemb1)
  enddo
120 format(8(F6.2,1x))
  
  ny=ink+2
  Yobs(1:2)=0
  Yens(1:2,:)=0.

  if(dPRData%snrRatioku(i, j) > 2.) then! .and. dPRData%srtrelPIAku(i,j) > 2.9) then  !SJM 12/8/2014 should not need to filter by relPIA if wind variance is accounted for
     if(wfractm .le. 99) then
       Yobs(ink+1) = dPRData%srtPIAku(i,j)
       Yens(ink+1,1:1*nMemb1)=1*dPRRet%pia13mod(i,j,1:nmemb1)
     else
       !Start SJM  12/4/2014
       !Yobs(ink+1) = dPRData%sigmaZeroKu(i,j)
       !Yens(ink+1,1:1*nMemb1)=dPRRet%simSigmaZeroKu(i,j,1:nmemb1)
       Yobs(ink+1) = dPRData%srtPIAku(i,j)
       Yens(ink+1,1:1*nMemb1)=1*dPRRet%pia13mod(i,j,1:nmemb1)
       
     endif
     
   else
     Yobs(ink+1) = 0
     Yens(ink+1,1:1*nMemb1)=0
     !end SJM 12/4/2014
  endif
     Yobs(ink+2)=-9999.
     Yens(ink+2,:)=-9999.

!end    WSO 9/10/13

  ipias=(/ink+1,ink+2/)
  ink=ink+2

  ifreq(1:9)=(/1,1,2,2,3,4,4,5,5/)
!begin  MG 9/18/13
  ipol(1:9)=(/1,2,1,2,1,1,2,1,2/)
!end    MG 9/18/13
  iobs(1:9)=(/1,2,3,4,5,6,7,8,9/)  !SJm 12/8/2014 - fixed 6->7 in 5th element (should be 37H not 37V)

  do k=1,9!noTb
     if( dPRRet%tb(i,j,ipol(k),ifreq(k),nmemb1)>50 .and. &
          tb(i,j,iobs(k))>0  ) then
        ymean=sum(dPRRet%tb(i,j,ipol(k),ifreq(k),1:1*nmemb1))/nmemb1
        ystd=sqrt(sum((dPRRet%tb(i,j,ipol(k),ifreq(k),1:1*nmemb1)-&
             ymean)**2)/nmemb1)
        if(wfractm>80) then
           if(abs(tb(i,j,iobs(k))-ymean) .lt. 3.5*ystd) then
              ink=ink+1
              Yens(ink,1:1*nMemb1) = dPRRet%tb(i,j,ipol(k),ifreq(k),&
                   1:1*nmemb1)
              Yobs(ink)=tb(i,j,iobs(k))
           endif
        else
           if(abs(tb(i,j,iobs(k))-ymean) .lt. 3.5*ystd) then
              ink=ink+1
              Yens(ink,1:1*nMemb1) = dPRRet%tb(i,j,ipol(k),ifreq(k),&
                   1:1*nmemb1)
              Yobs(ink)=tb(i,j,iobs(k))
           endif
        endif
     endif
  enddo
  ipolS2=(/1,2,1,1/)
  ifreqS2=(/6,6,7,8/)
  do k=1,4!noTb
     if(hFreqTb(k)>0 .and. dPRRet%tb(i,j,ipolS2(k),ifreqS2(k),1)>0) then
        ink=ink+1
        Yobs(ink)=hFreqTb(k)
        Yens(ink,1:nMemb1)=dPRRet%tb(i,j,ipolS2(k),ifreqS2(k),1:1*nmemb1)
     endif
  enddo
  ny=ink
  !tbout(10)=sum(dPRRet%tb(i,j,1,6,1:1*nmemb1))/(nmemb1)
  !tbout(11)=sum(dPRRet%tb(i,j,2,6,1:1*nmemb1))/(nmemb1)
  !tbout(12)=sum(dPRRet%tb(i,j,1,7,1:1*nmemb1))/(nmemb1)
  !tbout(13)=sum(dPRRet%tb(i,j,1,8,1:1*nmemb1))/(nmemb1)
  !do k=10,13
  !   tbout(k)=hFreqPRg(i,j,k-9)
  !enddo
  if(ink>=2) then
     call enkF1d(Xens(1:nx,1:nmemb1-1),Yens(1:ink,1:nmemb1-1),&
          Yobs(1:ink),nx,ny,nMemb1-1,xup,ipias,s0KuVar(i,j),s0KaVar(i,j),s0Cov(i,j))
     sfcRain(i,j)=xup(1)
     dPRRet%sfcRainEns(i,j,1:1*nmemb1)=Xens(1,1:1*nMemb1)
  endif
  
 
  nx=0  
  nx=nx+1
  dPRRet%sfcRainEns(i,j,1:1*nmemb1)=Xens(1,1:1*nMemb1)

!begin  MG 9/18/13 remove spurious minus sign
  do k=dPRData%node(1,i,j)+1,dPRData%node(5,i,j)+1
!end    MG 9/18/13
!  SFM  begin  07/29/2014; for M.Grecu eliminate NANs
!  SFM  begin  06/22/2014; for M.Grecu (unknown justification)
!begin  WSO 9/11/4 change zka1c21 to zku1c21
     if( dPRData%zku1c21(k,i,j)>12 .and. dPRRet%rrate(1,k,i,j)>0) then
!end    WSO 9/11/14
!  SFM  end    06/22/2014
!  SFM  end    07/29/2014
        nx=nx+1
        dPRRet%rrate(1:nmemb1,k,i,j)=Xens(nx,1:nMemb1)
        nx=nx+1
        dPRRet%pwc(1:nmemb1,k,i,j)=Xens(nx,1:nMemb1)
        nx=nx+1
        dPRRet%d0(1:nmemb1,k,i,j)=Xens(nx,1:nMemb1)
        nx=nx+1
        dPRRet%log10dNw(1:nmemb1,k,i,j)=Xens(nx,1:nMemb1)
        nx=nx+1
        dPRRet%z13c(1:nmemb1,k,i,j)=Xens(nx,1:nMemb1)
     endif
  enddo
  !Start SJM  12/4/2014
  nx=nx+1
  dPRRet%pia13mod(i,j,1:nmemb1)=Xens(nx,1:nMemb1)
  if(wfractm .gt. 99) then  
     !print*, dPRData%srtPIAku(i,j), sum(dPRRet%pia13mod(i,j,1:nmemb1))/nmemb1
  endif
  if(wfractm .gt. 99) then
    nx=nx+1
    dPRRet%simSigmaZeroKu(i,j,1:nmemb1)=Xens(nx,1:nMemb1)
    nx=nx+1
    !print*, stddev(dPRRet%sfcWindEns(i,j,1:nmemb1),nmemb1)
    !print*, i, j,  sum(dPRRet%sfcWindEns(i,j,1:nmemb1)-Xens(nx,1:nMemb1))/nmemb1
    dPRRet%sfcWindEns(i,j,1:nmemb1)=Xens(nx,1:nMemb1)
    do k=1,nmemb1
      if(dPRRet%sfcWindEns(i,j,k) .lt. 0.) dPRRet%sfcWindEns(i,j,k) = 0.
      if(dPRRet%sfcWindEns(i,j,k) .gt. 50.) dPRRet%sfcWindEns(i,j,k) = 50.
    end do
  endif
  !end SJM 12/4/2014
  !nx=nx+1
  !dPRRet%pia13mod(i,j,1:nmemb1)=Xens(nx,1:nMemb1)
  do k=1,10
     nx=nx+1
     dPRRet%cldwcoeff(i,j,k,1:nmemb1)=Xens(nx,1:nmemb1)
  enddo
  !print*, 'nx2=',nx
  do k=1,8
     nx=nx+1
     !print*,  i,j,k,1,&
     !     sum(dPRRet%tb(i,j,1,k,1:1*nmemb1)-Xens(nx,1:nmemb1))/nmemb1
     if(k==7) then
        !print*,  sum(dPRRet%tb(i,j,1,k,1:1*nmemb1))/nmemb1, &
        !     sum(Xens(nx,1:nmemb1))/nmemb1, hFreqTb(3)
     endif
     dPRRet%tb(i,j,1,k,1:1*nmemb1)= Xens(nx,1:nmemb1)

     nx=nx+1
     dPRRet%tb(i,j,2,k,1:1*nmemb1)= Xens(nx,1:nmemb1)
  enddo
50 format(20(F8.3,1x))

end subroutine filterUpNS

subroutine filterUpNSLand(dPRData,dPRRet, Xens,Yens,Yobs,Xup,tb,dprRain,sfcRain,nmemb1,ic,i,j,&
     nxu,nyu,wfractm,s0KuVar,s0KaVar,s0Cov, hfreqTb)
!  SFM  end    07/29/2014
  use f90DataTypes
  use f90Types
  implicit none
  type (dPRDataType)     :: dPRData
  type (dPRRetType)     :: dPRRet
  real :: Xens(nxu,nmemb1),Yens(nyu,nmemb1),Yobs(nyu), Xup(nxu)
  integer :: nmemb1,ic
  integer :: nxu,nyu
  integer :: ink,nx, ny, i, j, k
  integer :: ipias(2)
  integer :: ipol(15), ifreq(15), iobs(15)
  real :: tb(49,300,15), dprrain(49,300)!, tbRgrid(9,49,9300), 
  real :: sfcRain(49,300),sfcRainStd(49,300)
  real :: s0KuVar(49,300), s0KaVar(49,300), s0Cov(49,300), hFreqTb(4)
!  SFM  begin  07/29/2014; for M.Grecu  eliminate NANs
  real :: ymean, ystd, wfractm
  integer :: imemb, ifreqS2(4), ipolS2(4)
  real ::  stddev
!  SFM  end    07/29/2014

!!  integer*4 :: ob_flag(100)
!!  ob_flag = -9

  ink=0 

  Xens(1,1:1*nMemb1)=(dPRRet%sfcRainEns(i,j,1:1*nmemb1))
  nx=1 
  do k=dPRData%node(1,i,j)+1,dPRData%node(5,i,j)+1
!  SFM  begin  07/29/2014; for M.Grecu  eliminate NANs
!  SFM  begin  06/22/2014; for M.Grecu (unknonw justification)
!begin  WSO 9/11/4 change zka1c21 to zku1c21
     if( dPRData%zku1c21(k,i,j)>12 .and. dPRRet%rrate(1,k,i,j)>0) then
!end    WSO 9/11/14
!  SFM  end    06/22/2014
!  SFM  end    07/29/2014
        nx=nx+1
        Xens(nx,1:nMemb1)=dPRRet%rrate(1:nmemb1,k,i,j)
        nx=nx+1
        Xens(nx,1:nMemb1)=dPRRet%pwc(1:nmemb1,k,i,j)
        nx=nx+1
        Xens(nx,1:nMemb1)=dPRRet%d0(1:nmemb1,k,i,j)
        nx=nx+1
        Xens(nx,1:nMemb1)=dPRRet%log10dNw(1:nmemb1,k,i,j)
        nx=nx+1
        Xens(nx,1:nMemb1)=dPRRet%z13c(1:nmemb1,k,i,j)
     endif
  enddo
  !Start SJM  12/4/2014
  nx=nx+1
  Xens(nx,1:nMemb1)=dPRRet%pia13mod(i,j,1:nmemb1)

  do k=1,10
     nx=nx+1
     Xens(nx,1:nmemb1)=dPRRet%cldwcoeff(i,j,k,1:nmemb1)
  enddo
  !print*, 'nx1=',nx
  do k=1,8
     nx=nx+1
     Xens(nx,1:nmemb1)=dPRRet%tb(i,j,1,k,1:1*nmemb1)
     nx=nx+1
     Xens(nx,1:nmemb1)=dPRRet%tb(i,j,2,k,1:1*nmemb1)
  enddo
120 format(8(F6.2,1x))
  
  ny=ink+2
  Yobs(1:2)=0
  Yens(1:2,:)=0.

  if(dPRData%snrRatioku(i, j) > 2.) then! .and. dPRData%srtrelPIAku(i,j) > 2.9) then  !SJM 12/8/2014 should not need to filter by relPIA if wind variance is accounted for
     if(dPRData%dsrtrelPIA(i,j)>4) then
        Yobs(ink+1) = dPRData%dsrtPIAKa(i,j)-dPRData%dsrtPIAKu(i,j)
        Yens(ink+1,1:1*nMemb1)=dPRRet%pia35mod(i,j,1:nmemb1)-&
             dPRRet%pia13mod(i,j,1:nmemb1)
     endif
     if(dPRData%srtrelPIAku(i,j) > 4) then
        Yobs(ink+1) = dPRData%srtPIAKu(i,j)
        Yens(ink+1,1:1*nMemb1)=dPRRet%pia13mod(i,j,1:nmemb1)
     endif
     
  else
     Yobs(ink+1) = 0
     Yens(ink+1,1:1*nMemb1)=0
     !end SJM 12/4/2014
  endif
  Yobs(ink+2)=-9999.
  Yens(ink+2,:)=-9999.
  
!end    WSO 9/10/13

  ipias=(/ink+1,ink+2/)
  ink=ink+2

  ifreq(1:9)=(/1,1,2,2,3,4,4,5,5/)
!begin  MG 9/18/13
  ipol(1:9)=(/1,2,1,2,1,1,2,1,2/)
!end    MG 9/18/13
  iobs(1:9)=(/1,2,3,4,5,6,7,8,9/)  !SJm 12/8/2014 - fixed 6->7 in 5th element (should be 37H not 37V)

!!write(unit=98,fmt=*) ' =================== '
!!write(unit=98,fmt=*) ' index ',j+ic
!!write(unit=98,fmt=2211) tbRgrid(:,:,j+ic)
!!2211 format(10f11.4)
!!write(unit=98,fmt=*) ' =================== '
 ! print*, i, j, ink
  do k=1,9 !noTb
     if( dPRRet%tb(i,j,ipol(k),ifreq(k),nmemb1)>50 .and. &
          tb(i,j,iobs(k))>0  ) then
!  SFM  begin  07/29/2014; for M.Grecu eliminate NANs
        ymean=sum(dPRRet%tb(i,j,ipol(k),ifreq(k),1:1*nmemb1))/nmemb1
        ystd=sqrt(sum((dPRRet%tb(i,j,ipol(k),ifreq(k),1:1*nmemb1)-&
             ymean)**2)/nmemb1)

        if(wfractm>80) then
           if(abs(tb(i,j,iobs(k))-ymean) .lt. 3.5*ystd) then
              ink=ink+1
              Yens(ink,1:1*nMemb1) = dPRRet%tb(i,j,ipol(k),ifreq(k),&
                   1:1*nmemb1)
              Yobs(ink)=tb(i,j,iobs(k))
           endif
        else
           if(abs(tb(i,j,iobs(k))-ymean) .lt. 3.5*ystd) then
              ink=ink+1
              Yens(ink,1:1*nMemb1) = dPRRet%tb(i,j,ipol(k),ifreq(k),&
                   1:1*nmemb1)
              Yobs(ink)=tb(i,j,iobs(k))
           endif
        endif
!  SFM  end  07/29/2014
     endif
  enddo

  ipolS2=(/1,2,1,1/)
  ifreqS2=(/6,6,7,8/)
  
  do k=1,4!noTb
     if(hFreqTb(k)>0 .and. dPRRet%tb(i,j,ipolS2(k),ifreqS2(k),1)>0 ) then
        ink=ink+1
        Yobs(ink)=hFreqTb(k)
        Yens(ink,1:nMemb1)=dPRRet%tb(i,j,ipolS2(k),ifreqS2(k),1:1*nmemb1)
     endif
  enddo
  ny=ink
  !print*, 'tbInc ',i, j, ink, ny
  !print*, ny
 
  !print*, Xens(1,1), ink
  if(ink>=2) then
     call enkF1d(Xens(1:nx,1:nmemb1-1),Yens(1:ink,1:nmemb1-1),&
          Yobs(1:ink),nx,ny,nMemb1-1,xup,ipias,25.,25.,0.)
     sfcRain(i,j)=xup(1)
     dPRRet%sfcRainEns(i,j,1:1*nmemb1)=Xens(1,1:1*nMemb1)
  endif
  
 
  nx=0  
  nx=nx+1
  dPRRet%sfcRainEns(i,j,1:1*nmemb1)=Xens(1,1:1*nMemb1)

!begin  MG 9/18/13 remove spurious minus sign
  do k=dPRData%node(1,i,j)+1,dPRData%node(5,i,j)+1
!end    MG 9/18/13
!  SFM  begin  07/29/2014; for M.Grecu eliminate NANs
!  SFM  begin  06/22/2014; for M.Grecu (unknown justification)
!begin  WSO 9/11/4 change zka1c21 to zku1c21
     if( dPRData%zku1c21(k,i,j)>12 .and. dPRRet%rrate(1,k,i,j)>0) then
        !end    WSO 9/11/14
        !  SFM  end    06/22/2014
        !  SFM  end    07/29/2014
        nx=nx+1
        dPRRet%rrate(1:nmemb1,k,i,j)=Xens(nx,1:nMemb1)
        nx=nx+1
        dPRRet%pwc(1:nmemb1,k,i,j)=Xens(nx,1:nMemb1)
        nx=nx+1
        dPRRet%d0(1:nmemb1,k,i,j)=Xens(nx,1:nMemb1)
        nx=nx+1
        dPRRet%log10dNw(1:nmemb1,k,i,j)=Xens(nx,1:nMemb1)
        nx=nx+1
        dPRRet%z13c(1:nmemb1,k,i,j)=Xens(nx,1:nMemb1)
     endif
  enddo
  !Start SJM  12/4/2014
  nx=nx+1
  dPRRet%pia13mod(i,j,1:nmemb1)=Xens(nx,1:nMemb1)
 
  do k=1,10
     nx=nx+1
     dPRRet%cldwcoeff(i,j,k,1:nmemb1)=Xens(nx,1:nmemb1)
  enddo
  !print*, 'nx2=',nx
  do k=1,8
     nx=nx+1
     dPRRet%tb(i,j,1,k,1:1*nmemb1)= Xens(nx,1:nmemb1)
     nx=nx+1
     dPRRet%tb(i,j,2,k,1:1*nmemb1)= Xens(nx,1:nmemb1)
  enddo
  
50 format(20(F8.3,1x))

end subroutine filterUpNSLand

!    Medium Swath filtering

!  SFM  begin  07/29/2014; for M.Grecu  eliminate NANs
subroutine filterUpMS(dPRData,dPRRet, Xens,Yens,Yobs,Xup,tb,&
     dprRain,sfcRain,nmemb1,ic,i,j,&
     nxu,nyu,wfractm,s0KuVar,s0KaVar,s0Cov, hFreqTb)
!  SFM  end    07/29/2014
  use f90DataTypes
  use f90Types
  use writeENKF
  implicit none
  type (dPRDataType)     :: dPRData
  type (dPRRetType)     :: dPRRet
  real :: Xens(nxu,nmemb1),Yens(nyu,nmemb1),Yobs(nyu), Xup(nxu)
  integer :: nmemb1,ic
  integer :: nxu,nyu
  integer :: ink,nx, ny, i, j, k, iprint
  integer :: ipias(2)
  integer :: ipol(15), ifreq(15), iobs(15)
  real :: tb(49,300,15), dprrain(49,300)
  real :: sfcRain(49,300),sfcRainStd(49,300)
  real :: s0KuVar(49,300), s0KaVar(49,300), s0Cov(49,300), hFreqTb(4)
!  SFM  begin  07/29/2014; for M.Grecu  eliminate NANs
  real :: ymean, ystd, wfractm
  integer :: imemb, ifreqS2(4), ipolS2(4)
  integer :: nZKa, iZKa(88),k1
!  SFM  end    07/29/2014

!!  integer*4 :: ob_flag(100)
!!  ob_flag = -9

  ink=0 
  if(iwENKF==1) then
     call enkfwi1(nmemb1)
     call enkfwi1(dPRData%rainType(i,j))
     call enkfwi5(dPRData%node(:,i,j))
     call enkfwi1(i)
     call enkfwi1(j)
  endif
  nZKa=0
  do k=dPRData%node(1,i,j)+1,dPRData%node(5,i,j)+1
!  SFM  begin  06/22/2014; for M.Grecu (unknown justification)
     !print*, dPRData%zka1c21(k,i,j), dPRRet%z35mod0(1,k,i,j)
     if( dPRData%zka1c21(k,i,j)>15 .and. minval( dPRRet%z35mod0(1:1*nmemb1,k,i,j))>0) then
!  SFM  end    06/22/2014
        ink=ink+1
        Yens(ink,1:1*nMemb1) = dPRRet%z35mod0(1:1*nmemb1,k,i,j)
        Yobs(ink)=dPRData%zka1c21(k,i,j)
        iZKa(ink)=k
!!        ob_flag(ink) = 5
     endif
  enddo
  nZKa=ink
  if(iwENKF==1) then
     call enkfwi1(nZka)
     call enkfwf(wfractm)
     call enkfwf(dPRData%freezH(i,j) /1000.)
     do k=dPRData%node(1,i,j)+1,dPRData%node(5,i,j)+1
        call enkfwf(dPRData%zku1c21(k,i,j))
     enddo
     do k=1,nZka
        call enkfwi1(iZka(k))
        do k1=1,nmemb1
           call enkfwf(Yens(k,k1))
        enddo
        call enkfwf(Yobs(k))
     enddo
     do k=dPRData%node(1,i,j)+1,dPRData%node(5,i,j)+1
        do k1=1,nmemb1
           call enkfwf(dPRRet%MS%zKuEns(k1,k,i,j))
        enddo
        do k1=1,nmemb1
           call enkfwf(dPRRet%MS%log10dnw(k1,k,i,j))
        enddo
        do k1=1,nmemb1
           call enkfwf(dPRRet%MS%rrate(k1,k,i,j))
        enddo
     enddo
     do k1=1,nmemb1
        call enkfwf(dPRRet%MS%sfcRainEns(i,j,k1))
     enddo
     call enkfwf(dPRData%dsrtPIAku(i, j))
     call enkfwf(dPRData%dsrtPIAka(i, j))
     do k1=1,nmemb1
        call enkfwf(dPRRet%MS%pia13mod(i, j, k1))
        call enkfwf(dPRRet%MS%pia35mod(i, j, k1))
     enddo
     call enkfwf(dPRData%srtPIAku(i, j))
     call enkfwf(dPRData%srtrelPIAku(i, j))
     call enkfwf(dPRData%dsrtrelPIA(i, j))
  endif
     
  nx=0

  nx=nx+1
  Xens(1,1:1*nMemb1)=(dPRRet%MS%sfcRainEns(i,j,1:1*nmemb1))

  do k=dPRData%node(1,i,j)+1,dPRData%node(5,i,j)+1
!  SFM  begin  07/29/2014; for M.Grecu  eliminate NANs
!  SFM  begin  06/22/2014; for M.Grecu (unknown justification)
     if( dPRData%zku1c21(k,i,j)>12 .and. dPRRet%MS%rrate(1,k,i,j)>0) then
!  SFM  end    06/22/2014
!  SFM  end    07/29/2014
        nx=nx+1
        Xens(nx,1:nMemb1)=dPRRet%MS%rrate(1:nmemb1,k,i,j)
        nx=nx+1
        Xens(nx,1:nMemb1)=dPRRet%MS%pwc(1:nmemb1,k,i,j)
        nx=nx+1
        Xens(nx,1:nMemb1)=dPRRet%MS%d0(1:nmemb1,k,i,j)
        nx=nx+1
        Xens(nx,1:nMemb1)=dPRRet%MS%log10dNw(1:nmemb1,k,i,j)
        nx=nx+1
        Xens(nx,1:nMemb1)=dPRRet%MS%zKuEns(1:nmemb1,k,i,j)
        nx=nx+1
        Xens(nx,1:nMemb1)=dPRRet%MS%zKaEns(1:nmemb1,k,i,j)
        nx=nx+1
        Xens(nx,1:nMemb1)=dPRRet%z35mod0(1:nmemb1,k,i,j)
     endif
  enddo
  !Start SJM  12/9/2014
  nx=nx+1
  Xens(nx,1:nMemb1)=dPRRet%MS%pia13mod(i,j,1:nmemb1)
  nx=nx+1
  Xens(nx,1:nMemb1)=dPRRet%MS%pia35mod(i,j,1:nmemb1)
  if(wfractm .gt. 99) then
    nx=nx+1
    Xens(nx,1:nMemb1)=dPRRet%MS%simSigmaZeroKu(i,j,1:nmemb1)
    nx=nx+1
    Xens(nx,1:nMemb1)=dPRRet%MS%simSigmaZeroKa(i,j,1:nmemb1)
    nx=nx+1
    Xens(nx,1:nMemb1)=dPRRet%MS%sfcWindEns(i,j,1:nmemb1)
  endif
  !end SJM 12/9/2014
  !nx=nx+1
  !Xens(nx,1:nMemb1)=dPRRet%MS%pia13mod(i,j,1:nmemb1)
  !nx=nx+1
  !Xens(nx,1:nMemb1)=dPRRet%MS%pia35mod(i,j,1:nmemb1)

!begin  MG 9/18/13 addition
!  do k=1,10
!     nx=nx+1
!     Xens(nx,1:nmemb1)=dPRRet%MS%cldwcoeff(i,j,k,1:nmemb1)
!  enddo
  do k=1,8
     !print*, nx, ' in'
     nx=nx+1
     Xens(nx,1:nmemb1)=dPRRet%MS%tb(i,j,1,k,1:1*nmemb1)
     nx=nx+1
     Xens(nx,1:nmemb1)=dPRRet%MS%tb(i,j,2,k,1:1*nmemb1)
  enddo
!end    MG 9/18/13 addition
120 format(8(F6.2,1x))
  Yobs(ink+1:ink+2)=0
!!       ob_flag(ink+1:ink+2) = 1
  Yens(ink+1:ink+2,:)=0.
  ny=ink+2
!begin  WSO 9/5/13 uses DSRT for dual-wavelength applications
!  if(dPRData%srtrelPIAku(i,j)>2.9 .and. &
!       dPRData%zku1c21(1+dPRData%node(5,i,j),i,j)>0) then
!     Yobs(ink+1) = dPRData%srtPIAku(i,j)
!     Yens(ink+1,1:1*nMemb1)=1*dPRRet%pia13mod(i,j,1:nmemb1)
!     if(dPRData%dsrtPIAka(i,j)>0) then    
!        Yobs(ink+2)=dPRData%dsrtPIAka(i,j)
!        Yens(ink+2,:)=dPRRet%pia35mod(i,j,1:nmemb1)
!     else
!        Yobs(ink+2)=0.
!        Yens(ink+2,:)=0.
!     endif
!  endif
!end    WSO 9/5/13
!begin  WSO 9/10/13 uses DSRT for dual-wavelength applications
  if(dPRData%snrRatioka(i, j) > 2.) then  !substantial surface signal from Ka and Ku
    
   ! if(dPRData%dsrtrelPIAka(i,j) > 2.9) then
      if(wfractm .le. 99) then
        !delta-PIA is substantial relative
                                             !to standard deviation of non-raining cross section difference
        Yobs(ink+1) = dPRData%dsrtPIAka(i, j)-dPRData%dsrtPIAku(i, j)
!!      ob_flag(ink+1) = 2
        Yens(ink+1,1:1*nMemb1) =  dPRRet%MS%pia35mod(i, j, 1:nmemb1)-dPRRet%MS%pia13mod(i, j, 1:nmemb1)
!  SFM  begin  06/22/2014; for M.Grecu (unknown justification)
        Yobs(ink+2) = 0
!!      ob_flag(ink+2) = 2
!  SFM  begin  06/22/2014; for M.Grecu (unknown justification)
        Yens(ink+2,:) = 0
!  SFM  end    06/22/2014
        
        !Yobs(ink+1) = dPRData%dsrtrelPIA(i,j)
        !Yens(ink+1,1:1*nMemb1)=dPRRet%pia35mod(i,j,1:nmemb1)-&
        !     dPRRet%pia13mod(i,j,1:nmemb1)
        !Yobs(ink+2)=0
        !Yens(ink+2,:)=0.
      else
        !Start SJM  12/9/2014
        Yobs(ink+1) = dPRData%sigmaZeroKu(i,j)
        Yens(ink+1,1:1*nMemb1)=dPRRet%MS%simSigmaZeroKu(i,j,1:nmemb1)
        Yobs(ink+2) = dPRData%sigmaZeroKa(i,j)
        Yens(ink+2,1:1*nMemb1)=dPRRet%MS%simSigmaZeroKa(i,j,1:nmemb1)
        !end SJM 12/9/2014
      endif
      
   else if(dPRData%snrRatioku(i, j) > 2.) then   !substantial surface signal from Ku only
      if(dPRData%srtrelPIAku(i, j) > 2.9) then    !PIA at Ku is substantial
         if(wfractm .le. 99) then
                                                !relative to standard deviation of non-raining cross section

        Yobs(ink+1) = dPRData%srtPIAku(i, j)
!!      ob_flag(ink+1) = 2
        Yens(ink+1,1:1*nMemb1) = 1*dPRRet%MS%pia13mod(i, j, 1:nmemb1)
        Yobs(ink+2) = 0.
!!      ob_flag(ink+2) = 3
        Yens(ink+2,:) = 0.
      else
        !Start SJM  12/9/2014
        Yobs(ink+1) = dPRData%sigmaZeroKu(i,j)
        Yens(ink+1,1:1*nMemb1)=dPRRet%MS%simSigmaZeroKu(i,j,1:nmemb1)
        Yobs(ink+2) = 0.!dPRData%sigmaZeroKu(i,j)
        Yens(ink+2,1:1*nMemb1)=0.!dPRRet%simSigmaZeroKu(i,j,1:nmemb1)
        !end SJM 12/9/2014
      endif
    endif

  else !no substantial surface signal from Ku or Ka

!   punt!
    Yobs(ink+1) = 0.
!!      ob_flag(ink+1) = 1
    Yens(ink+1,1:1*nMemb1) = 0.
    Yobs(ink+2) = 0.
!!      ob_flag(ink+2) = 1
    Yens(ink+2,:) = 0.

  endif
  if(Yobs(ink+1)<0) then
     Yobs(ink+1)=0
     Yens(ink+1,1:1*nMemb1) = 0.
  endif
  if(Yobs(ink+2)<0) then
     Yobs(ink+2)=0
     Yens(ink+2,1:1*nMemb1) = 0.
  endif
!end    WSO 9/10/13
  Yobs(ink+1) = 0.
  !!      ob_flag(ink+1) = 1
  Yens(ink+1,1:1*nMemb1) = 0.
  Yobs(ink+2) = 0.
  !!      ob_flag(ink+2) = 1
  Yens(ink+2,:) = 0.

  ipias=(/ink+1,ink+2/)
  ink=ink+2

  ifreq(1:9)=(/1,1,2,2,3,4,4,5,5/)
!begin  MG 9/18/13
  ipol(1:9)=(/1,2,1,2,1,1,2,1,2/)
!end    MG 9/18/13
  iobs(1:9)=(/1,2,3,4,5,6,7,8,9/)  !SJm 12/8/2014 - fixed 6->7 in 5th element (should be 37H not 37V)
  
!!write(unit=99,fmt=*) ' =================== '
!!write(unit=99,fmt=*) ' index ',j+ic
!!write(unit=99,fmt=2211) tbRgrid(:,:,j+ic)
!!2211 format(10f11.4)
!!write(unit=99,fmt=*) ' =================== '
  if(iwEnkF==1) then
     do k=1,9
        call enkfwf(tb(i,j,k)) 
     enddo
     do k=1,4
        call enkfwf(hFreqTb(k))
     enddo
  endif

  do k=1,-9 !noTb
     if( dPRRet%MS%tb(i,j,ipol(k),ifreq(k),nmemb1)>0 .and. &
          tb(i,j,iobs(k))>0  ) then
        ymean=sum(dPRRet%MS%tb(i,j,ipol(k),ifreq(k),1:1*nmemb1))/nmemb1
        ystd=sqrt(sum((dPRRet%MS%tb(i,j,ipol(k),ifreq(k),1:1*nmemb1)-&
             ymean)**2)/nmemb1)
!  SFM  begin  07/29/2014; for M.Grecu  eliminate NANs
        if(wfractm>80) then
           if(abs(tb(i,j,iobs(k))-ymean) .lt. 5.5*ystd) then
              ink=ink+1
              Yens(ink,1:1*nMemb1) = dPRRet%MS%tb(i,j,ipol(k),ifreq(k),1:1*nmemb1)
              Yobs(ink)=tb(i,j,iobs(k))
           endif
        else
           if(abs(tb(i,j,iobs(k))-ymean) .lt. 3.5*ystd) then
              ink=ink+1
              Yens(ink,1:1*nMemb1) = dPRRet%MS%tb(i,j,&
                   ipol(k),ifreq(k),1:1*nmemb1)
              Yobs(ink)=tb(i,j,iobs(k))
           endif
!  SFM  end    07/29/2014
        endif
     endif
  enddo
  ipolS2=(/1,2,1,1/)
  ifreqS2=(/6,6,7,8/)
  do k=1,-4!noTb
     if(hFreqTb(k)>0 .and. dPRRet%tb(i,j,ipolS2(k),ifreqS2(k),1)>0) then
        ink=ink+1
        Yobs(ink)=hFreqTb(k)
        Yens(ink,1:nMemb1)=dPRRet%MS%tb(i,j,ipolS2(k),ifreqS2(k),1:1*nmemb1)
     endif
  enddo

  if(iwEnkF==1)then
     do k=1,9
        do k1=1,nmemb1
           call enkfwf(dPRRet%MS%tb(i,j,ipol(k),ifreq(k),k1))
        enddo
     enddo
     do k=1,4
        do k1=1,nmemb1
           call enkfwf(dPRRet%MS%tb(i,j,ipolS2(k),ifreqS2(k),k1))
        enddo
     enddo
  endif
  ny=ink


  if(Xens(1,1)>-5 .and. ink>=2) then
!begin  MG 9/18/13 changed arguments of Yobs
     !print*, nx, ny
     call enkF1d(Xens(1:nx,1:nmemb1-1),Yens(1:ink,1:nmemb1-1),&
          Yobs(1:ink),nx,ny,nMemb1-1,xup(1:nx),ipias,25.,s0KaVar(i,j),s0Cov(i,j))
!end    MG 9/18/13
     !sfcRain(i,j)=xup(1)
!begin  WSO 9/5/13 rename
     !if(dPRData%node(5,i,j)>dPRData%node(4,i,j)+2) &
          !write(*,50) (xup(1)), sfcRain(i,j), &!sum(dPRRet%sfcRainEns(i,j,1:1*nmemb1))/nmemb1, &
          !sum(dPRRet%MS%sfcRainEns(i,j,1:1*nmemb1))/nmemb1,dprrain(i,j),  &
          !sum(dPRRet%pia35mod(i,j,1:nmemb1))/nmemb1, dPRData%dsrtPIAka(i,j),&
          !sum(dPRRet%pia13mod(i,j,1:nmemb1))/nmemb1, dPRData%srtPIAku(i,j),& !,&
          !(sum(dPRRet%convtb(i,j,ipol(k),ifreq(k),1:1*nmemb1))/nmemb1,k=1,4), &
          !(tbRgrid(iobs(k),i,j+ic),k=1,4),  minval(Xens(1,:nmemb1)), maxval(Xens(1,:nmemb1))
!end    WSO 9/5/13

  endif
  nx=0
  nx=nx+1
  dPRRet%MS%sfcRainEns(i,j,1:1*nmemb1)=Xens(1,1:1*nMemb1)
  iprint=0
  do k=dPRData%node(1,i,j)+1,dPRData%node(5,i,j)+1
!  SFM  begin  07/29/2014; for M.Grecu  eliminate NANs
!  SFM  begin  06/22/2014; for M.Grecu (unknown justification)
     if( dPRData%zku1c21(k,i,j)>12 .and. dPRRet%MS%rrate(1,k,i,j)>0) then
!  SFM  end    06/22/2014
        nx=nx+1
        if(sum(Xens(nx,1:nMemb1))/nMemb1>300) then
          iprint=1
        endif
        dPRRet%MS%rrate(1:nmemb1,k,i,j)=Xens(nx,1:nMemb1)
         
        do imemb=1,nmemb1
           if(dPRRet%MS%rrate(imemb,k,i,j)<0) dPRRet%MS%rrate(imemb,k,i,j)=0
        enddo
        nx=nx+1
        dPRRet%MS%pwc(1:nmemb1,k,i,j)=Xens(nx,1:nMemb1)
        do imemb=1,nmemb1
           if(dPRRet%MS%pwc(imemb,k,i,j)<0) dPRRet%MS%pwc(imemb,k,i,j)=0
        enddo
!  SFM  end    07/29/2014
        nx=nx+1
        dPRRet%MS%d0(1:nmemb1,k,i,j)=Xens(nx,1:nMemb1)
        nx=nx+1
        dPRRet%MS%log10dNw(1:nmemb1,k,i,j)=Xens(nx,1:nMemb1)
        nx=nx+1
        dPRRet%MS%zKuEns(1:nmemb1,k,i,j)=Xens(nx,1:nMemb1)
        nx=nx+1
        dPRRet%MS%zKaEns(1:nmemb1,k,i,j)=Xens(nx,1:nMemb1)
        nx=nx+1
        dPRRet%z35mod0(1:nmemb1,k,i,j)=Xens(nx,1:nMemb1)
     endif
  enddo
  if(iprint==11) then
     do k=dPRData%node(1,i,j)+1,dPRData%node(5,i,j)+1
        if( dPRData%zku1c21(k,i,j)>12) then
           !print*, dPRData%zku1c21(k,i,j), dPRData%zka1c21(k,i,j)
        endif
        
     enddo
     print*, dPRData%dsrtPIAku(i, j), dPRData%dsrtPIAka(i, j)
     print*, Yobs(ipias)
     print*, dPRData%srtrelPIAku(i, j) 
     print*, dPRData%srtrelPIAku(i, j) 
     print*, dprdata%snrRatioku(i,j)
     print*, dprdata%snrRatioka(i,j)
  endif
  !Start SJM  12/4/2014
  nx=nx+1
  dPRRet%MS%pia13mod(i,j,1:nmemb1)=Xens(nx,1:nMemb1)
  nx=nx+1
  dPRRet%MS%pia35mod(i,j,1:nmemb1)=Xens(nx,1:nMemb1)
  if(wfractm .gt. 99) then
    nx=nx+1
    dPRRet%MS%simSigmaZeroKu(i,j,1:nmemb1)=Xens(nx,1:nMemb1)
    nx=nx+1
    dPRRet%MS%simSigmaZeroKa(i,j,1:nmemb1)=Xens(nx,1:nMemb1)
    nx=nx+1
    dPRRet%MS%sfcWindEns(i,j,1:nmemb1)=Xens(nx,1:nMemb1)
    do k=1,nmemb1
      if(dPRRet%MS%sfcWindEns(i,j,k) .lt. 0.) dPRRet%MS%sfcWindEns(i,j,k) = 0.
      if(dPRRet%MS%sfcWindEns(i,j,k) .gt. 50.) dPRRet%MS%sfcWindEns(i,j,k) = 50.
    end do
  endif
  !end SJM 12/4/2014
  !nx=nx+1
  !dPRRet%MS%pia13mod(i,j,1:nmemb1)=Xens(nx,1:nMemb1)
  !nx=nx+1
  !dPRRet%MS%pia35mod(i,j,1:nmemb1)=Xens(nx,1:nMemb1)
!begin  MG 9/18/13 addition
!  do k=1,10
!     nx=nx+1
!     dPRRet%MS%cldwcoeff(i,j,k,1:nmemb1)=Xens(nx,1:nmemb1)
!  enddo

  do k=1,8
     !print*, nx, ' out'
     nx=nx+1
     dPRRet%MS%tb(i,j,1,k,1:1*nmemb1)= Xens(nx,1:nmemb1)
     nx=nx+1
     dPRRet%MS%tb(i,j,2,k,1:1*nmemb1)= Xens(nx,1:nmemb1)
  enddo
!end MG 9/18/13 

!!50 format(20(F8.3,1x))
!!if (ny .ge. 1) write(unit=99,fmt=1234) 'MS YOBS ',ic+j, ny, Yobs(1:ny)
!!1234 format(a8,2i5,100f11.4)
!!if (ny .ge. 1) write(unit=99,fmt=1235) 'MS FLGS ',ob_flag(1:ny)
!1235 format(a8,10x,100i11)
end subroutine filterUpMS
