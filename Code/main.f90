program main
  use f90Types   ! allocates memory for the 1D 
  implicit none
  integer :: nmu,nmfreq
  type (radarRetType)    :: radarRet
  type (radarDataType)   :: radarData
  type (stormStructType) :: stormStruct
  integer :: ngates, nmemb1
  ngates=88
  nmemb1=5
  nmfreq=8
  nmu=5

  call readtablesLiang2(nmu,nmfreq) 
  call allocateDPRProfRet(radarRet,nmfreq,nmemb1,ngates, 9)
  call allocateDPRProfData(radarData, ngates)              
  call allocateStormStructData(stormStruct)
  ! allocates memory for the 5-node 
  ! allocates memory for the 
  radarData%ngates=ngates
  radarData%dr=0.25

  stormStruct%nodes  = nodes!dPRData%node(:,i,j)
  radarData%z13obs   = zku
  radarData%z35obs   = zka
  radarData%pia13srt = srtPIAku
  !radarData%relpia13srt =! dPRData%srtrelPIAku(i,j)
  !radarData%pia35srt = -99 
  !radarData%pia35srt = !dPRData%dsrtPIAka(i,j)
  radarData%hfreez   = freezH /1000. 
  stormStruct%iSurf = binRealSurface+1

  do k=1,nbin
     if(radarData%z13obs(k)<-99) radarData%z13obs(k)=-99
  enddo
  
  stormStruct%rainType=dPRData%rainType(i,j)
  stormStruct%rainType=stormStruct%rainType/100
  
  
  radarRet%logdNw(k*9+1:(k+1)*9)= & !Sept 17, 2015 MG
       logdNwf(k*9+1:(k+1)*9)
  call  ensRadRetStCvKu(radarData,stormStruct,                  &
       retParam, nmu,radarRet, itop, rms1, rms2, sysdN, iit, &
       xscalev, randemiss, dPRData%localZenithAngle(i,j), &
       wfractPix, ichunk, i, j, dZms(i,j), msFlag(i, j)) !! MS&WS
  call deallocateStormStructData(stormStruct)
  call deallocateDPRProfRet(radarRet)
  call deallocateDPRProfData(radarData)
  
  !This option reads in Liang's mu=2 table
end program
