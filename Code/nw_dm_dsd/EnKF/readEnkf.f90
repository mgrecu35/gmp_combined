subroutine openfile(orbNumber)
  integer:: orbNumber
  print*, orbNumber
  call openenkffile(orbNumber)
end subroutine openfile

subroutine closefile()
  call closeenkffile()
end subroutine closefile
  
subroutine read1(nmemb1,rainType,node5,nZka,wfractm,freezH)
  integer, intent(out):: nmemb1,rainType,node5(5),nZka
  real, intent(out):: wfractm,freezH

  call enkfri1(nmemb1)
  call enkfri1(rainType)
  call enkfri5(node5)
  call enkfri1(i)
  call enkfri1(j)
  call enkfri1(nZka)
  call enkfrf(wfractm)
  call enkfrf(freezH)

end subroutine read1

subroutine read2(nmemb1,node5,nZka,zKuObs,iZka,YEns,YObs,&
     zKuEns,log10dnw,rrate,sfcRainEns,dsrtPIAKu,dsrtPIAKa,&
     pia13mod,pia35mod,srtPIAKu,srtRelPIAKu,dsrtRelPIA,tbObs,&
     tbSim,msrflag,nsrflag)
  implicit none
  integer, intent(in):: nmemb1, node5(5),nZka
  real, intent(out):: zKuObs(88)
  integer, intent(out):: iZka(nZka)
  real, intent(out):: Yens(nZKa,nmemb1), YObs(nZKa)
  integer :: k,k1
  real, intent(out):: zKuEns(88,nmemb1),log10dnw(88,nmemb1),&
       rrate(88,nmemb1),sfcRainEns(nmemb1),dsrtPIAKu,dsrtPIAKa,&
       pia13mod(nmemb1),pia35mod(nmemb1)
  real, intent(out) :: srtPIAKu,srtRelPIAKu,dsrtRelPIA,tbObs(13),&
     tbSim(13,nmemb1) 
  integer, intent(out):: msrflag,nsrflag

  do k=node5(1)+1,node5(5)+1
     call enkfrf(zKuObs(k))
  enddo
  do k=1,nZka
     call enkfri1(iZka(k))
     do k1=1,nmemb1
        call enkfrf(Yens(k,k1))
     enddo
     call enkfrf(Yobs(k))
  enddo

  do k=node5(1)+1,node5(5)+1
     do k1=1,nmemb1
        call enkfrf(zKuEns(k,k1))
     enddo
     do k1=1,nmemb1
        call enkfrf(log10dnw(k,k1))
     enddo
     do k1=1,nmemb1
        call enkfrf(rrate(k,k1))
     enddo
  enddo
  do k1=1,nmemb1
     call enkfrf(sfcRainEns(k1))
  enddo
  call enkfrf(dsrtPIAku)
  call enkfrf(dsrtPIAka)
  !print*, dsrtPIAku, dsrtPIAKa
  do k1=1,nmemb1
     call enkfrf(pia13mod(k1))
     call enkfrf(pia35mod(k1))
  enddo
  call enkfrf(srtPIAku)
  call enkfrf(srtRelPIAku)
  call enkfrf(dsrtRelPIA)
  call enkfri1(msrflag)
  call enkfri1(nsrflag)
  do k=1,13
     call enkfrf(tbobs(k))
  enddo
  do k=1,9
     do k1=1,nmemb1
        call enkfrf(tbSim(k,k1))
     enddo
  enddo
  do k=1,4
     do k1=1,nmemb1
        call enkfrf(tbSim(9+k,k1))
     enddo
  enddo

end subroutine read2



subroutine read3(nmemb1,node5,nZka,zKuObs,&
     zKuEns,log10dnw,rrate,sfcRainEns,dsrtPIAKu,dsrtPIAKa,&
     pia13mod,pia35mod,srtPIAKu,srtRelPIAKu,dsrtRelPIA,tbObs,&
     tbSim,msrflag,nsrflag)
  implicit none
  integer, intent(in):: nmemb1, node5(5),nZka
  real, intent(out):: zKuObs(88)
  integer :: k,k1
  real, intent(out):: zKuEns(88,nmemb1),log10dnw(88,nmemb1),&
       rrate(88,nmemb1),sfcRainEns(nmemb1),dsrtPIAKu,dsrtPIAKa,&
       pia13mod(nmemb1),pia35mod(nmemb1)

  integer, intent(out):: msrflag,nsrflag

  real, intent(out) :: srtPIAKu,srtRelPIAKu,dsrtRelPIA,tbObs(13),&
     tbSim(13,nmemb1) 

  do k=node5(1)+1,node5(5)+1
     call enkfrf(zKuObs(k))
  enddo

  do k=node5(1)+1,node5(5)+1
     do k1=1,nmemb1
        call enkfrf(zKuEns(k,k1))
     enddo
     do k1=1,nmemb1
        call enkfrf(log10dnw(k,k1))
     enddo
     do k1=1,nmemb1
        call enkfrf(rrate(k,k1))
     enddo
  enddo
  do k1=1,nmemb1
     call enkfrf(sfcRainEns(k1))
  enddo
  call enkfrf(dsrtPIAku)
  call enkfrf(dsrtPIAka)
  !print*, dsrtPIAku, dsrtPIAKa
  do k1=1,nmemb1
     call enkfrf(pia13mod(k1))
     call enkfrf(pia35mod(k1))
  enddo
  
    call enkfrf(srtPIAku)
  call enkfrf(srtRelPIAku)
  call enkfrf(dsrtRelPIA)
  call enkfri1(msrflag)
  call enkfri1(nsrflag)

  do k=1,13
     call enkfrf(tbobs(k))
  enddo
  do k=1,9
     do k1=1,nmemb1
        call enkfrf(tbSim(k,k1))
     enddo
  enddo
  do k=1,4
     do k1=1,nmemb1
        call enkfrf(tbSim(k,k1))
     enddo
  enddo

end subroutine read3
