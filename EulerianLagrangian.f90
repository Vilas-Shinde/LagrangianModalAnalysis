SUBROUTINE EUL_to_LAG(ni,nj,nk,nr,grdnw,varnw,delT,fn,Nt,omL,varsn,volsn)
 
  implicit none
  
  integer,intent(in) :: ni,nj,nk,nr,Nt,fn
  real*8,intent(in) :: grdnw(ni,nj,nk,3),varnw(ni,nj,nk,fn,nr),delT
  real*8,intent(inout) :: omL(ni,nj,nk,Nt,3),varsn(ni,nj,nk,Nt,nr),volsn(ni,nj,nk,Nt,2)
  
  integer :: i,j,k,iNt,ii,jj,kk
  real*8  :: omLL(ni,nj,nk,3),varsnL(ni,nj,nk,nr),volsnL(ni,nj,nk)
   
  ! Inital omL
  omL(1:ni,1:nj,1:nk,1,1:3)=grdnw(1:ni,1:nj,1:nk,1:3)
  omLL(1:ni,1:nj,1:nk,1:3)=omL(1:ni,1:nj,1:nk,1,1:3)
  varsn(1:ni,1:nj,1:nk,1,1:nr)=varnw(1:ni,1:nj,1:nk,1,1:nr)
  varsnL(1:ni,1:nj,1:nk,1:nr)=varsn(1:ni,1:nj,1:nk,1,1:nr)
  volsn(1:ni,1:nj,1:nk,1,1)=volsn(1:ni,1:nj,1:nk,1,2)
  volsnL(1:ni,1:nj,1:nk)=volsn(1:ni,1:nj,1:nk,1,1)
  
  i=1;j=1;k=1
  
  do iNt=2,Nt
     ! Calculate omL surface position
     omLL(1:ni,1:nj,1:nk,1:3)=omLL(1:ni,1:nj,1:nk,1:3)+&
          delT*varsnL(1:ni,1:nj,1:nk,1:3)
     ! Calculate velocity on omL surface
     do ii=1,ni
        do jj=1,nj
           do kk=1,nk
              do i=1,ni-1
                 if(omLL(ii,jj,kk,1).ge.grdnw(i,j,k,1).and.omLL(ii,jj,kk,1).lt.grdnw(i+1,j,k,1)) exit
              enddo
              do j=1,nj-1
                 if(omLL(ii,jj,kk,2).ge.grdnw(i,j,k,2).and.omLL(ii,jj,kk,2).lt.grdnw(i,j+1,k,2)) exit
              enddo
              do k=1,nk-1
                 if(omLL(ii,jj,kk,3).ge.grdnw(i,j,k,3).and.omLL(ii,jj,kk,3).lt.grdnw(i,j,k+1,3)) exit
              enddo
              if(fn.eq.1)then
                 varsnL(ii,jj,kk,1:nr)=varnw(i,j,k,1,1:nr)
                 volsnL(ii,jj,kk)=volsn(i,j,k,1,2)
              else
                 varsnL(ii,jj,kk,1:nr)=varnw(i,j,k,iNt,1:nr)
                 volsnL(ii,jj,kk)=volsn(i,j,k,iNt,2)
              endif
           enddo
        enddo
     enddo
     ! Store the omL and varsn
     omL(1:ni,1:nj,1:nk,iNt,1:3)=omLL(1:ni,1:nj,1:nk,1:3)
     varsn(1:ni,1:nj,1:nk,iNt,1:nr)=varsnL(1:ni,1:nj,1:nk,1:nr)
     volsn(1:ni,1:nj,1:nk,iNt,1)=volsnL(1:ni,1:nj,1:nk)/sum(volsnL(1:ni,1:nj,1:nk))
     print*, iNt,Nt
  enddo
 

END SUBROUTINE EUL_to_LAG
