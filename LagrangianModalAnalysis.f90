! Lagrangian modal analysis in terms of proper orthogonal decomposition (POD), considering a simple double-gyre flow 
PROGRAM LMA
  
  IMPLICIT NONE
  
  integer :: bs,stream
  character(len=100) :: xqfile,fdl3difile
  integer :: ni,nj,nk,n,ng,nr
  character :: scrap
  real*8, parameter :: gamma=1.4
  real*8, allocatable :: grd(:,:,:,:),sol(:,:,:,:)
  real*8 :: dt,dtau
  real*8, allocatable :: grdnw(:,:,:,:),varnw(:,:,:,:,:),grdnwr(:,:,:,:)
  
  integer :: i,iNt,Nt,j,k
  integer :: fnn,fn,fi,td
  real*8, allocatable :: omL(:,:,:,:,:),varsn(:,:,:,:,:),volsn(:,:,:,:,:)
  character(len=1000) :: case, plane, fname(4), ch, command, pd
  
  integer :: nmode
  real*8, allocatable :: vol(:,:,:),lam(:,:),phi(:,:,:,:,:),psi(:,:,:),lambda(:,:,:)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Eulerian Inputs
  case='EULsnap'
  plane='snap'
! number of snapshots
  fn=101
! Eulerian time step
  dt=0.15d0
! Number of variables
  nr=4

!!!!!!!!!!!!!!!!!!!!! Lagrangian inputs
  ! Number of snapshots
  Nt=fn
  ! Lagrangian time-step
  dtau=dt
  ! Time direction: 1 -> forward, -1 -> backward
  td=1
  ! Reduced number of LPOD modes for postprocessing
  nmode=min(Nt/2,200)
  
  print*,' '
  print*, 'Generating Eulerian snapshots of Double Gyre flow.....START'
  print*,' '
  ! Double gyre flow
  
  CALL DOUBLE_GYRE(fn,dt,nr,case,plane)
  
  print*,' '
  print*, 'Generating Eulerian snapshots for Double Gyre flow.....END'
  print*,' '
  
  command=trim(plane)//'list.txt' 
  fname(1)=command
  open(10,file=trim(fname(1)))  

  print*,' '
  print*, 'Reading Eulerian snapshots...........'
  print*,' '

  do fi=1,fn
     read(10,*) fname(2)                     !!!! solution file
     read(10,*) fname(3)                     !!!! grid file
     print*, fi, fn, trim(fname(2)),' ' ,trim(fname(3))   
!!!!!!!!!!!!!!!!!!!!!!!! reading grid
     fname(4)=trim(case)//'/'//trim(fname(3))
     OPEN(11,FILE=trim(fname(4)),FORM='UNFORMATTED',ACTION='READ',CONVERT='LITTLE_ENDIAN')
     READ(11) ng
     READ(11) ni, nj, nk
     if(.not.allocated(grd))then
        allocate(grd(ni,nj,nk,3))
     endif
     READ(11) grd
     CLOSE(11)
!!!!!!!!!!!!!!!!!!!!!!!! reading sol
     fname(4)=trim(case)//'/'//trim(fname(2))
     OPEN(11,FILE=trim(fname(4)),FORM='UNFORMATTED',ACTION='READ',CONVERT='LITTLE_ENDIAN')
     READ(11) ng
     READ(11) ni, nj, nk
     if(.not.allocated(sol))then
        allocate(sol(ni,nj,nk,nr))
     endif
     READ(11) sol
     CLOSE(11)
!!!!!!!!!!!!!!!!!!!!!!!!!! Variable conversion if any
     if(.not.allocated(grdnw))then
        allocate(grdnw(ni,nj,nk,3))
     endif
     ! grid
     grdnw(:,:,:,1:3)=grd(:,:,:,1:3)
     
     if(.not.allocated(varnw))then
        allocate(varnw(ni,nj,nk,fn,nr))
     endif
     ! u
     varnw(:,:,:,fi,1)=sol(:,:,:,1)
     ! v
     varnw(:,:,:,fi,2)=sol(:,:,:,2)
     ! absolute velocity
     varnw(:,:,:,fi,3)=sol(:,:,:,3)
     ! vorticity
     varnw(:,:,:,fi,4)=sol(:,:,:,4)
  enddo
  close(10)
  deallocate(grd,sol)

  print*,' '
  print*, 'Reading Eulerian snapshots...........DONE'
  print*,' '

    
  allocate(grdnwr(ni,nj,nk,3))
  grdnwr(:,:,:,:)=grdnw(:,:,:,:)

  ! Compute the weight matrix for a non-uniform grid
  allocate(vol(ni,nj,nk),volsn(ni,nj,nk,Nt,2))
  vol=1.d0; volsn=0.d0
  call vol_wt(ni,nj,nk,grdnw(:,:,:,:),vol)
  vol=vol/sum(vol(:,:,:))
  do fi=1,Nt
     volsn(:,:,:,fi,2)=vol(:,:,:)
  enddo
  deallocate(vol)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Lagrangian snapshots

  case='LAGsnap'
  
  allocate(omL(ni,nj,nk,Nt,3),varsn(ni,nj,nk,Nt,nr))
  omL=0.d0; varsn=0.d0

  print*,' '
  print*, 'Generating Lagrangian snapshots.....START'
  print*,' '

  ! Eulerian snapshots to Lagrangian snapshots
  CALL EUL_to_LAG(ni,nj,nk,nr,grdnw,varnw,td*dtau,fn,Nt,omL,varsn,volsn)
  deallocate(grdnw,varnw)
  
  ! Output directories
  command='mkdir -p '//trim(case)
  call system(command)
  if(td.gt.0)then
     command='mkdir -p '//trim(case)//'/'//'forward'
     call system(command)
     command='rm -rf '//trim(case)//'/'//'forward'//trim(plane)//'*'
     call  system(command)
  else
     command='mkdir -p '//trim(case)//'/'//'backward'
     call system(command)
     command='rm -rf '//trim(case)//'/'//'backward'//trim(plane)//'*'
     call  system(command)
  endif

  ! Lagrangian snapshots
  do i=1,Nt
     call  wtec(ng,ni,nj,nk,grdnwr(:,:,:,:),varsn(:,:,:,i,1:nr),nr,i,case,plane,td)
  enddo

  print*,' '
  print*, 'Generating Lagrangian snapshots.....END'
  print*,' '


!!!!!!!!!!!!!!!!!!!! Lagrangian Proper Orthogonal Decomposition

  print*,' '
  print*, 'Computing LPOD modes.................'
  print*,' '  
  
  case='LPOD'
  
  command='mkdir -p '//trim(case)
  call system(command)
  
  if(td.gt.0)then
     command='mkdir -p '//trim(case)//'/'//'forward'
     call system(command)
     command='rm -rf '//trim(case)//'/'//'forward'//trim(plane)//'*'
     call  system(command)
  else
     command='mkdir -p '//trim(case)//'/'//'backward'
     call system(command)
     command='rm -rf '//trim(case)//'/'//'backward'//trim(plane)//'*'
     call  system(command)
  endif
  
  allocate(lam(Nt,nr),phi(ni,nj,nk,nmode,nr),psi(Nt,nmode,nr))
  CALL POD(ni,nj,nk,Nt,nr,nmode,varsn,volsn(:,:,:,:,1),lam,phi,psi)

  do i=1,nr
!!!! Psi
     write(ch,*) i
     if(td.gt.0)then
        fname(1)=trim(case)//'/forward/'//trim(plane)//'_'//trim(adjustl(ch))//'_psi.txt'
        fname(2)=trim(case)//'/forward/'//trim(plane)//'_'//trim(adjustl(ch))//'_lam.txt'
     else
        fname(1)=trim(case)//'/backward/'//trim(plane)//'_'//trim(adjustl(ch))//'_psi.txt'
        fname(2)=trim(case)//'/backward/'//trim(plane)//'_'//trim(adjustl(ch))//'_lam.txt'
     endif
     open(10,file=fname(1));  open(11,file=fname(2))
     do iNt=1,Nt
        write(10,100) iNt, iNt*dtau, psi(iNt,1:nmode,i)
        if(iNt.le.nmode)then
           write(11,100) iNt, lam(iNt,i), sum(psi(:,iNt,i)**2)/Nt
        else
           write(11,100) iNt, lam(iNt,i)
        endif
     enddo
     close(10);close(11)
!!!! Phi
     call wtec(ng,ni,nj,nk,grdnwr(:,:,:,1:3),phi(:,:,:,1:nmode,i),nmode,i,case,plane,td)
  enddo
  deallocate(lam,phi,psi)

  print*,' '
  print*, 'Computing LPOD modes.................DONE'
  print*,' '
    
  deallocate(volsn)
  deallocate(omL,varsn)
  deallocate(grdnwr)
   
100 format (i8,30000(3x,es16.8e3))
  
END PROGRAM LMA


SUBROUTINE wtec(ng,ni,nj,nk,grd,sol,nr,num,case,plane,td)

  implicit none
  
  integer, intent(in) :: ng,ni,nj,nk,nr,num,td
  double precision, intent(in) :: grd(ni,nj,nk,3), sol(ni,nj,nk,nr)
  character(len=1000), intent(in) :: case, plane

  integer :: bs,stream,n
  character(len=1000) :: ch,fgrd,fsol

  write(ch,*) num
  if(td.gt.0)then
     fgrd=trim(case)//'/forward/'//trim(plane)//'_'//trim(adjustl(ch))//'.x'
     fsol=trim(case)//'/forward/'//trim(plane)//'_'//trim(adjustl(ch))//'.q'
  else
     fgrd=trim(case)//'/backward/'//trim(plane)//'_'//trim(adjustl(ch))//'.x'
     fsol=trim(case)//'/backward/'//trim(plane)//'_'//trim(adjustl(ch))//'.q'
  endif
  
  OPEN(10,FILE=TRIM(fgrd),FORM='UNFORMATTED',ACTION='WRITE',CONVERT='LITTLE_ENDIAN',STATUS='REPLACE')
  write(10) ng
  DO n = 1, ng
     write(10) ni, nj, nk
  ENDDO
  DO n=1,ng
     write(10) grd(1:ni,1:nj,1:nk,1:3)
  ENDDO
  close(10)

  OPEN(10,FILE=TRIM(fsol),FORM='UNFORMATTED',ACTION='WRITE',CONVERT='LITTLE_ENDIAN',STATUS='REPLACE')
  write(10) ng
  DO n = 1, ng
     write(10) ni, nj, nk, nr
  ENDDO
  DO n=1,ng
     write(10) sol(1:ni,1:nj,1:nk,1:nr)
  ENDDO
  close(10)

END SUBROUTINE wtec


SUBROUTINE vol_wt(ni,nj,nk,grdnw,vol)
  
  implicit none
  
  integer,intent(in)::ni,nj,nk
  real*8,intent(in)::grdnw(ni,nj,nk,3)
  real*8,intent(inout)::vol(ni,nj,nk)
  
  integer :: i,j,k
  real*8 :: dx,dy,dz,dxx,dyy,dzz
  dx=1.d0;dy=1.d0;dz=1.d0
  dxx=1.d0;dyy=1.d0;dzz=1.d0
  if(ni.eq.1)then
     i=ni
     do j=1,nj-1
        do k=1,nk-1
           dyy=abs(grdnw(i,j+1,k,2)-grdnw(i,j,k,2))
           dzz=abs(grdnw(i,j+1,k,3)-grdnw(i,j,k,3))
           dy=sqrt(dyy**2+dzz**2)
           
           dyy=abs(grdnw(i,j,k+1,2)-grdnw(i,j,k,2))
           dzz=abs(grdnw(i,j,k+1,3)-grdnw(i,j,k,3))
           dz=sqrt(dyy**2+dzz**2)
           vol(i,j,k)=dx*dy*dz
        enddo
        vol(i,j,k)=dx*dy*dz
     enddo
     vol(i,j,1:nk)=dx*dy*dz
  elseif(nj.eq.1)then
     j=nj
     do i=1,ni-1
        do k=1,nk-1
           dxx=abs(grdnw(i+1,j,k,1)-grdnw(i,j,k,1))
           dzz=abs(grdnw(i+1,j,k,3)-grdnw(i,j,k,3))
           dx=sqrt(dxx**2+dzz**2)
           
           dxx=abs(grdnw(i,j,k+1,1)-grdnw(i,j,k,1))
           dzz=abs(grdnw(i,j,k+1,3)-grdnw(i,j,k,3))
           dz=sqrt(dxx**2+dzz**2)
           vol(i,j,k)=dx*dy*dz
        enddo
        vol(i,j,k)=dx*dy*dz
     enddo
     vol(i,j,1:nk)=dx*dy*dz
  elseif(nk.eq.1)then
     k=nk
     do i=1,ni-1
        do j=1,nj-1
           dxx=abs(grdnw(i+1,j,k,1)-grdnw(i,j,k,1))
           dyy=abs(grdnw(i+1,j,k,2)-grdnw(i,j,k,2))
           dx=sqrt(dxx**2+dyy**2)

           dxx=abs(grdnw(i,j+1,k,1)-grdnw(i,j,k,1))
           dyy=abs(grdnw(i,j+1,k,2)-grdnw(i,j,k,2))
           dy=sqrt(dxx**2+dyy**2)
           vol(i,j,k)=dx*dy*dz
        enddo
        vol(i,j,k)=dx*dy*dz
     enddo
     vol(i,1:nj,k)=dx*dy*dz
  else
     do i=1,ni-1
        do j=1,nj-1
           do k=1,nk-1
              dxx=abs(grdnw(i+1,j,k,1)-grdnw(i,j,k,1))
              dyy=abs(grdnw(i+1,j,k,2)-grdnw(i,j,k,2))
              dx=sqrt(dxx**2+dyy**2)
              
              dyy=abs(grdnw(i,j+1,k,2)-grdnw(i,j,k,2))
              dzz=abs(grdnw(i,j+1,k,3)-grdnw(i,j,k,3))
              dy=sqrt(dyy**2+dzz**2)

              dzz=abs(grdnw(i,j,k+1,3)-grdnw(i,j,k,3))
              dxx=abs(grdnw(i,j,k+1,1)-grdnw(i,j,k,1))
              dz=sqrt(dzz**2+dxx**2)
              vol(i,j,k)=dx*dy*dz
           enddo
           vol(i,j,k)=dx*dy*dz
        enddo
        vol(i,j,1:nk)=dx*dy*dz
     enddo
     vol(i,1:nj,1:nk)=dx*dy*dz
  endif

END SUBROUTINE vol_wt
