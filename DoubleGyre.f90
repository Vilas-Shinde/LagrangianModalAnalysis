SUBROUTINE DOUBLE_GYRE(fn,dt,nr,case,plane)
  
  implicit none

  integer, intent(in) :: fn,nr
  double precision, intent(in) :: dt
  character(len=1000),intent(in) :: case,plane

  character(len=1000) :: fname(10),command,ch
  integer :: i,j,k,fi
  integer :: ni,nj,nk,ng
  double precision, allocatable :: grd(:,:,:,:), sol(:,:,:,:)
  double precision, parameter :: pi=atan(1.0d0)*4.0d0
  double precision ::  xmin,xmax,ymin,ymax,zcoord,dx,dy
  double precision :: A,ep,om,tau
  double precision :: f, dfdx,d2fdx2
  
  command='mkdir -p '//trim(case)
  call system(command)
  command='rm -rf '//trim(case)//'/'//trim(plane)//'*'
  call  system(command)

  command=trim(plane)//'list.txt'  
  fname(1)=command
  open(10,file=trim(fname(1)))
  ! Grid inputs 
  ng=1; ni=601; nj=301; nk=1
   
  ! Stream function psi(x,y,t) = A sin(pi f(x,t)) cos(pi y)
  ! f(x,t) = ep sin(om t) x^2 + x - 2 ep sin(om t)x
  ! u = - pi A sin(pi f(x,t)) cos(pi y)
  ! v = pi A cos(pi f(x,t)) sin(pi y) df/dx
  ! Model parameters
  A = 0.1; ep = 0.1; om=2.0*pi/10.0 
  
  ! GRID
  ! xmin, xmax, ymin, ymax, zcoord
  xmin=0.0; xmax=2.0; ymin=0.0; ymax=1.0; zcoord=0.0
  dx=(xmax-xmin)/(ni-1);   dy=(ymax-ymin)/(nj-1)
  allocate(grd(ni,nj,nk,3))
  do k=1,nk
     do j=1,nj
        do i=1,ni
           grd(i,j,k,1)=xmin+(i-1)*dx
           grd(i,j,k,2)=ymin+(j-1)*dy
           grd(i,j,k,3)=zcoord
        enddo
     enddo
  enddo
  ! SOLUTION
  allocate(sol(ni,nj,nk,nr))
  
  do fi=1,fn
     tau=(fi-1)*dt
     do k=1,nk
        do j=1,nj
           do i=1,ni
              ! u
              sol(i,j,k,1)=-pi*A*sin(pi*f(ep,om,tau,grd(i,j,k,1)))*cos(pi*grd(i,j,k,2))
              ! v
              sol(i,j,k,2)=pi*A*cos(pi*f(ep,om,tau,grd(i,j,k,1)))*sin(pi*grd(i,j,k,2))*dfdx(ep,om,tau,grd(i,j,k,1))
              ! Abs velocity
              sol(i,j,k,3)=sqrt(sol(i,j,k,1)**2+sol(i,j,k,2)**2)
              ! Vorticity  
              sol(i,j,k,4)=pi*A*sin(pi*grd(i,j,k,2))*(cos(pi*f(ep,om,tau,grd(i,j,k,1)))*&
                   d2fdx2(ep,om,tau,grd(i,j,k,1))-pi*dfdx(ep,om,tau,grd(i,j,k,1))**2*sin(pi*f(ep,om,tau,grd(i,j,k,1))))&
                   -pi**2*A*sin(pi*f(ep,om,tau,grd(i,j,k,1)))*sin(pi*grd(i,j,k,2))
           enddo
        enddo
     enddo
     
     write(ch,*) fi
     fname(3)=trim(plane)//'-'//trim(adjustl(ch))//'.x'
     fname(2)=trim(plane)//'-'//trim(adjustl(ch))//'.q'
     
     write(10,'(a)') trim(fname(2))                     !!!! solution file
     write(10,'(a)') trim(fname(3))                     !!!! grid file
     print*, fi, fn, trim(fname(2)),' ' ,trim(fname(3))
     
!!!!!!!!!!!!!!!!!!!!!!!! writing  grid
     fname(4)=trim(case)//'/'//trim(fname(3))
     OPEN(11,FILE=trim(fname(4)),FORM='UNFORMATTED',ACTION='WRITE',CONVERT='LITTLE_ENDIAN',STATUS='REPLACE')
     WRITE(11) ng
     WRITE(11) ni, nj, nk
     WRITE(11) grd
     CLOSE(11)
!!!!!!!!!!!!!!!!!!!!!!!! reading sol
     fname(4)=trim(case)//'/'//trim(fname(2))
     OPEN(11,FILE=trim(fname(4)),FORM='UNFORMATTED',ACTION='WRITE',CONVERT='LITTLE_ENDIAN',STATUS='REPLACE')
     WRITE(11) ng
     WRITE(11) ni, nj, nk, nr
     WRITE(11) sol
     CLOSE(11)
     
  enddo
  
  close(10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  deallocate(grd,sol)
  
END SUBROUTINE  DOUBLE_GYRE

DOUBLE PRECISION FUNCTION f(ep,om,t,x)

  implicit none
  
  double precision, intent(in) :: ep, om, t, x
  
  f=ep*sin(om*t)*x**2+x-2.0*ep*sin(om*t)*x
  
END FUNCTION f

DOUBLE PRECISION FUNCTION dfdx(ep,om,t,x)

  implicit none

  double precision, intent(in) :: ep, om, t, x
  
  dfdx=2.*ep*sin(om*t)*x+1.0-2.0*ep*sin(om*t)
  
END FUNCTION dfdx

DOUBLE PRECISION FUNCTION d2fdx2(ep,om,t,x)

  implicit none
  
  double precision, intent(in) :: ep, om, t, x
  
  d2fdx2=2.*ep*sin(om*t)

END FUNCTION d2fdx2
