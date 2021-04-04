SUBROUTINE POD(ni,nj,nk,Nt,Nvar,nmode,var,vol,lambda,phi,at)

  ! !!! INPUT !!!!!!!!
  ! Nemsh: number of mesh (integer)
  ! Nt: number of snapshots (integer)
  ! Nvar: number of variables (integer)
  ! var: snapshots for the variables (double precision)
  ! !!! OUTPUT !!!!!!
  ! lambda: eigen values (double precision)
  ! phi: pod basis (double precision)

  implicit none
  
  integer,intent(in)::ni,nj,nk,Nt,Nvar,nmode
  real*8,intent(in)::var(ni,nj,nk,Nt,Nvar),vol(ni,nj,nk,Nt)
  real*8,intent(out)::lambda(Nt,Nvar),phi(ni,nj,nk,nmode,Nvar),at(Nt,nmode,Nvar)

  ! Local variables
  double precision, parameter::epsilon=1.0d-8
  integer:: i,j,ivar
  double precision::Ct(Nt,Nt),Eval(Nt),Evec(Nt,Nt),tmp(Nt),norm(Nt,Nvar), maxx
  
  do ivar=1,Nvar

     print*, 'Computing correlation matrix for', ivar, 'of', Nvar

     ! Correlation matrix Ct
     do j=1,Nt
        do i=j,Nt
           Ct(i,j)=sum(var(1:ni,1:nj,1:nk,i,ivar)*sqrt(vol(1:ni,1:nj,1:nk,i))*&
                var(1:ni,1:nj,1:nk,j,ivar)*sqrt(vol(1:ni,1:nj,1:nk,j)))
           Ct(j,i)=Ct(i,j)
        enddo
     enddo
     Ct=Ct/Nt
     
     maxx=maxval(Ct(1:Nt,1:Nt))
     if(maxx.gt.0.d0)then
        Ct=Ct/maxx
     endif

     ! Eigen values and eigen vectors    
     Eval=0.d0;Evec=0.d0
     call jacobi(Nt,Ct(1:Nt,1:Nt),Eval(1:Nt),Evec(1:Nt,1:Nt),epsilon)  
     ! arranging eigen values and vectors in descending order  
     do i=1,Nt
        do j=i+1,Nt
           if(Eval(j).gt.Eval(i)) then
              tmp(1)=Eval(j)
              Eval(j)=Eval(i)
              Eval(i)=tmp(1)
              
              tmp(1:Nt)=Evec(1:Nt,j)
              Evec(1:Nt,j)=Evec(1:Nt,i)
              Evec(1:Nt,i)=tmp(1:Nt)
           endif
        enddo
        ! a check for -ve eigen values
        if(Eval(i).lt.0.d0)then
           Eval(i)=0.d0
           Evec(1:Nt,i)=0.d0
        endif
     enddo

     lambda(1:Nt,ivar)=Eval(1:Nt)*maxx
         
     ! POD MODES
     
     do j=1,nmode
        
        phi(1:ni,1:nj,1:nk,j,ivar)=0.d0
        do i=1,Nt
           phi(1:ni,1:nj,1:nk,j,ivar)=phi(1:ni,1:nj,1:nk,j,ivar)+(var(1:ni,1:nj,1:nk,i,ivar)*Evec(i,j))
        enddo
        
        if(lambda(j,ivar).gt.0.d0)then
           phi(1:ni,1:nj,1:nk,j,ivar)=phi(1:ni,1:nj,1:nk,j,ivar)/&
                sqrt(sum(phi(1:ni,1:nj,1:nk,j,ivar)**2*vol(1:ni,1:nj,1:nk,j)))
        else
           phi(1:ni,1:nj,1:nk,j,ivar)=0.d0
        endif
        
        ! POD time coefficients
        do i=1,Nt
           at(i,j,ivar)=sum(phi(1:ni,1:nj,1:nk,j,ivar)*sqrt(vol(1:ni,1:nj,1:nk,j))*&
                sqrt(vol(1:ni,1:nj,1:nk,i))*var(1:ni,1:nj,1:nk,i,ivar))
        enddo
        
     enddo
     
  enddo
  
  
END SUBROUTINE POD



SUBROUTINE  Jacobi(n,a,e,x,abserr)
!===========================================================
! Evaluate eigenvalues and eigenvectors
! of a real symmetric matrix a(n,n): a*x = lambda*x 
! method: Jacoby method for symmetric matrices 
! Alex G. (December 2009)
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - number of equations
! abserr - abs tolerance [sum of (off-diagonal elements)^2]
! output ...
! e(i) - eigenvalues
! x(i,j) - eigenvectors
! comments ...
!===========================================================
implicit none
integer,intent(in)::n
real(8),intent(in)::abserr
real(8),intent(inout)::a(n,n),x(n,n),e(n)
integer::i,j,k
real(8)::b2,bar,beta,coeff,c,s,cs,sc,b2old

! initialize x(i,j)=0, x(i,i)=1
! *** the array operation x=0.0 is specific for Fortran 90/95
x = 0.d0
do i=1,n
  x(i,i) = 1.d0
end do

! find the sum of all off-diagonal elements (squared)
b2 = 0.d0
do i=1,n
   do j=1,n
      if (i.ne.j) b2 = b2 + a(i,j)**2
   end do
end do
! to set a relative error criterian
b2old=b2
if (b2 <= abserr) return

! average for off-diagonal elements /2
bar = 0.5d0*b2/(n*n)

do while ((b2.gt.abserr).and.(b2old.gt.abserr))
   b2old=b2
   do i=1,n-1
      do j=i+1,n
         if (a(j,i)**2 <= bar) cycle  ! do not touch small elements
         b2 = b2 - 2.d0*a(j,i)**2
         bar = 0.5d0*b2/(n*n)
         ! calculate coefficient c and s for Givens matrix
         beta = (a(j,j)-a(i,i))/(2.d0*a(j,i))
         coeff = 0.5d0*beta/sqrt(1.d0+beta**2)
         s = sqrt(max(0.5d0+coeff,0.d0))
         c = sqrt(max(0.5d0-coeff,0.d0))
         ! recalculate rows i and j
         do k=1,n
            cs =  c*a(i,k)+s*a(j,k)
            sc = -s*a(i,k)+c*a(j,k)
            a(i,k) = cs
            a(j,k) = sc
         end do
         ! new matrix a_{k+1} from a_{k}, and eigenvectors 
         do k=1,n
            cs =  c*a(k,i)+s*a(k,j)
            sc = -s*a(k,i)+c*a(k,j)
            a(k,i) = cs
            a(k,j) = sc
            cs =  c*x(k,i)+s*x(k,j)
            sc = -s*x(k,i)+c*x(k,j)
            x(k,i) = cs
            x(k,j) = sc
         end do
      end do
   end do
   b2old=b2old-b2
end do

! rearrange the output
do i=1,n
   e(i)=a(i,i)
end do

END SUBROUTINE Jacobi
