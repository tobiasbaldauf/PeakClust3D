program pksampler
implicit none
INTEGER(4), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
!input correlators
integer, parameter :: nr=40
real(8), dimension(nr) :: r,xik0j0,xik2j0,xik4j0,xik2j2,xik4j2,xik4j4,xik1j1,xik3j1,xik3j3
real(8) :: sigma0a,sigma1a,sigma2a
!
integer(4), dimension(nr) :: cnt
real(8), dimension(nr) :: xipk,xipktmp

real(8), parameter :: pi=3.14159265d0

real(8)  :: xmin
integer :: i,nrand
real(8) :: nu1,nu2,abun
real(4) :: u,x1,x2,nubar1,nubar2,nusq1,nusq2

logical :: xweight,xcut

!upcrossing constraint, x*sigma2 weight in abundance and Tr(H) for clustering
xweight=.False.


!load input correlators and variances
open(20,file='sigs.dat',status='old')
	read(20,*) sigma0a,sigma1a,sigma2a
close(20)


open(20,file='corrs.dat',status='old')
do i=1,nr
	read(20,*) r(i),xik0j0(i),xik2j0(i),xik4j0(i),xik2j2(i),xik4j2(i),xik4j4(i),xik1j1(i),xik3j1(i),xik3j3(i)
enddo
close(20)

!correlation for fixed peak height difference
nu1=2.33
nu2=2.57
nrand=5000000


xmin=0.0d0
abun=abunf(nu1,xweight)*abunf(nu2,xweight)
print*, abun,abunf(nu1,xweight),abunf(nu2,xweight)

xcut=.False.
xmin=0.0d0
abun=abunf(nu1,xweight)*abunf(nu2,xweight)
print*, abun


call sampnr(nu1,nu2,nrand,cnt,xipk,.False.,xweight,xcut)


open(20,file='pkcorr.dat',status='replace')
do i=1,nr
	write(20,'(f5.1,es15.4,i10,es15.5)') r(i),xipk(i)/abun-1.0d0,cnt(i),dble(cnt(i))/dble(nrand)
enddo
close(20)

stop
!lognormal sample of peak heights
xipktmp=0.0d0
xipk=0.0d0

nubar1=0
nubar2=0
nusq1=0
nusq2=0

do i=1,400
	print*,i/400.0d0
	call random_number(x1)
	call random_number(x2)
	nu1=exp(dble(sqrt(-2.0*log(x1))*cos(2.0*pi*x2))*sqrt(0.022)+0.552)/sigma0a
	nu2=exp(dble(sqrt(-2.0*log(x1))*sin(2.0*pi*x2))*sqrt(0.022)+0.552)/sigma0a
	!nu1=dble(sqrt(-2.0*log(x1))*cos(2.0*pi*x2))*0.25+1.95
	!nu2=dble(sqrt(-2.0*log(x1))*sin(2.0*pi*x2))*0.25+1.95		
	
	print*, nu1,nu2
	nubar1=nubar1+nu1
	nubar2=nubar2+nu2
	nusq1=nusq1+nu1**2.0
	nusq2=nusq2+nu2**2.0
	abun=abun+abunf(nu1,xweight)*abunf(nu2,xweight)
	call sampnr(nu1,nu2,5000,cnt,xipktmp,.True.,xweight,xcut)
	xipk=xipk+xipktmp
enddo

print*, nubar1/400.0
print*, nubar2/400.0
print*, sqrt(nusq1/400.0-(nubar1/400.0)**2.0)
print*, sqrt(nusq2/400.0-(nubar2/400.0)**2.0)



open(20,file='pkcorr_lognorm_mnu0p55_snu0p022.dat',status='replace')
do i=1,nr
	write(20,'(f5.1,es15.4,i10,es15.5)') r(i),xipk(i)/abun-1.0d0,cnt(i),dble(cnt(i))/dble(nrand)
enddo
close(20)

stop

!gaussian sample of peak heights
xipktmp=0.0d0
xipk=0.0d0
do i=1,400
	print*,i/400.0d0
	call random_number(x1)
	call random_number(x2)
	nu1=dble(sqrt(-2.0*log(x1))*cos(2.0*pi*x2))*0.4+2.6
	nu2=dble(sqrt(-2.0*log(x1))*sin(2.0*pi*x2))*0.4+2.6
	!nu1=dble(sqrt(-2.0*log(x1))*cos(2.0*pi*x2))*0.25+1.95
	!nu2=dble(sqrt(-2.0*log(x1))*sin(2.0*pi*x2))*0.25+1.95		
	
	print*, nu1,nu2
	abun=abun+abunf(nu1,xweight)*abunf(nu2,xweight)
	call sampnr(nu1,nu2,5000,cnt,xipktmp,.True.,xweight,xcut)
	xipk=xipk+xipktmp
enddo


open(20,file='pkcorr_gauss_mnu2p6_snu0p4.dat',status='replace')
do i=1,nr
	write(20,'(f5.1,es15.4,i10,es15.5)') r(i),xipk(i)/abun-1.0d0,cnt(i),dble(cnt(i))/dble(nrand)
enddo
close(20)



!Threshold Sample
xipktmp=0.0d0
xipk=0.0d0
do i=1,400
	print*,i/400.0d0
	call random_number(u)
	nu1=2.4d0+3.0d0*u
	call random_number(u)
	nu2=2.4d0+3.0d0*u
	!print*, nu1,nu2
	abun=abun+abunf(nu1,xweight)*abunf(nu2,xweight)
	call sampnr(nu1,nu2,5000,cnt,xipktmp,.True.,xweight,xcut)
	xipk=xipk+xipktmp
enddo


open(20,file='pkcorr_thresh.dat',status='replace')
do i=1,nr
	write(20,'(f5.1,es15.4,i10,es15.5)') r(i),xipk(i)/abun-1.0d0,cnt(i),dble(cnt(i))/dble(nrand)
enddo
close(20)

!=======================================================================================
!=======================================================================================
!=======================================================================================
contains

function abunf(nu,xweight)
implicit none
real(8) :: nu
real(8) :: abunf
real(8) :: xmax,xmin,dx
integer :: i,nn
logical :: xweight

nn=1000
xmax=10
xmin=0
dx=(xmax-xmin)/dble(nn)
abunf=0.0d0
do i=1,nn
	if (xweight) then
		abunf=abunf+wfull(nu,(i-0.5d0)*dx)*dx*(i-0.5d0)*dx*sigma2a
	else
		abunf=abunf+wfull(nu,(i-0.5d0)*dx)*dx
	endif
enddo


end function abunf


!=======================================================================================
function wfull(nu,x)
implicit none
real(8) :: nu,x,gamma,da,pre
real(8) :: wfull

gamma=sigma1a**2/sigma0a/sigma2a
da=(Sigma1a**6*Sigma2a**6*(-60*Sigma1a**4*Sigma2a**4+60*Sigma0a**2*Sigma2a**6))/dble(922640625)
pre=(8*Pi**2*2.0/8.0*Sigma2a**9)/sqrt((2.0*pi)**10.0*da)
!print*,pre
wfull=pre*(2*exp((-5*x**2)/2. - Nu**2/2. - (x - gamma*Nu)**2/(2.*(1 - gamma**2)))*(-32 + 32*exp((15*x**2)/8.) + 10*x**2 + 155*exp((15*x**2)/8.)*x**2 + 5*exp((5*x**2)/2.)*Sqrt(10*Pi)*x*(-3 + x**2)*Erf((Sqrt(2.5)*x)/2.) + 5*exp((5*x**2)/2.)*Sqrt(10*Pi)*x*(-3 + x**2)*Erf(Sqrt(2.5)*x)))/455625.
if (xcut) then
if (x>3.0*xmin) then
wfull=pre*(exp(-nu**2/2. - (-(gamma*nu) + x)**2/(2.*(1 - gamma**2)) -(5*(9*xmin**2 + x**2))/2.)*(exp((15*xmin*x)/4.)*(-8*exp((45*xmin*x)/4.)*(16 +90*xmin**2 - 15*xmin*x - 5*x**2) + exp((15*(9*xmin**2 + x**2))/8.)*(128 +6075*xmin**4 - 14175*xmin**3*x + 620*x**2 - 15*xmin*x*(92 + 135*x**2) +45*xmin**2*(-32 + 225*x**2))) - 20*exp((5*(9*xmin**2 +x**2))/2.)*Sqrt(10*Pi)*x*(-3 + x**2)*Erf((Sqrt(2.5)*(3*xmin - x))/2.) -20*exp((5*(9*xmin**2 + x**2))/2.)*Sqrt(10*Pi)*x*(-3 +x**2)*Erf(Sqrt(2.5)*(3*xmin - x))))/911250.
else
wfull=0.0d0
endif
endif
end function wfull


!=======================================================================================
subroutine sampnr(nu1,nu2,nrand,cnt,xipk,silent,xweight,xcut)
!calculates the pk-pk correlation function
implicit none
integer :: i,j,nrand
real(8) :: nu1,nu2
real(8), dimension(:) :: xipk
integer, dimension(:) :: cnt

real(8), dimension(1,1) :: norm
real(8), dimension(8,8) :: ainvmat
real(8), dimension(12,8) :: bamat
real(8), dimension(12,12) :: vinvmat
real(8), dimension(8,1) :: Yh
real(8) :: da
real(8), dimension(12,1) :: w,mu
logical :: silent,xweight,xcut

real(4), dimension(3,3) :: mat,evmat
real(4), dimension(3) :: evals
integer :: nrot

cnt=0
xipk=0.0d0
do i=1,nr
	
	call setmat(i,ainvmat,bamat,vinvmat,da)
	call choldc(vinvmat)
	yh=reshape((/nu1*sigma0a,0.0d0,0.0d0,0.0d0,nu2*sigma0a,0.0d0,0.0d0,0.0d0/),(/8,1/))
	mu=matmul(BAMat,Yh)
	norm=1.0d0/sqrt((2.0d0*pi)**8*da)*exp(-1.0/2.0*matmul(matmul(transpose(Yh),ainvmat),Yh))

	do j=1,nrand
		call normalrandvec(w)
		w=mu+matmul(vinvmat,w)		
		if (w(7,1)>0.) cycle
		if (w(1,1)>0.) cycle
		if ((w(7,1)*w(8,1)-w(10,1)**2)<0.) cycle
		if ((w(1,1)*w(2,1)-w(4,1)**2)<0.) cycle
		if (det(w(7:12,1))>0.) cycle
		if (det(w(1:6,1))>0.) cycle
		if (xcut .eqv. .True.) then
			mat(1,:)=(/w(1,1),w(4,1),w(5,1)/)
			mat(2,:)=(/w(4,1),w(2,1),w(6,1)/)
			mat(3,:)=(/w(6,1),w(6,1),w(3,1)/)
			call jacobi(mat,evals,evmat,nrot)
			if (minval(abs(evals))<xmin*sigma2a) cycle
			mat(1,:)=(/w(7,1),w(10,1),w(11,1)/)
			mat(2,:)=(/w(10,1),w(8,1),w(12,1)/)
			mat(3,:)=(/w(11,1),w(12,1),w(9,1)/)	
			call jacobi(mat,evals,evmat,nrot)
			if (minval(abs(evals))<xmin*sigma2a) cycle
		endif
		if (xweight) then
			xipk(i)=xipk(i)+abs(det(w(1:6,1))*det(w(7:12,1)))*(w(7,1)+w(8,1)+w(9,1))*(w(1,1)+w(2,1)+w(3,1))
		else
			xipk(i)=xipk(i)+abs(det(w(1:6,1))*det(w(7:12,1)))
		endif
		cnt(i)=cnt(i)+1
	enddo
	if (silent .eqv. .False.) write(*,'(f5.1,es15.4,i10,es15.5)') r(i),norm(1,1)*xipk(i)/dble(nrand)/abun-1.0d0,cnt(i),dble(cnt(i))/dble(nrand)

	xipk(i)=norm(1,1)*xipk(i)/dble(nrand)
enddo

return
end subroutine sampnr

!=======================================================================================
SUBROUTINE choldc(a)
!Given an N × N positive-definite symmetric matrix a, this routine constructs its Cholesky decomposition,
!A = L · LT . On input, only the upper triangle of a need be given; it is not modified. The Cholesky factor L
!is returned in the lower triangle of a, except for its diagonal elements, which are returned in p, a vector of
!length N.


IMPLICIT NONE
REAL(8), DIMENSION(:,:) :: a
REAL(8), DIMENSION(:), allocatable :: p
INTEGER(4) :: i,n
REAL(8) :: summ
n=size(a,1)
allocate(p(n))
do i=1,n
	summ=a(i,i)-dot_product(a(i,1:i-1),a(i,1:i-1))
	if (summ <= 0.0) print*, ('choldc failed')
	p(i)=sqrt(summ)
	a(i+1:n,i)=(a(i,i+1:n)-matmul(a(i+1:n,1:i-1),a(i,1:i-1)))/p(i)
end do

do i=1,n
	a(i,i)=p(i)
	a(i,i+1:n)=0.0d0
enddo
deallocate(p)

END SUBROUTINE choldc

!=======================================================================================
function det(w)
implicit none
real(8), dimension(:) :: w
real(8) :: det

det=w(1)*w(2)*w(3)-w(3)*w(4)**2-w(2)*w(5)**2+2.*w(4)*w(5)*w(6)-w(1)*w(6)**2

return 
end function det

!=======================================================================================
subroutine setmat(i,ainvmat,bamat,vinvmat,da)
implicit none
real(8), dimension(:,:) :: ainvmat
real(8), dimension(:,:) :: bamat
real(8), dimension(:,:) :: vinvmat
real(8) :: da
integer :: i

da=((xik2j0(i)**2 + 2*xik2j0(i)*xik2j2(i) + xik2j2(i)**2 - sigma1a**4)**2*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))/729.

bamat(1,1)=-(27*xik1j1(i)**3*(xik3j1(i) + xik3j3(i)) + 9*xik1j1(i)*(xik3j1(i) + xik3j3(i))*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + 15*xik1j1(i)**2*(xik2j0(i)**2 - xik2j0(i)*xik2j2(i) - 2*xik2j2(i)**2 - sigma1a**4) + 5*(xik0j0(i)*(xik2j0(i) + xik2j2(i)) - sigma0a**2*sigma1a**2)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4))/(15.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
bamat(1,2)=0.0d0
bamat(1,3)=0.0d0
bamat(1,4)=(-15*xik1j1(i)**3*(xik2j0(i) + xik2j2(i)) + 9*xik0j0(i)*xik1j1(i)**2*(xik3j1(i) + xik3j3(i)) + 3*(xik2j0(i) - 2*xik2j2(i))*(xik3j1(i) + xik3j3(i))*(xik0j0(i)**2 - sigma0a**4) - 5*xik1j1(i)*((-2*xik2j0(i) + xik2j2(i))*sigma0a**2*sigma1a**2 + xik0j0(i)*(xik2j0(i)**2 - xik2j0(i)*xik2j2(i) - 2*xik2j2(i)**2 + sigma1a**4)))/(5.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
bamat(1,5)=(5*xik2j0(i)**3*sigma0a**2 + 20*xik2j2(i)**3*sigma0a**2 - 20*xik0j0(i)*xik2j2(i)**2*sigma1a**2 - 5*xik2j0(i)**2*(3*xik2j2(i)*sigma0a**2 + xik0j0(i)*sigma1a**2) + xik0j0(i)*sigma1a**2*(-9*xik1j1(i)*(xik3j1(i) + xik3j3(i)) + 5*sigma1a**4) + xik2j2(i)*(-18*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 45*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4) + xik2j0(i)*(9*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 20*xik0j0(i)*xik2j2(i)*sigma1a**2 - 5*sigma0a**2*sigma1a**4))/(15.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
bamat(1,6)=0.0d0
bamat(1,7)=0.0d0
bamat(1,8)=-(9*xik1j1(i)**2*(xik3j1(i) + xik3j3(i))*sigma0a**2 - 15*xik1j1(i)**3*sigma1a**2 + 3*(xik3j1(i) + xik3j3(i))*(xik0j0(i)**2 - sigma0a**4)*sigma1a**2 + 5*xik1j1(i)*(xik2j0(i)**2*sigma0a**2 - 2*xik2j2(i)**2*sigma0a**2 + xik0j0(i)*xik2j2(i)*sigma1a**2 + sigma0a**2*sigma1a**4 - xik2j0(i)*(xik2j2(i)*sigma0a**2 + 2*xik0j0(i)*sigma1a**2)))/(5.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
bamat(2,1)=-(27*xik1j1(i)**3*(xik3j1(i) + xik3j3(i)) + 9*xik1j1(i)*(xik3j1(i) + xik3j3(i))*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + 15*xik1j1(i)**2*(xik2j0(i)**2 - xik2j0(i)*xik2j2(i) - 2*xik2j2(i)**2 - sigma1a**4) + 5*(xik0j0(i)*(xik2j0(i) + xik2j2(i)) - sigma0a**2*sigma1a**2)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4))/(15.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
bamat(2,2)=0.0d0
bamat(2,3)=0.0d0
bamat(2,4)=(-15*xik1j1(i)**3*(xik2j0(i) + xik2j2(i)) + 9*xik0j0(i)*xik1j1(i)**2*(xik3j1(i) + xik3j3(i)) + 3*(xik2j0(i) - 2*xik2j2(i))*(xik3j1(i) + xik3j3(i))*(xik0j0(i)**2 - sigma0a**4) - 5*xik1j1(i)*((-2*xik2j0(i) + xik2j2(i))*sigma0a**2*sigma1a**2 + xik0j0(i)*(xik2j0(i)**2 - xik2j0(i)*xik2j2(i) - 2*xik2j2(i)**2 + sigma1a**4)))/(5.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
bamat(2,5)=(5*xik2j0(i)**3*sigma0a**2 + 20*xik2j2(i)**3*sigma0a**2 - 20*xik0j0(i)*xik2j2(i)**2*sigma1a**2 - 5*xik2j0(i)**2*(3*xik2j2(i)*sigma0a**2 + xik0j0(i)*sigma1a**2) + xik0j0(i)*sigma1a**2*(-9*xik1j1(i)*(xik3j1(i) + xik3j3(i)) + 5*sigma1a**4) + xik2j2(i)*(-18*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 45*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4) + xik2j0(i)*(9*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 20*xik0j0(i)*xik2j2(i)*sigma1a**2 - 5*sigma0a**2*sigma1a**4))/(15.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
bamat(2,6)=0.0d0
bamat(2,7)=0.0d0
bamat(2,8)=-(9*xik1j1(i)**2*(xik3j1(i) + xik3j3(i))*sigma0a**2 - 15*xik1j1(i)**3*sigma1a**2 + 3*(xik3j1(i) + xik3j3(i))*(xik0j0(i)**2 - sigma0a**4)*sigma1a**2 + 5*xik1j1(i)*(xik2j0(i)**2*sigma0a**2 - 2*xik2j2(i)**2*sigma0a**2 + xik0j0(i)*xik2j2(i)*sigma1a**2 + sigma0a**2*sigma1a**4 - xik2j0(i)*(xik2j2(i)*sigma0a**2 + 2*xik0j0(i)*sigma1a**2)))/(5.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
bamat(3,1)=-((3*xik1j1(i)**2 + xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2)*(5*xik2j0(i)**2 - 20*xik2j0(i)*xik2j2(i) + 20*xik2j2(i)**2 + 27*xik1j1(i)*xik3j1(i) - 18*xik1j1(i)*xik3j3(i) - 5*sigma1a**4))/(15.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
bamat(3,2)=0.0d0
bamat(3,3)=0.0d0
bamat(3,4)=(-15*xik1j1(i)**3*(xik2j0(i) - 2*xik2j2(i)) + 9*xik0j0(i)*xik1j1(i)**2*(3*xik3j1(i) - 2*xik3j3(i)) + 3*(xik2j0(i) - 2*xik2j2(i))*(3*xik3j1(i) - 2*xik3j3(i))*(xik0j0(i)**2 - sigma0a**4) - 5*xik1j1(i)*(-2*(xik2j0(i) - 2*xik2j2(i))*sigma0a**2*sigma1a**2 + xik0j0(i)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 + sigma1a**4)))/(5.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
bamat(3,5)=((xik2j0(i)*sigma0a**2 - 2*xik2j2(i)*sigma0a**2 - xik0j0(i)*sigma1a**2)*(5*xik2j0(i)**2 - 20*xik2j0(i)*xik2j2(i) + 20*xik2j2(i)**2 + 27*xik1j1(i)*xik3j1(i) - 18*xik1j1(i)*xik3j3(i) - 5*sigma1a**4))/(15.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
bamat(3,6)=0.0d0
bamat(3,7)=0.0d0
bamat(3,8)=-(9*xik1j1(i)**2*(3*xik3j1(i) - 2*xik3j3(i))*sigma0a**2 - 15*xik1j1(i)**3*sigma1a**2 + 3*(3*xik3j1(i) - 2*xik3j3(i))*(xik0j0(i)**2 - sigma0a**4)*sigma1a**2 + 5*xik1j1(i)*(xik2j0(i)**2*sigma0a**2 + 4*xik2j2(i)**2*sigma0a**2 + 4*xik0j0(i)*xik2j2(i)*sigma1a**2 + sigma0a**2*sigma1a**4 - 2*xik2j0(i)*(2*xik2j2(i)*sigma0a**2 + xik0j0(i)*sigma1a**2)))/(5.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
bamat(4,1)=0.0d0
bamat(4,2)=0.0d0
bamat(4,3)=0.0d0
bamat(4,4)=0.0d0
bamat(4,5)=0.0d0
bamat(4,6)=0.0d0
bamat(4,7)=0.0d0
bamat(4,8)=0.0d0
bamat(5,1)=0.0d0
bamat(5,2)=(3*(xik2j0(i) + xik2j2(i))*(xik3j1(i) + xik3j3(i)))/(5.*(xik2j0(i)**2 + 2*xik2j0(i)*xik2j2(i) + xik2j2(i)**2 - sigma1a**4))
bamat(5,3)=0.0d0
bamat(5,4)=0.0d0
bamat(5,5)=0.0d0
bamat(5,6)=(-3*(xik3j1(i) + xik3j3(i))*sigma1a**2)/(5.*(xik2j0(i)**2 + 2*xik2j0(i)*xik2j2(i) + xik2j2(i)**2 - sigma1a**4))
bamat(5,7)=0.0d0
bamat(5,8)=0.0d0
bamat(6,1)=0.0d0
bamat(6,2)=0.0d0
bamat(6,3)=(3*(xik2j0(i) + xik2j2(i))*(xik3j1(i) + xik3j3(i)))/(5.*(xik2j0(i)**2 + 2*xik2j0(i)*xik2j2(i) + xik2j2(i)**2 - sigma1a**4))
bamat(6,4)=0.0d0
bamat(6,5)=0.0d0
bamat(6,6)=0.0d0
bamat(6,7)=(-3*(xik3j1(i) + xik3j3(i))*sigma1a**2)/(5.*(xik2j0(i)**2 + 2*xik2j0(i)*xik2j2(i) + xik2j2(i)**2 - sigma1a**4))
bamat(6,8)=0.0d0
bamat(7,1)=(5*xik2j0(i)**3*sigma0a**2 + 20*xik2j2(i)**3*sigma0a**2 - 20*xik0j0(i)*xik2j2(i)**2*sigma1a**2 - 5*xik2j0(i)**2*(3*xik2j2(i)*sigma0a**2 + xik0j0(i)*sigma1a**2) + xik0j0(i)*sigma1a**2*(-9*xik1j1(i)*(xik3j1(i) + xik3j3(i)) + 5*sigma1a**4) + xik2j2(i)*(-18*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 45*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4) + xik2j0(i)*(9*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 20*xik0j0(i)*xik2j2(i)*sigma1a**2 - 5*sigma0a**2*sigma1a**4))/(15.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
bamat(7,2)=0.0d0
bamat(7,3)=0.0d0
bamat(7,4)=(9*xik1j1(i)**2*(xik3j1(i) + xik3j3(i))*sigma0a**2 - 15*xik1j1(i)**3*sigma1a**2 + 3*(xik3j1(i) + xik3j3(i))*(xik0j0(i)**2 - sigma0a**4)*sigma1a**2 + 5*xik1j1(i)*(xik2j0(i)**2*sigma0a**2 - 2*xik2j2(i)**2*sigma0a**2 + xik0j0(i)*xik2j2(i)*sigma1a**2 + sigma0a**2*sigma1a**4 - xik2j0(i)*(xik2j2(i)*sigma0a**2 + 2*xik0j0(i)*sigma1a**2)))/(5.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
bamat(7,5)=-(27*xik1j1(i)**3*(xik3j1(i) + xik3j3(i)) + 9*xik1j1(i)*(xik3j1(i) + xik3j3(i))*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + 15*xik1j1(i)**2*(xik2j0(i)**2 - xik2j0(i)*xik2j2(i) - 2*xik2j2(i)**2 - sigma1a**4) + 5*(xik0j0(i)*(xik2j0(i) + xik2j2(i)) - sigma0a**2*sigma1a**2)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4))/(15.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
bamat(7,6)=0.0d0
bamat(7,7)=0.0d0
bamat(7,8)=(15*xik1j1(i)**3*(xik2j0(i) + xik2j2(i)) - 9*xik0j0(i)*xik1j1(i)**2*(xik3j1(i) + xik3j3(i)) - 3*(xik2j0(i) - 2*xik2j2(i))*(xik3j1(i) + xik3j3(i))*(xik0j0(i)**2 - sigma0a**4) + 5*xik1j1(i)*((-2*xik2j0(i) + xik2j2(i))*sigma0a**2*sigma1a**2 + xik0j0(i)*(xik2j0(i)**2 - xik2j0(i)*xik2j2(i) - 2*xik2j2(i)**2 + sigma1a**4)))/(5.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
bamat(8,1)=(5*xik2j0(i)**3*sigma0a**2 + 20*xik2j2(i)**3*sigma0a**2 - 20*xik0j0(i)*xik2j2(i)**2*sigma1a**2 - 5*xik2j0(i)**2*(3*xik2j2(i)*sigma0a**2 + xik0j0(i)*sigma1a**2) + xik0j0(i)*sigma1a**2*(-9*xik1j1(i)*(xik3j1(i) + xik3j3(i)) + 5*sigma1a**4) + xik2j2(i)*(-18*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 45*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4) + xik2j0(i)*(9*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 20*xik0j0(i)*xik2j2(i)*sigma1a**2 - 5*sigma0a**2*sigma1a**4))/(15.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
bamat(8,2)=0.0d0
bamat(8,3)=0.0d0
bamat(8,4)=(9*xik1j1(i)**2*(xik3j1(i) + xik3j3(i))*sigma0a**2 - 15*xik1j1(i)**3*sigma1a**2 + 3*(xik3j1(i) + xik3j3(i))*(xik0j0(i)**2 - sigma0a**4)*sigma1a**2 + 5*xik1j1(i)*(xik2j0(i)**2*sigma0a**2 - 2*xik2j2(i)**2*sigma0a**2 + xik0j0(i)*xik2j2(i)*sigma1a**2 + sigma0a**2*sigma1a**4 - xik2j0(i)*(xik2j2(i)*sigma0a**2 + 2*xik0j0(i)*sigma1a**2)))/(5.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
bamat(8,5)=-(27*xik1j1(i)**3*(xik3j1(i) + xik3j3(i)) + 9*xik1j1(i)*(xik3j1(i) + xik3j3(i))*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + 15*xik1j1(i)**2*(xik2j0(i)**2 - xik2j0(i)*xik2j2(i) - 2*xik2j2(i)**2 - sigma1a**4) + 5*(xik0j0(i)*(xik2j0(i) + xik2j2(i)) - sigma0a**2*sigma1a**2)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4))/(15.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
bamat(8,6)=0.0d0
bamat(8,7)=0.0d0
bamat(8,8)=(15*xik1j1(i)**3*(xik2j0(i) + xik2j2(i)) - 9*xik0j0(i)*xik1j1(i)**2*(xik3j1(i) + xik3j3(i)) - 3*(xik2j0(i) - 2*xik2j2(i))*(xik3j1(i) + xik3j3(i))*(xik0j0(i)**2 - sigma0a**4) + 5*xik1j1(i)*((-2*xik2j0(i) + xik2j2(i))*sigma0a**2*sigma1a**2 + xik0j0(i)*(xik2j0(i)**2 - xik2j0(i)*xik2j2(i) - 2*xik2j2(i)**2 + sigma1a**4)))/(5.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
bamat(9,1)=((xik2j0(i)*sigma0a**2 - 2*xik2j2(i)*sigma0a**2 - xik0j0(i)*sigma1a**2)*(5*xik2j0(i)**2 - 20*xik2j0(i)*xik2j2(i) + 20*xik2j2(i)**2 + 27*xik1j1(i)*xik3j1(i) - 18*xik1j1(i)*xik3j3(i) - 5*sigma1a**4))/(15.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
bamat(9,2)=0.0d0
bamat(9,3)=0.0d0
bamat(9,4)=(9*xik1j1(i)**2*(3*xik3j1(i) - 2*xik3j3(i))*sigma0a**2 - 15*xik1j1(i)**3*sigma1a**2 + 3*(3*xik3j1(i) - 2*xik3j3(i))*(xik0j0(i)**2 - sigma0a**4)*sigma1a**2 + 5*xik1j1(i)*(xik2j0(i)**2*sigma0a**2 + 4*xik2j2(i)**2*sigma0a**2 + 4*xik0j0(i)*xik2j2(i)*sigma1a**2 + sigma0a**2*sigma1a**4 - 2*xik2j0(i)*(2*xik2j2(i)*sigma0a**2 + xik0j0(i)*sigma1a**2)))/(5.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
bamat(9,5)=-((3*xik1j1(i)**2 + xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2)*(5*xik2j0(i)**2 - 20*xik2j0(i)*xik2j2(i) + 20*xik2j2(i)**2 + 27*xik1j1(i)*xik3j1(i) - 18*xik1j1(i)*xik3j3(i) - 5*sigma1a**4))/(15.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
bamat(9,6)=0.0d0
bamat(9,7)=0.0d0
bamat(9,8)=(15*xik1j1(i)**3*(xik2j0(i) - 2*xik2j2(i)) - 9*xik0j0(i)*xik1j1(i)**2*(3*xik3j1(i) - 2*xik3j3(i)) - 3*(xik2j0(i) - 2*xik2j2(i))*(3*xik3j1(i) - 2*xik3j3(i))*(xik0j0(i)**2 - sigma0a**4) + 5*xik1j1(i)*(-2*(xik2j0(i) - 2*xik2j2(i))*sigma0a**2*sigma1a**2 + xik0j0(i)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 + sigma1a**4)))/(5.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
bamat(10,1)=0.0d0
bamat(10,2)=0.0d0
bamat(10,3)=0.0d0
bamat(10,4)=0.0d0
bamat(10,5)=0.0d0
bamat(10,6)=0.0d0
bamat(10,7)=0.0d0
bamat(10,8)=0.0d0
bamat(11,1)=0.0d0
bamat(11,2)=(3*(xik3j1(i) + xik3j3(i))*sigma1a**2)/(5.*(xik2j0(i)**2 + 2*xik2j0(i)*xik2j2(i) + xik2j2(i)**2 - sigma1a**4))
bamat(11,3)=0.0d0
bamat(11,4)=0.0d0
bamat(11,5)=0.0d0
bamat(11,6)=(-3*(xik2j0(i) + xik2j2(i))*(xik3j1(i) + xik3j3(i)))/(5.*(xik2j0(i)**2 + 2*xik2j0(i)*xik2j2(i) + xik2j2(i)**2 - sigma1a**4))
bamat(11,7)=0.0d0
bamat(11,8)=0.0d0
bamat(12,1)=0.0d0
bamat(12,2)=0.0d0
bamat(12,3)=(3*(xik3j1(i) + xik3j3(i))*sigma1a**2)/(5.*(xik2j0(i)**2 + 2*xik2j0(i)*xik2j2(i) + xik2j2(i)**2 - sigma1a**4))
bamat(12,4)=0.0d0
bamat(12,5)=0.0d0
bamat(12,6)=0.0d0
bamat(12,7)=(-3*(xik2j0(i) + xik2j2(i))*(xik3j1(i) + xik3j3(i)))/(5.*(xik2j0(i)**2 + 2*xik2j0(i)*xik2j2(i) + xik2j2(i)**2 - sigma1a**4))
bamat(12,8)=0.0d0
ainvmat(1,1)=(-(xik2j0(i)**2*sigma0a**2) + 4*xik2j0(i)*xik2j2(i)*sigma0a**2 - 4*xik2j2(i)**2*sigma0a**2 - 3*xik1j1(i)**2*sigma1a**2 + sigma0a**2*sigma1a**4)/(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4))
ainvmat(1,2)=0.0d0
ainvmat(1,3)=0.0d0
ainvmat(1,4)=(3*xik1j1(i)*(-(xik2j0(i)*sigma0a**2) + 2*xik2j2(i)*sigma0a**2 + xik0j0(i)*sigma1a**2))/(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4))
ainvmat(1,5)=(3*xik1j1(i)**2*(xik2j0(i) - 2*xik2j2(i)) + xik0j0(i)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4))/(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4))
ainvmat(1,6)=0.0d0
ainvmat(1,7)=0.0d0
ainvmat(1,8)=(3*xik1j1(i)*(-3*xik1j1(i)**2 - xik0j0(i)*xik2j0(i) + 2*xik0j0(i)*xik2j2(i) + sigma0a**2*sigma1a**2))/(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4))
ainvmat(2,1)=0.0d0
ainvmat(2,2)=(-3*sigma1a**2)/(xik2j0(i)**2 + 2*xik2j0(i)*xik2j2(i) + xik2j2(i)**2 - sigma1a**4)
ainvmat(2,3)=0.0d0
ainvmat(2,4)=0.0d0
ainvmat(2,5)=0.0d0
ainvmat(2,6)=(3*(xik2j0(i) + xik2j2(i)))/(xik2j0(i)**2 + 2*xik2j0(i)*xik2j2(i) + xik2j2(i)**2 - sigma1a**4)
ainvmat(2,7)=0.0d0
ainvmat(2,8)=0.0d0
ainvmat(3,1)=0.0d0
ainvmat(3,2)=0.0d0
ainvmat(3,3)=(-3*sigma1a**2)/(xik2j0(i)**2 + 2*xik2j0(i)*xik2j2(i) + xik2j2(i)**2 - sigma1a**4)
ainvmat(3,4)=0.0d0
ainvmat(3,5)=0.0d0
ainvmat(3,6)=0.0d0
ainvmat(3,7)=(3*(xik2j0(i) + xik2j2(i)))/(xik2j0(i)**2 + 2*xik2j0(i)*xik2j2(i) + xik2j2(i)**2 - sigma1a**4)
ainvmat(3,8)=0.0d0
ainvmat(4,1)=(3*xik1j1(i)*(-(xik2j0(i)*sigma0a**2) + 2*xik2j2(i)*sigma0a**2 + xik0j0(i)*sigma1a**2))/(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4))
ainvmat(4,2)=0.0d0
ainvmat(4,3)=0.0d0
ainvmat(4,4)=(-9*xik1j1(i)**2*sigma0a**2 + 3*(-xik0j0(i)**2 + sigma0a**4)*sigma1a**2)/(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4))
ainvmat(4,5)=(3*xik1j1(i)*(3*xik1j1(i)**2 + xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2))/(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4))
ainvmat(4,6)=0.0d0
ainvmat(4,7)=0.0d0
ainvmat(4,8)=(9*xik0j0(i)*xik1j1(i)**2 + 3*xik0j0(i)**2*(xik2j0(i) - 2*xik2j2(i)) - 3*(xik2j0(i) - 2*xik2j2(i))*sigma0a**4)/(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4))
ainvmat(5,1)=(3*xik1j1(i)**2*(xik2j0(i) - 2*xik2j2(i)) + xik0j0(i)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4))/(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4))
ainvmat(5,2)=0.0d0
ainvmat(5,3)=0.0d0
ainvmat(5,4)=(3*xik1j1(i)*(3*xik1j1(i)**2 + xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2))/(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4))
ainvmat(5,5)=(-(xik2j0(i)**2*sigma0a**2) + 4*xik2j0(i)*xik2j2(i)*sigma0a**2 - 4*xik2j2(i)**2*sigma0a**2 - 3*xik1j1(i)**2*sigma1a**2 + sigma0a**2*sigma1a**4)/(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4))
ainvmat(5,6)=0.0d0
ainvmat(5,7)=0.0d0
ainvmat(5,8)=(3*xik1j1(i)*(xik2j0(i)*sigma0a**2 - 2*xik2j2(i)*sigma0a**2 - xik0j0(i)*sigma1a**2))/(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4))
ainvmat(6,1)=0.0d0
ainvmat(6,2)=(3*(xik2j0(i) + xik2j2(i)))/(xik2j0(i)**2 + 2*xik2j0(i)*xik2j2(i) + xik2j2(i)**2 - sigma1a**4)
ainvmat(6,3)=0.0d0
ainvmat(6,4)=0.0d0
ainvmat(6,5)=0.0d0
ainvmat(6,6)=(-3*sigma1a**2)/(xik2j0(i)**2 + 2*xik2j0(i)*xik2j2(i) + xik2j2(i)**2 - sigma1a**4)
ainvmat(6,7)=0.0d0
ainvmat(6,8)=0.0d0
ainvmat(7,1)=0.0d0
ainvmat(7,2)=0.0d0
ainvmat(7,3)=(3*(xik2j0(i) + xik2j2(i)))/(xik2j0(i)**2 + 2*xik2j0(i)*xik2j2(i) + xik2j2(i)**2 - sigma1a**4)
ainvmat(7,4)=0.0d0
ainvmat(7,5)=0.0d0
ainvmat(7,6)=0.0d0
ainvmat(7,7)=(-3*sigma1a**2)/(xik2j0(i)**2 + 2*xik2j0(i)*xik2j2(i) + xik2j2(i)**2 - sigma1a**4)
ainvmat(7,8)=0.0d0
ainvmat(8,1)=(3*xik1j1(i)*(-3*xik1j1(i)**2 - xik0j0(i)*xik2j0(i) + 2*xik0j0(i)*xik2j2(i) + sigma0a**2*sigma1a**2))/(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4))
ainvmat(8,2)=0.0d0
ainvmat(8,3)=0.0d0
ainvmat(8,4)=(9*xik0j0(i)*xik1j1(i)**2 + 3*xik0j0(i)**2*(xik2j0(i) - 2*xik2j2(i)) - 3*(xik2j0(i) - 2*xik2j2(i))*sigma0a**4)/(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4))
ainvmat(8,5)=(3*xik1j1(i)*(xik2j0(i)*sigma0a**2 - 2*xik2j2(i)*sigma0a**2 - xik0j0(i)*sigma1a**2))/(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4))
ainvmat(8,6)=0.0d0
ainvmat(8,7)=0.0d0
ainvmat(8,8)=(-9*xik1j1(i)**2*sigma0a**2 + 3*(-xik0j0(i)**2 + sigma0a**4)*sigma1a**2)/(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4))
vinvmat(1,1)=(25*xik2j0(i)**4*sigma0a**2 + 100*xik2j2(i)**4*sigma0a**2 + 81*xik1j1(i)**2*xik3j1(i)**2*sigma0a**2 + 162*xik1j1(i)**2*xik3j1(i)*xik3j3(i)*sigma0a**2 + 81*xik1j1(i)**2*xik3j3(i)**2*sigma0a**2 - 200*xik0j0(i)*xik2j2(i)**3*sigma1a**2 - 270*xik1j1(i)**3*xik3j1(i)*sigma1a**2 + 27*xik0j0(i)**2*xik3j1(i)**2*sigma1a**2 - 270*xik1j1(i)**3*xik3j3(i)*sigma1a**2 + 54*xik0j0(i)**2*xik3j1(i)*xik3j3(i)*sigma1a**2 + 27*xik0j0(i)**2*xik3j3(i)**2*sigma1a**2 - 27*xik3j1(i)**2*sigma0a**4*sigma1a**2 - 54*xik3j1(i)*xik3j3(i)*sigma0a**4*sigma1a**2 - 27*xik3j3(i)**2*sigma0a**4*sigma1a**2 + 90*xik1j1(i)*xik3j1(i)*sigma0a**2*sigma1a**4 + 90*xik1j1(i)*xik3j3(i)*sigma0a**2*sigma1a**4 + 75*xik1j1(i)**2*sigma1a**6 - 25*sigma0a**2*sigma1a**8 - 50*xik2j0(i)**3*(xik2j2(i)*sigma0a**2 + xik0j0(i)*sigma1a**2) + 405*xik1j1(i)**4*sigma2a**2 - 270*xik1j1(i)**2*sigma0a**2*sigma1a**2*sigma2a**2 - 45*xik0j0(i)**2*sigma1a**4*sigma2a**2 + 45*sigma0a**4*sigma1a**4*sigma2a**2 + 10*xik0j0(i)*xik2j2(i)*(9*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma1a**2 + 5*sigma1a**6 - 54*xik1j1(i)**2*sigma2a**2) + 15*xik2j2(i)**2*(-12*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 25*xik1j1(i)**2*sigma1a**2 + 5*sigma0a**2*sigma1a**4 + 12*xik0j0(i)**2*sigma2a**2 - 12*sigma0a**4*sigma2a**2) - 15*xik2j0(i)**2*(5*xik2j2(i)**2*sigma0a**2 - 6*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 5*xik1j1(i)**2*sigma1a**2 - 10*xik0j0(i)*xik2j2(i)*sigma1a**2 + 3*(-xik0j0(i)**2 + sigma0a**4)*sigma2a**2) + 10*xik2j0(i)*(10*xik2j2(i)**3*sigma0a**2 + xik0j0(i)*(-18*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 27*xik1j1(i)**2*sigma2a**2) + 3*xik2j2(i)*(-3*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 10*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4 - 6*xik0j0(i)**2*sigma2a**2 + 6*sigma0a**4*sigma2a**2)))/(225.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(1,2)=(25*xik2j0(i)**4*sigma0a**2 + 100*xik2j2(i)**4*sigma0a**2 + 81*xik1j1(i)**2*xik3j1(i)**2*sigma0a**2 + 162*xik1j1(i)**2*xik3j1(i)*xik3j3(i)*sigma0a**2 + 81*xik1j1(i)**2*xik3j3(i)**2*sigma0a**2 - 200*xik0j0(i)*xik2j2(i)**3*sigma1a**2 - 270*xik1j1(i)**3*xik3j1(i)*sigma1a**2 + 27*xik0j0(i)**2*xik3j1(i)**2*sigma1a**2 - 270*xik1j1(i)**3*xik3j3(i)*sigma1a**2 + 54*xik0j0(i)**2*xik3j1(i)*xik3j3(i)*sigma1a**2 + 27*xik0j0(i)**2*xik3j3(i)**2*sigma1a**2 - 27*xik3j1(i)**2*sigma0a**4*sigma1a**2 - 54*xik3j1(i)*xik3j3(i)*sigma0a**4*sigma1a**2 - 27*xik3j3(i)**2*sigma0a**4*sigma1a**2 + 90*xik1j1(i)*xik3j1(i)*sigma0a**2*sigma1a**4 + 90*xik1j1(i)*xik3j3(i)*sigma0a**2*sigma1a**4 + 75*xik1j1(i)**2*sigma1a**6 - 25*sigma0a**2*sigma1a**8 - 50*xik2j0(i)**3*(xik2j2(i)*sigma0a**2 + xik0j0(i)*sigma1a**2) + 135*xik1j1(i)**4*sigma2a**2 - 90*xik1j1(i)**2*sigma0a**2*sigma1a**2*sigma2a**2 - 15*xik0j0(i)**2*sigma1a**4*sigma2a**2 + 15*sigma0a**4*sigma1a**4*sigma2a**2 + 10*xik0j0(i)*xik2j2(i)*(9*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma1a**2 + 5*sigma1a**6 - 18*xik1j1(i)**2*sigma2a**2) + 15*xik2j2(i)**2*(-12*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 25*xik1j1(i)**2*sigma1a**2 + 5*sigma0a**2*sigma1a**4 + 4*xik0j0(i)**2*sigma2a**2 - 4*sigma0a**4*sigma2a**2) - 15*xik2j0(i)**2*(5*xik2j2(i)**2*sigma0a**2 - 6*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 5*xik1j1(i)**2*sigma1a**2 - 10*xik0j0(i)*xik2j2(i)*sigma1a**2 + (-xik0j0(i)**2 + sigma0a**4)*sigma2a**2) + 10*xik2j0(i)*(10*xik2j2(i)**3*sigma0a**2 + xik0j0(i)*(-18*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 9*xik1j1(i)**2*sigma2a**2) + 3*xik2j2(i)*(-3*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 10*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4 - 2*xik0j0(i)**2*sigma2a**2 + 2*sigma0a**4*sigma2a**2)))/(225.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(1,3)=(25*xik2j0(i)**4*sigma0a**2 - 200*xik2j2(i)**4*sigma0a**2 + 243*xik1j1(i)**2*xik3j1(i)**2*sigma0a**2 + 81*xik1j1(i)**2*xik3j1(i)*xik3j3(i)*sigma0a**2 - 162*xik1j1(i)**2*xik3j3(i)**2*sigma0a**2 + 100*xik0j0(i)*xik2j2(i)**3*sigma1a**2 - 540*xik1j1(i)**3*xik3j1(i)*sigma1a**2 + 81*xik0j0(i)**2*xik3j1(i)**2*sigma1a**2 + 135*xik1j1(i)**3*xik3j3(i)*sigma1a**2 + 27*xik0j0(i)**2*xik3j1(i)*xik3j3(i)*sigma1a**2 - 54*xik0j0(i)**2*xik3j3(i)**2*sigma1a**2 - 81*xik3j1(i)**2*sigma0a**4*sigma1a**2 - 27*xik3j1(i)*xik3j3(i)*sigma0a**4*sigma1a**2 + 54*xik3j3(i)**2*sigma0a**4*sigma1a**2 + 180*xik1j1(i)*xik3j1(i)*sigma0a**2*sigma1a**4 - 45*xik1j1(i)*xik3j3(i)*sigma0a**2*sigma1a**4 + 75*xik1j1(i)**2*sigma1a**6 - 25*sigma0a**2*sigma1a**8 - 25*xik2j0(i)**3*(5*xik2j2(i)*sigma0a**2 + 2*xik0j0(i)*sigma1a**2) + 135*xik1j1(i)**4*sigma2a**2 - 90*xik1j1(i)**2*sigma0a**2*sigma1a**2*sigma2a**2 - 15*xik0j0(i)**2*sigma1a**4*sigma2a**2 + 15*sigma0a**4*sigma1a**4*sigma2a**2 - 5*xik0j0(i)*xik2j2(i)*(-9*xik1j1(i)*(7*xik3j1(i) + 2*xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 36*xik1j1(i)**2*sigma2a**2) - 30*xik2j2(i)**2*(3*xik1j1(i)*(xik3j1(i) - 4*xik3j3(i))*sigma0a**2 + 10*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4 - 2*xik0j0(i)**2*sigma2a**2 + 2*sigma0a**4*sigma2a**2) + 15*xik2j0(i)**2*(10*xik2j2(i)**2*sigma0a**2 + 3*xik1j1(i)*(4*xik3j1(i) - xik3j3(i))*sigma0a**2 - 5*xik1j1(i)**2*sigma1a**2 + 15*xik0j0(i)*xik2j2(i)*sigma1a**2 + (xik0j0(i)**2 - sigma0a**4)*sigma2a**2) + 5*xik2j0(i)*(20*xik2j2(i)**3*sigma0a**2 - 60*xik0j0(i)*xik2j2(i)**2*sigma1a**2 + 2*xik0j0(i)*(9*xik1j1(i)*(-4*xik3j1(i) + xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 9*xik1j1(i)**2*sigma2a**2) + 3*xik2j2(i)*(-3*xik1j1(i)*(7*xik3j1(i) + 2*xik3j3(i))*sigma0a**2 + 20*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4 - 4*xik0j0(i)**2*sigma2a**2 + 4*sigma0a**4*sigma2a**2)))/(225.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(1,4)=0.0d0
vinvmat(1,5)=0.0d0
vinvmat(1,6)=0.0d0
vinvmat(1,7)=(-1890*xik1j1(i)**3*(xik2j0(i) + xik2j2(i))*(xik3j1(i) + xik3j3(i)) + 405*xik1j1(i)**4*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + 9*xik0j0(i)**2*(-42*xik2j2(i)*(xik3j1(i) + xik3j3(i))**2 + 5*xik2j0(i)**2*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + 20*xik2j2(i)**2*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + xik2j0(i)*(21*xik3j1(i)**2 + 42*xik3j1(i)*xik3j3(i) + 21*xik3j3(i)**2 - 20*xik2j2(i)*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))) - 5*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma1a**4) - 175*xik0j0(i)*(xik2j0(i)**4 - 2*xik2j0(i)**3*xik2j2(i) - 3*xik2j0(i)**2*xik2j2(i)**2 + 4*xik2j0(i)*xik2j2(i)**3 + 4*xik2j2(i)**4 - 6*xik2j0(i)*xik2j2(i)*sigma1a**4 + 3*xik2j2(i)**2*sigma1a**4 - sigma1a**8) - 630*xik1j1(i)*(xik3j1(i) + xik3j3(i))*((-2*xik2j0(i) + xik2j2(i))*sigma0a**2*sigma1a**2 + xik0j0(i)*(xik2j0(i)**2 - xik2j0(i)*xik2j2(i) - 2*xik2j2(i)**2 + sigma1a**4)) - 3*xik1j1(i)**2*(175*xik2j0(i)**3 - 350*xik2j2(i)**3 - 9*(21*xik0j0(i)*(xik3j1(i) + xik3j3(i))**2 - 10*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma0a**2*sigma1a**2) + 20*xik2j2(i)*(9*xik0j0(i)*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) - 35*sigma1a**4) - 5*xik2j0(i)*(105*xik2j2(i)**2 + 18*xik0j0(i)*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + 35*sigma1a**4)) + sigma0a**2*(-180*xik2j2(i)**2*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma0a**2 + 350*xik2j0(i)**3*sigma1a**2 + 1400*xik2j2(i)**3*sigma1a**2 + 45*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma0a**2*sigma1a**4 - 15*xik2j0(i)**2*(21*xik4j0(i)*sigma0a**2 + 30*xik4j2(i)*sigma0a**2 + 9*xik4j4(i)*sigma0a**2 + 70*xik2j2(i)*sigma1a**2) + xik2j0(i)*(-189*xik3j1(i)**2*sigma0a**2 - 378*xik3j1(i)*xik3j3(i)*sigma0a**2 - 189*xik3j3(i)**2*sigma0a**2 + 1260*xik2j2(i)*xik4j0(i)*sigma0a**2 + 1800*xik2j2(i)*xik4j2(i)*sigma0a**2 + 540*xik2j2(i)*xik4j4(i)*sigma0a**2 - 350*sigma1a**6) + 14*xik2j2(i)*(27*xik3j1(i)**2*sigma0a**2 + 54*xik3j1(i)*xik3j3(i)*sigma0a**2 + 27*xik3j3(i)**2*sigma0a**2 - 25*sigma1a**6)))/(1575.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(1,8)=(-1890*xik1j1(i)**3*(xik2j0(i) + xik2j2(i))*(xik3j1(i) + xik3j3(i)) + 135*xik1j1(i)**4*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + 3*xik0j0(i)**2*(-126*xik2j2(i)*(xik3j1(i) + xik3j3(i))**2 + 5*xik2j0(i)**2*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + 20*xik2j2(i)**2*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + xik2j0(i)*(63*xik3j1(i)**2 + 126*xik3j1(i)*xik3j3(i) + 63*xik3j3(i)**2 - 20*xik2j2(i)*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))) - 5*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma1a**4) - 175*xik0j0(i)*(xik2j0(i)**4 - 2*xik2j0(i)**3*xik2j2(i) - 3*xik2j0(i)**2*xik2j2(i)**2 + 4*xik2j0(i)*xik2j2(i)**3 + 4*xik2j2(i)**4 - 6*xik2j0(i)*xik2j2(i)*sigma1a**4 + 3*xik2j2(i)**2*sigma1a**4 - sigma1a**8) - 630*xik1j1(i)*(xik3j1(i) + xik3j3(i))*((-2*xik2j0(i) + xik2j2(i))*sigma0a**2*sigma1a**2 + xik0j0(i)*(xik2j0(i)**2 - xik2j0(i)*xik2j2(i) - 2*xik2j2(i)**2 + sigma1a**4)) - 3*xik1j1(i)**2*(175*xik2j0(i)**3 - 350*xik2j2(i)**3 + 3*(-63*xik0j0(i)*(xik3j1(i) + xik3j3(i))**2 + 10*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma0a**2*sigma1a**2) + 20*xik2j2(i)*(3*xik0j0(i)*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) - 35*sigma1a**4) - 5*xik2j0(i)*(105*xik2j2(i)**2 + 6*xik0j0(i)*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + 35*sigma1a**4)) + sigma0a**2*(-60*xik2j2(i)**2*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma0a**2 + 350*xik2j0(i)**3*sigma1a**2 + 1400*xik2j2(i)**3*sigma1a**2 + 15*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma0a**2*sigma1a**4 - 15*xik2j0(i)**2*(7*xik4j0(i)*sigma0a**2 + 10*xik4j2(i)*sigma0a**2 + 3*xik4j4(i)*sigma0a**2 + 70*xik2j2(i)*sigma1a**2) + xik2j0(i)*(-189*xik3j1(i)**2*sigma0a**2 - 378*xik3j1(i)*xik3j3(i)*sigma0a**2 - 189*xik3j3(i)**2*sigma0a**2 + 420*xik2j2(i)*xik4j0(i)*sigma0a**2 + 600*xik2j2(i)*xik4j2(i)*sigma0a**2 + 180*xik2j2(i)*xik4j4(i)*sigma0a**2 - 350*sigma1a**6) + 14*xik2j2(i)*(27*xik3j1(i)**2*sigma0a**2 + 54*xik3j1(i)*xik3j3(i)*sigma0a**2 + 27*xik3j3(i)**2*sigma0a**2 - 25*sigma1a**6)))/(1575.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(1,9)=(-945*xik1j1(i)**3*(4*xik2j0(i)*xik3j1(i) + xik2j2(i)*xik3j1(i) - xik2j0(i)*xik3j3(i) - 4*xik2j2(i)*xik3j3(i)) + 135*xik1j1(i)**4*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + 3*xik0j0(i)**2*(-126*xik2j2(i)*(3*xik3j1(i)**2 + xik3j1(i)*xik3j3(i) - 2*xik3j3(i)**2) + 5*xik2j0(i)**2*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + 20*xik2j2(i)**2*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + xik2j0(i)*(189*xik3j1(i)**2 + 63*xik3j1(i)*xik3j3(i) + 2*(-63*xik3j3(i)**2 + 10*xik2j2(i)*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i)))) + 5*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i))*sigma1a**4) - 175*xik0j0(i)*(xik2j0(i)**4 - 5*xik2j0(i)**3*xik2j2(i) + 6*xik2j0(i)**2*xik2j2(i)**2 + 4*xik2j0(i)*xik2j2(i)**3 - 8*xik2j2(i)**4 - 3*xik2j0(i)*xik2j2(i)*sigma1a**4 + 6*xik2j2(i)**2*sigma1a**4 - sigma1a**8) - 3*xik1j1(i)**2*(175*xik2j0(i)**3 - 525*xik2j0(i)**2*xik2j2(i) + 700*xik2j2(i)**3 - 3*(63*xik0j0(i)*(3*xik3j1(i)**2 + xik3j1(i)*xik3j3(i) - 2*xik3j3(i)**2) + 10*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i))*sigma0a**2*sigma1a**2) + 5*xik2j2(i)*(12*xik0j0(i)*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) - 35*sigma1a**4) - 5*xik2j0(i)*(6*xik0j0(i)*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + 35*sigma1a**4)) - 315*xik1j1(i)*((-8*xik2j0(i)*xik3j1(i) + 7*xik2j2(i)*xik3j1(i) + 2*xik2j0(i)*xik3j3(i) + 2*xik2j2(i)*xik3j3(i))*sigma0a**2*sigma1a**2 - xik0j0(i)*(2*xik2j2(i)**2*(xik3j1(i) - 4*xik3j3(i)) + xik2j0(i)**2*(-4*xik3j1(i) + xik3j3(i)) + xik2j0(i)*xik2j2(i)*(7*xik3j1(i) + 2*xik3j3(i)) + (-4*xik3j1(i) + xik3j3(i))*sigma1a**4)) + sigma0a**2*(-60*xik2j2(i)**2*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i))*sigma0a**2 + 350*xik2j0(i)**3*sigma1a**2 - 700*xik2j2(i)**3*sigma1a**2 + 15*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i))*sigma0a**2*sigma1a**4 - 15*xik2j0(i)**2*(7*xik4j0(i)*sigma0a**2 - 5*xik4j2(i)*sigma0a**2 - 12*xik4j4(i)*sigma0a**2 + 105*xik2j2(i)*sigma1a**2) + 7*xik2j2(i)*(162*xik3j1(i)**2*sigma0a**2 + 54*xik3j1(i)*xik3j3(i)*sigma0a**2 - 108*xik3j3(i)**2*sigma0a**2 + 25*sigma1a**6) - xik2j0(i)*(567*xik3j1(i)**2*sigma0a**2 + 189*xik3j1(i)*xik3j3(i)*sigma0a**2 - 378*xik3j3(i)**2*sigma0a**2 + 10*(6*xik2j2(i)*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i))*sigma0a**2 - 210*xik2j2(i)**2*sigma1a**2 + 35*sigma1a**6))))/(1575.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(1,10)=0.0d0
vinvmat(1,11)=0.0d0
vinvmat(1,12)=0.0d0
vinvmat(2,1)=(25*xik2j0(i)**4*sigma0a**2 + 100*xik2j2(i)**4*sigma0a**2 + 81*xik1j1(i)**2*xik3j1(i)**2*sigma0a**2 + 162*xik1j1(i)**2*xik3j1(i)*xik3j3(i)*sigma0a**2 + 81*xik1j1(i)**2*xik3j3(i)**2*sigma0a**2 - 200*xik0j0(i)*xik2j2(i)**3*sigma1a**2 - 270*xik1j1(i)**3*xik3j1(i)*sigma1a**2 + 27*xik0j0(i)**2*xik3j1(i)**2*sigma1a**2 - 270*xik1j1(i)**3*xik3j3(i)*sigma1a**2 + 54*xik0j0(i)**2*xik3j1(i)*xik3j3(i)*sigma1a**2 + 27*xik0j0(i)**2*xik3j3(i)**2*sigma1a**2 - 27*xik3j1(i)**2*sigma0a**4*sigma1a**2 - 54*xik3j1(i)*xik3j3(i)*sigma0a**4*sigma1a**2 - 27*xik3j3(i)**2*sigma0a**4*sigma1a**2 + 90*xik1j1(i)*xik3j1(i)*sigma0a**2*sigma1a**4 + 90*xik1j1(i)*xik3j3(i)*sigma0a**2*sigma1a**4 + 75*xik1j1(i)**2*sigma1a**6 - 25*sigma0a**2*sigma1a**8 - 50*xik2j0(i)**3*(xik2j2(i)*sigma0a**2 + xik0j0(i)*sigma1a**2) + 135*xik1j1(i)**4*sigma2a**2 - 90*xik1j1(i)**2*sigma0a**2*sigma1a**2*sigma2a**2 - 15*xik0j0(i)**2*sigma1a**4*sigma2a**2 + 15*sigma0a**4*sigma1a**4*sigma2a**2 + 10*xik0j0(i)*xik2j2(i)*(9*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma1a**2 + 5*sigma1a**6 - 18*xik1j1(i)**2*sigma2a**2) + 15*xik2j2(i)**2*(-12*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 25*xik1j1(i)**2*sigma1a**2 + 5*sigma0a**2*sigma1a**4 + 4*xik0j0(i)**2*sigma2a**2 - 4*sigma0a**4*sigma2a**2) - 15*xik2j0(i)**2*(5*xik2j2(i)**2*sigma0a**2 - 6*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 5*xik1j1(i)**2*sigma1a**2 - 10*xik0j0(i)*xik2j2(i)*sigma1a**2 + (-xik0j0(i)**2 + sigma0a**4)*sigma2a**2) + 10*xik2j0(i)*(10*xik2j2(i)**3*sigma0a**2 + xik0j0(i)*(-18*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 9*xik1j1(i)**2*sigma2a**2) + 3*xik2j2(i)*(-3*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 10*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4 - 2*xik0j0(i)**2*sigma2a**2 + 2*sigma0a**4*sigma2a**2)))/(225.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(2,2)=(25*xik2j0(i)**4*sigma0a**2 + 100*xik2j2(i)**4*sigma0a**2 + 81*xik1j1(i)**2*xik3j1(i)**2*sigma0a**2 + 162*xik1j1(i)**2*xik3j1(i)*xik3j3(i)*sigma0a**2 + 81*xik1j1(i)**2*xik3j3(i)**2*sigma0a**2 - 200*xik0j0(i)*xik2j2(i)**3*sigma1a**2 - 270*xik1j1(i)**3*xik3j1(i)*sigma1a**2 + 27*xik0j0(i)**2*xik3j1(i)**2*sigma1a**2 - 270*xik1j1(i)**3*xik3j3(i)*sigma1a**2 + 54*xik0j0(i)**2*xik3j1(i)*xik3j3(i)*sigma1a**2 + 27*xik0j0(i)**2*xik3j3(i)**2*sigma1a**2 - 27*xik3j1(i)**2*sigma0a**4*sigma1a**2 - 54*xik3j1(i)*xik3j3(i)*sigma0a**4*sigma1a**2 - 27*xik3j3(i)**2*sigma0a**4*sigma1a**2 + 90*xik1j1(i)*xik3j1(i)*sigma0a**2*sigma1a**4 + 90*xik1j1(i)*xik3j3(i)*sigma0a**2*sigma1a**4 + 75*xik1j1(i)**2*sigma1a**6 - 25*sigma0a**2*sigma1a**8 - 50*xik2j0(i)**3*(xik2j2(i)*sigma0a**2 + xik0j0(i)*sigma1a**2) + 405*xik1j1(i)**4*sigma2a**2 - 270*xik1j1(i)**2*sigma0a**2*sigma1a**2*sigma2a**2 - 45*xik0j0(i)**2*sigma1a**4*sigma2a**2 + 45*sigma0a**4*sigma1a**4*sigma2a**2 + 10*xik0j0(i)*xik2j2(i)*(9*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma1a**2 + 5*sigma1a**6 - 54*xik1j1(i)**2*sigma2a**2) + 15*xik2j2(i)**2*(-12*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 25*xik1j1(i)**2*sigma1a**2 + 5*sigma0a**2*sigma1a**4 + 12*xik0j0(i)**2*sigma2a**2 - 12*sigma0a**4*sigma2a**2) - 15*xik2j0(i)**2*(5*xik2j2(i)**2*sigma0a**2 - 6*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 5*xik1j1(i)**2*sigma1a**2 - 10*xik0j0(i)*xik2j2(i)*sigma1a**2 + 3*(-xik0j0(i)**2 + sigma0a**4)*sigma2a**2) + 10*xik2j0(i)*(10*xik2j2(i)**3*sigma0a**2 + xik0j0(i)*(-18*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 27*xik1j1(i)**2*sigma2a**2) + 3*xik2j2(i)*(-3*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 10*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4 - 6*xik0j0(i)**2*sigma2a**2 + 6*sigma0a**4*sigma2a**2)))/(225.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(2,3)=(25*xik2j0(i)**4*sigma0a**2 - 200*xik2j2(i)**4*sigma0a**2 + 243*xik1j1(i)**2*xik3j1(i)**2*sigma0a**2 + 81*xik1j1(i)**2*xik3j1(i)*xik3j3(i)*sigma0a**2 - 162*xik1j1(i)**2*xik3j3(i)**2*sigma0a**2 + 100*xik0j0(i)*xik2j2(i)**3*sigma1a**2 - 540*xik1j1(i)**3*xik3j1(i)*sigma1a**2 + 81*xik0j0(i)**2*xik3j1(i)**2*sigma1a**2 + 135*xik1j1(i)**3*xik3j3(i)*sigma1a**2 + 27*xik0j0(i)**2*xik3j1(i)*xik3j3(i)*sigma1a**2 - 54*xik0j0(i)**2*xik3j3(i)**2*sigma1a**2 - 81*xik3j1(i)**2*sigma0a**4*sigma1a**2 - 27*xik3j1(i)*xik3j3(i)*sigma0a**4*sigma1a**2 + 54*xik3j3(i)**2*sigma0a**4*sigma1a**2 + 180*xik1j1(i)*xik3j1(i)*sigma0a**2*sigma1a**4 - 45*xik1j1(i)*xik3j3(i)*sigma0a**2*sigma1a**4 + 75*xik1j1(i)**2*sigma1a**6 - 25*sigma0a**2*sigma1a**8 - 25*xik2j0(i)**3*(5*xik2j2(i)*sigma0a**2 + 2*xik0j0(i)*sigma1a**2) + 135*xik1j1(i)**4*sigma2a**2 - 90*xik1j1(i)**2*sigma0a**2*sigma1a**2*sigma2a**2 - 15*xik0j0(i)**2*sigma1a**4*sigma2a**2 + 15*sigma0a**4*sigma1a**4*sigma2a**2 - 5*xik0j0(i)*xik2j2(i)*(-9*xik1j1(i)*(7*xik3j1(i) + 2*xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 36*xik1j1(i)**2*sigma2a**2) - 30*xik2j2(i)**2*(3*xik1j1(i)*(xik3j1(i) - 4*xik3j3(i))*sigma0a**2 + 10*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4 - 2*xik0j0(i)**2*sigma2a**2 + 2*sigma0a**4*sigma2a**2) + 15*xik2j0(i)**2*(10*xik2j2(i)**2*sigma0a**2 + 3*xik1j1(i)*(4*xik3j1(i) - xik3j3(i))*sigma0a**2 - 5*xik1j1(i)**2*sigma1a**2 + 15*xik0j0(i)*xik2j2(i)*sigma1a**2 + (xik0j0(i)**2 - sigma0a**4)*sigma2a**2) + 5*xik2j0(i)*(20*xik2j2(i)**3*sigma0a**2 - 60*xik0j0(i)*xik2j2(i)**2*sigma1a**2 + 2*xik0j0(i)*(9*xik1j1(i)*(-4*xik3j1(i) + xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 9*xik1j1(i)**2*sigma2a**2) + 3*xik2j2(i)*(-3*xik1j1(i)*(7*xik3j1(i) + 2*xik3j3(i))*sigma0a**2 + 20*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4 - 4*xik0j0(i)**2*sigma2a**2 + 4*sigma0a**4*sigma2a**2)))/(225.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(2,4)=0.0d0
vinvmat(2,5)=0.0d0
vinvmat(2,6)=0.0d0
vinvmat(2,7)=(-1890*xik1j1(i)**3*(xik2j0(i) + xik2j2(i))*(xik3j1(i) + xik3j3(i)) + 135*xik1j1(i)**4*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + 3*xik0j0(i)**2*(-126*xik2j2(i)*(xik3j1(i) + xik3j3(i))**2 + 5*xik2j0(i)**2*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + 20*xik2j2(i)**2*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + xik2j0(i)*(63*xik3j1(i)**2 + 126*xik3j1(i)*xik3j3(i) + 63*xik3j3(i)**2 - 20*xik2j2(i)*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))) - 5*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma1a**4) - 175*xik0j0(i)*(xik2j0(i)**4 - 2*xik2j0(i)**3*xik2j2(i) - 3*xik2j0(i)**2*xik2j2(i)**2 + 4*xik2j0(i)*xik2j2(i)**3 + 4*xik2j2(i)**4 - 6*xik2j0(i)*xik2j2(i)*sigma1a**4 + 3*xik2j2(i)**2*sigma1a**4 - sigma1a**8) - 630*xik1j1(i)*(xik3j1(i) + xik3j3(i))*((-2*xik2j0(i) + xik2j2(i))*sigma0a**2*sigma1a**2 + xik0j0(i)*(xik2j0(i)**2 - xik2j0(i)*xik2j2(i) - 2*xik2j2(i)**2 + sigma1a**4)) - 3*xik1j1(i)**2*(175*xik2j0(i)**3 - 350*xik2j2(i)**3 + 3*(-63*xik0j0(i)*(xik3j1(i) + xik3j3(i))**2 + 10*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma0a**2*sigma1a**2) + 20*xik2j2(i)*(3*xik0j0(i)*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) - 35*sigma1a**4) - 5*xik2j0(i)*(105*xik2j2(i)**2 + 6*xik0j0(i)*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + 35*sigma1a**4)) + sigma0a**2*(-60*xik2j2(i)**2*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma0a**2 + 350*xik2j0(i)**3*sigma1a**2 + 1400*xik2j2(i)**3*sigma1a**2 + 15*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma0a**2*sigma1a**4 - 15*xik2j0(i)**2*(7*xik4j0(i)*sigma0a**2 + 10*xik4j2(i)*sigma0a**2 + 3*xik4j4(i)*sigma0a**2 + 70*xik2j2(i)*sigma1a**2) + xik2j0(i)*(-189*xik3j1(i)**2*sigma0a**2 - 378*xik3j1(i)*xik3j3(i)*sigma0a**2 - 189*xik3j3(i)**2*sigma0a**2 + 420*xik2j2(i)*xik4j0(i)*sigma0a**2 + 600*xik2j2(i)*xik4j2(i)*sigma0a**2 + 180*xik2j2(i)*xik4j4(i)*sigma0a**2 - 350*sigma1a**6) + 14*xik2j2(i)*(27*xik3j1(i)**2*sigma0a**2 + 54*xik3j1(i)*xik3j3(i)*sigma0a**2 + 27*xik3j3(i)**2*sigma0a**2 - 25*sigma1a**6)))/(1575.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(2,8)=(-1890*xik1j1(i)**3*(xik2j0(i) + xik2j2(i))*(xik3j1(i) + xik3j3(i)) + 405*xik1j1(i)**4*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + 9*xik0j0(i)**2*(-42*xik2j2(i)*(xik3j1(i) + xik3j3(i))**2 + 5*xik2j0(i)**2*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + 20*xik2j2(i)**2*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + xik2j0(i)*(21*xik3j1(i)**2 + 42*xik3j1(i)*xik3j3(i) + 21*xik3j3(i)**2 - 20*xik2j2(i)*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))) - 5*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma1a**4) - 175*xik0j0(i)*(xik2j0(i)**4 - 2*xik2j0(i)**3*xik2j2(i) - 3*xik2j0(i)**2*xik2j2(i)**2 + 4*xik2j0(i)*xik2j2(i)**3 + 4*xik2j2(i)**4 - 6*xik2j0(i)*xik2j2(i)*sigma1a**4 + 3*xik2j2(i)**2*sigma1a**4 - sigma1a**8) - 630*xik1j1(i)*(xik3j1(i) + xik3j3(i))*((-2*xik2j0(i) + xik2j2(i))*sigma0a**2*sigma1a**2 + xik0j0(i)*(xik2j0(i)**2 - xik2j0(i)*xik2j2(i) - 2*xik2j2(i)**2 + sigma1a**4)) - 3*xik1j1(i)**2*(175*xik2j0(i)**3 - 350*xik2j2(i)**3 - 9*(21*xik0j0(i)*(xik3j1(i) + xik3j3(i))**2 - 10*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma0a**2*sigma1a**2) + 20*xik2j2(i)*(9*xik0j0(i)*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) - 35*sigma1a**4) - 5*xik2j0(i)*(105*xik2j2(i)**2 + 18*xik0j0(i)*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + 35*sigma1a**4)) + sigma0a**2*(-180*xik2j2(i)**2*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma0a**2 + 350*xik2j0(i)**3*sigma1a**2 + 1400*xik2j2(i)**3*sigma1a**2 + 45*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma0a**2*sigma1a**4 - 15*xik2j0(i)**2*(21*xik4j0(i)*sigma0a**2 + 30*xik4j2(i)*sigma0a**2 + 9*xik4j4(i)*sigma0a**2 + 70*xik2j2(i)*sigma1a**2) + xik2j0(i)*(-189*xik3j1(i)**2*sigma0a**2 - 378*xik3j1(i)*xik3j3(i)*sigma0a**2 - 189*xik3j3(i)**2*sigma0a**2 + 1260*xik2j2(i)*xik4j0(i)*sigma0a**2 + 1800*xik2j2(i)*xik4j2(i)*sigma0a**2 + 540*xik2j2(i)*xik4j4(i)*sigma0a**2 - 350*sigma1a**6) + 14*xik2j2(i)*(27*xik3j1(i)**2*sigma0a**2 + 54*xik3j1(i)*xik3j3(i)*sigma0a**2 + 27*xik3j3(i)**2*sigma0a**2 - 25*sigma1a**6)))/(1575.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(2,9)=(-945*xik1j1(i)**3*(4*xik2j0(i)*xik3j1(i) + xik2j2(i)*xik3j1(i) - xik2j0(i)*xik3j3(i) - 4*xik2j2(i)*xik3j3(i)) + 135*xik1j1(i)**4*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + 3*xik0j0(i)**2*(-126*xik2j2(i)*(3*xik3j1(i)**2 + xik3j1(i)*xik3j3(i) - 2*xik3j3(i)**2) + 5*xik2j0(i)**2*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + 20*xik2j2(i)**2*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + xik2j0(i)*(189*xik3j1(i)**2 + 63*xik3j1(i)*xik3j3(i) + 2*(-63*xik3j3(i)**2 + 10*xik2j2(i)*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i)))) + 5*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i))*sigma1a**4) - 175*xik0j0(i)*(xik2j0(i)**4 - 5*xik2j0(i)**3*xik2j2(i) + 6*xik2j0(i)**2*xik2j2(i)**2 + 4*xik2j0(i)*xik2j2(i)**3 - 8*xik2j2(i)**4 - 3*xik2j0(i)*xik2j2(i)*sigma1a**4 + 6*xik2j2(i)**2*sigma1a**4 - sigma1a**8) - 3*xik1j1(i)**2*(175*xik2j0(i)**3 - 525*xik2j0(i)**2*xik2j2(i) + 700*xik2j2(i)**3 - 3*(63*xik0j0(i)*(3*xik3j1(i)**2 + xik3j1(i)*xik3j3(i) - 2*xik3j3(i)**2) + 10*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i))*sigma0a**2*sigma1a**2) + 5*xik2j2(i)*(12*xik0j0(i)*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) - 35*sigma1a**4) - 5*xik2j0(i)*(6*xik0j0(i)*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + 35*sigma1a**4)) - 315*xik1j1(i)*((-8*xik2j0(i)*xik3j1(i) + 7*xik2j2(i)*xik3j1(i) + 2*xik2j0(i)*xik3j3(i) + 2*xik2j2(i)*xik3j3(i))*sigma0a**2*sigma1a**2 - xik0j0(i)*(2*xik2j2(i)**2*(xik3j1(i) - 4*xik3j3(i)) + xik2j0(i)**2*(-4*xik3j1(i) + xik3j3(i)) + xik2j0(i)*xik2j2(i)*(7*xik3j1(i) + 2*xik3j3(i)) + (-4*xik3j1(i) + xik3j3(i))*sigma1a**4)) + sigma0a**2*(-60*xik2j2(i)**2*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i))*sigma0a**2 + 350*xik2j0(i)**3*sigma1a**2 - 700*xik2j2(i)**3*sigma1a**2 + 15*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i))*sigma0a**2*sigma1a**4 - 15*xik2j0(i)**2*(7*xik4j0(i)*sigma0a**2 - 5*xik4j2(i)*sigma0a**2 - 12*xik4j4(i)*sigma0a**2 + 105*xik2j2(i)*sigma1a**2) + 7*xik2j2(i)*(162*xik3j1(i)**2*sigma0a**2 + 54*xik3j1(i)*xik3j3(i)*sigma0a**2 - 108*xik3j3(i)**2*sigma0a**2 + 25*sigma1a**6) - xik2j0(i)*(567*xik3j1(i)**2*sigma0a**2 + 189*xik3j1(i)*xik3j3(i)*sigma0a**2 - 378*xik3j3(i)**2*sigma0a**2 + 10*(6*xik2j2(i)*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i))*sigma0a**2 - 210*xik2j2(i)**2*sigma1a**2 + 35*sigma1a**6))))/(1575.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(2,10)=0.0d0
vinvmat(2,11)=0.0d0
vinvmat(2,12)=0.0d0
vinvmat(3,1)=(25*xik2j0(i)**4*sigma0a**2 - 200*xik2j2(i)**4*sigma0a**2 + 243*xik1j1(i)**2*xik3j1(i)**2*sigma0a**2 + 81*xik1j1(i)**2*xik3j1(i)*xik3j3(i)*sigma0a**2 - 162*xik1j1(i)**2*xik3j3(i)**2*sigma0a**2 + 100*xik0j0(i)*xik2j2(i)**3*sigma1a**2 - 540*xik1j1(i)**3*xik3j1(i)*sigma1a**2 + 81*xik0j0(i)**2*xik3j1(i)**2*sigma1a**2 + 135*xik1j1(i)**3*xik3j3(i)*sigma1a**2 + 27*xik0j0(i)**2*xik3j1(i)*xik3j3(i)*sigma1a**2 - 54*xik0j0(i)**2*xik3j3(i)**2*sigma1a**2 - 81*xik3j1(i)**2*sigma0a**4*sigma1a**2 - 27*xik3j1(i)*xik3j3(i)*sigma0a**4*sigma1a**2 + 54*xik3j3(i)**2*sigma0a**4*sigma1a**2 + 180*xik1j1(i)*xik3j1(i)*sigma0a**2*sigma1a**4 - 45*xik1j1(i)*xik3j3(i)*sigma0a**2*sigma1a**4 + 75*xik1j1(i)**2*sigma1a**6 - 25*sigma0a**2*sigma1a**8 - 25*xik2j0(i)**3*(5*xik2j2(i)*sigma0a**2 + 2*xik0j0(i)*sigma1a**2) + 135*xik1j1(i)**4*sigma2a**2 - 90*xik1j1(i)**2*sigma0a**2*sigma1a**2*sigma2a**2 - 15*xik0j0(i)**2*sigma1a**4*sigma2a**2 + 15*sigma0a**4*sigma1a**4*sigma2a**2 - 5*xik0j0(i)*xik2j2(i)*(-9*xik1j1(i)*(7*xik3j1(i) + 2*xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 36*xik1j1(i)**2*sigma2a**2) - 30*xik2j2(i)**2*(3*xik1j1(i)*(xik3j1(i) - 4*xik3j3(i))*sigma0a**2 + 10*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4 - 2*xik0j0(i)**2*sigma2a**2 + 2*sigma0a**4*sigma2a**2) + 15*xik2j0(i)**2*(10*xik2j2(i)**2*sigma0a**2 + 3*xik1j1(i)*(4*xik3j1(i) - xik3j3(i))*sigma0a**2 - 5*xik1j1(i)**2*sigma1a**2 + 15*xik0j0(i)*xik2j2(i)*sigma1a**2 + (xik0j0(i)**2 - sigma0a**4)*sigma2a**2) + 5*xik2j0(i)*(20*xik2j2(i)**3*sigma0a**2 - 60*xik0j0(i)*xik2j2(i)**2*sigma1a**2 + 2*xik0j0(i)*(9*xik1j1(i)*(-4*xik3j1(i) + xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 9*xik1j1(i)**2*sigma2a**2) + 3*xik2j2(i)*(-3*xik1j1(i)*(7*xik3j1(i) + 2*xik3j3(i))*sigma0a**2 + 20*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4 - 4*xik0j0(i)**2*sigma2a**2 + 4*sigma0a**4*sigma2a**2)))/(225.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(3,2)=(25*xik2j0(i)**4*sigma0a**2 - 200*xik2j2(i)**4*sigma0a**2 + 243*xik1j1(i)**2*xik3j1(i)**2*sigma0a**2 + 81*xik1j1(i)**2*xik3j1(i)*xik3j3(i)*sigma0a**2 - 162*xik1j1(i)**2*xik3j3(i)**2*sigma0a**2 + 100*xik0j0(i)*xik2j2(i)**3*sigma1a**2 - 540*xik1j1(i)**3*xik3j1(i)*sigma1a**2 + 81*xik0j0(i)**2*xik3j1(i)**2*sigma1a**2 + 135*xik1j1(i)**3*xik3j3(i)*sigma1a**2 + 27*xik0j0(i)**2*xik3j1(i)*xik3j3(i)*sigma1a**2 - 54*xik0j0(i)**2*xik3j3(i)**2*sigma1a**2 - 81*xik3j1(i)**2*sigma0a**4*sigma1a**2 - 27*xik3j1(i)*xik3j3(i)*sigma0a**4*sigma1a**2 + 54*xik3j3(i)**2*sigma0a**4*sigma1a**2 + 180*xik1j1(i)*xik3j1(i)*sigma0a**2*sigma1a**4 - 45*xik1j1(i)*xik3j3(i)*sigma0a**2*sigma1a**4 + 75*xik1j1(i)**2*sigma1a**6 - 25*sigma0a**2*sigma1a**8 - 25*xik2j0(i)**3*(5*xik2j2(i)*sigma0a**2 + 2*xik0j0(i)*sigma1a**2) + 135*xik1j1(i)**4*sigma2a**2 - 90*xik1j1(i)**2*sigma0a**2*sigma1a**2*sigma2a**2 - 15*xik0j0(i)**2*sigma1a**4*sigma2a**2 + 15*sigma0a**4*sigma1a**4*sigma2a**2 - 5*xik0j0(i)*xik2j2(i)*(-9*xik1j1(i)*(7*xik3j1(i) + 2*xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 36*xik1j1(i)**2*sigma2a**2) - 30*xik2j2(i)**2*(3*xik1j1(i)*(xik3j1(i) - 4*xik3j3(i))*sigma0a**2 + 10*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4 - 2*xik0j0(i)**2*sigma2a**2 + 2*sigma0a**4*sigma2a**2) + 15*xik2j0(i)**2*(10*xik2j2(i)**2*sigma0a**2 + 3*xik1j1(i)*(4*xik3j1(i) - xik3j3(i))*sigma0a**2 - 5*xik1j1(i)**2*sigma1a**2 + 15*xik0j0(i)*xik2j2(i)*sigma1a**2 + (xik0j0(i)**2 - sigma0a**4)*sigma2a**2) + 5*xik2j0(i)*(20*xik2j2(i)**3*sigma0a**2 - 60*xik0j0(i)*xik2j2(i)**2*sigma1a**2 + 2*xik0j0(i)*(9*xik1j1(i)*(-4*xik3j1(i) + xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 9*xik1j1(i)**2*sigma2a**2) + 3*xik2j2(i)*(-3*xik1j1(i)*(7*xik3j1(i) + 2*xik3j3(i))*sigma0a**2 + 20*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4 - 4*xik0j0(i)**2*sigma2a**2 + 4*sigma0a**4*sigma2a**2)))/(225.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(3,3)=(25*xik2j0(i)**4*sigma0a**2 + 400*xik2j2(i)**4*sigma0a**2 + 729*xik1j1(i)**2*xik3j1(i)**2*sigma0a**2 - 972*xik1j1(i)**2*xik3j1(i)*xik3j3(i)*sigma0a**2 + 324*xik1j1(i)**2*xik3j3(i)**2*sigma0a**2 + 400*xik0j0(i)*xik2j2(i)**3*sigma1a**2 - 810*xik1j1(i)**3*xik3j1(i)*sigma1a**2 + 243*xik0j0(i)**2*xik3j1(i)**2*sigma1a**2 + 540*xik1j1(i)**3*xik3j3(i)*sigma1a**2 - 324*xik0j0(i)**2*xik3j1(i)*xik3j3(i)*sigma1a**2 + 108*xik0j0(i)**2*xik3j3(i)**2*sigma1a**2 - 243*xik3j1(i)**2*sigma0a**4*sigma1a**2 + 324*xik3j1(i)*xik3j3(i)*sigma0a**4*sigma1a**2 - 108*xik3j3(i)**2*sigma0a**4*sigma1a**2 + 270*xik1j1(i)*xik3j1(i)*sigma0a**2*sigma1a**4 - 180*xik1j1(i)*xik3j3(i)*sigma0a**2*sigma1a**4 + 75*xik1j1(i)**2*sigma1a**6 - 25*sigma0a**2*sigma1a**8 - 50*xik2j0(i)**3*(4*xik2j2(i)*sigma0a**2 + xik0j0(i)*sigma1a**2) + 405*xik1j1(i)**4*sigma2a**2 - 270*xik1j1(i)**2*sigma0a**2*sigma1a**2*sigma2a**2 - 45*xik0j0(i)**2*sigma1a**4*sigma2a**2 + 45*sigma0a**4*sigma1a**4*sigma2a**2 - 20*xik0j0(i)*xik2j2(i)*(18*xik1j1(i)*(-3*xik3j1(i) + 2*xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 27*xik1j1(i)**2*sigma2a**2) + 15*xik2j0(i)**2*(40*xik2j2(i)**2*sigma0a**2 + 6*xik1j1(i)*(3*xik3j1(i) - 2*xik3j3(i))*sigma0a**2 - 5*xik1j1(i)**2*sigma1a**2 + 20*xik0j0(i)*xik2j2(i)*sigma1a**2 + 3*(xik0j0(i)**2 - sigma0a**4)*sigma2a**2) - 60*xik2j2(i)**2*(-6*xik1j1(i)*(3*xik3j1(i) - 2*xik3j3(i))*sigma0a**2 + 5*xik1j1(i)**2*sigma1a**2 + 3*(-xik0j0(i)**2 + sigma0a**4)*sigma2a**2) - 10*xik2j0(i)*(80*xik2j2(i)**3*sigma0a**2 + 60*xik0j0(i)*xik2j2(i)**2*sigma1a**2 - xik0j0(i)*(18*xik1j1(i)*(-3*xik3j1(i) + 2*xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 27*xik1j1(i)**2*sigma2a**2) - 6*xik2j2(i)*(-6*xik1j1(i)*(3*xik3j1(i) - 2*xik3j3(i))*sigma0a**2 + 5*xik1j1(i)**2*sigma1a**2 + 3*(-xik0j0(i)**2 + sigma0a**4)*sigma2a**2)))/(225.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(3,4)=0.0d0
vinvmat(3,5)=0.0d0
vinvmat(3,6)=0.0d0
vinvmat(3,7)=(-945*xik1j1(i)**3*(4*xik2j0(i)*xik3j1(i) + xik2j2(i)*xik3j1(i) - xik2j0(i)*xik3j3(i) - 4*xik2j2(i)*xik3j3(i)) + 135*xik1j1(i)**4*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + 3*xik0j0(i)**2*(-126*xik2j2(i)*(3*xik3j1(i)**2 + xik3j1(i)*xik3j3(i) - 2*xik3j3(i)**2) + 5*xik2j0(i)**2*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + 20*xik2j2(i)**2*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + xik2j0(i)*(189*xik3j1(i)**2 + 63*xik3j1(i)*xik3j3(i) + 2*(-63*xik3j3(i)**2 + 10*xik2j2(i)*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i)))) + 5*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i))*sigma1a**4) - 175*xik0j0(i)*(xik2j0(i)**4 - 5*xik2j0(i)**3*xik2j2(i) + 6*xik2j0(i)**2*xik2j2(i)**2 + 4*xik2j0(i)*xik2j2(i)**3 - 8*xik2j2(i)**4 - 3*xik2j0(i)*xik2j2(i)*sigma1a**4 + 6*xik2j2(i)**2*sigma1a**4 - sigma1a**8) - 3*xik1j1(i)**2*(175*xik2j0(i)**3 - 525*xik2j0(i)**2*xik2j2(i) + 700*xik2j2(i)**3 - 3*(63*xik0j0(i)*(3*xik3j1(i)**2 + xik3j1(i)*xik3j3(i) - 2*xik3j3(i)**2) + 10*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i))*sigma0a**2*sigma1a**2) + 5*xik2j2(i)*(12*xik0j0(i)*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) - 35*sigma1a**4) - 5*xik2j0(i)*(6*xik0j0(i)*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + 35*sigma1a**4)) - 315*xik1j1(i)*((-8*xik2j0(i)*xik3j1(i) + 7*xik2j2(i)*xik3j1(i) + 2*xik2j0(i)*xik3j3(i) + 2*xik2j2(i)*xik3j3(i))*sigma0a**2*sigma1a**2 - xik0j0(i)*(2*xik2j2(i)**2*(xik3j1(i) - 4*xik3j3(i)) + xik2j0(i)**2*(-4*xik3j1(i) + xik3j3(i)) + xik2j0(i)*xik2j2(i)*(7*xik3j1(i) + 2*xik3j3(i)) + (-4*xik3j1(i) + xik3j3(i))*sigma1a**4)) + sigma0a**2*(-60*xik2j2(i)**2*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i))*sigma0a**2 + 350*xik2j0(i)**3*sigma1a**2 - 700*xik2j2(i)**3*sigma1a**2 + 15*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i))*sigma0a**2*sigma1a**4 - 15*xik2j0(i)**2*(7*xik4j0(i)*sigma0a**2 - 5*xik4j2(i)*sigma0a**2 - 12*xik4j4(i)*sigma0a**2 + 105*xik2j2(i)*sigma1a**2) + 7*xik2j2(i)*(162*xik3j1(i)**2*sigma0a**2 + 54*xik3j1(i)*xik3j3(i)*sigma0a**2 - 108*xik3j3(i)**2*sigma0a**2 + 25*sigma1a**6) - xik2j0(i)*(567*xik3j1(i)**2*sigma0a**2 + 189*xik3j1(i)*xik3j3(i)*sigma0a**2 - 378*xik3j3(i)**2*sigma0a**2 + 10*(6*xik2j2(i)*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i))*sigma0a**2 - 210*xik2j2(i)**2*sigma1a**2 + 35*sigma1a**6))))/(1575.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(3,8)=(-945*xik1j1(i)**3*(4*xik2j0(i)*xik3j1(i) + xik2j2(i)*xik3j1(i) - xik2j0(i)*xik3j3(i) - 4*xik2j2(i)*xik3j3(i)) + 135*xik1j1(i)**4*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + 3*xik0j0(i)**2*(-126*xik2j2(i)*(3*xik3j1(i)**2 + xik3j1(i)*xik3j3(i) - 2*xik3j3(i)**2) + 5*xik2j0(i)**2*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + 20*xik2j2(i)**2*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + xik2j0(i)*(189*xik3j1(i)**2 + 63*xik3j1(i)*xik3j3(i) + 2*(-63*xik3j3(i)**2 + 10*xik2j2(i)*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i)))) + 5*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i))*sigma1a**4) - 175*xik0j0(i)*(xik2j0(i)**4 - 5*xik2j0(i)**3*xik2j2(i) + 6*xik2j0(i)**2*xik2j2(i)**2 + 4*xik2j0(i)*xik2j2(i)**3 - 8*xik2j2(i)**4 - 3*xik2j0(i)*xik2j2(i)*sigma1a**4 + 6*xik2j2(i)**2*sigma1a**4 - sigma1a**8) - 3*xik1j1(i)**2*(175*xik2j0(i)**3 - 525*xik2j0(i)**2*xik2j2(i) + 700*xik2j2(i)**3 - 3*(63*xik0j0(i)*(3*xik3j1(i)**2 + xik3j1(i)*xik3j3(i) - 2*xik3j3(i)**2) + 10*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i))*sigma0a**2*sigma1a**2) + 5*xik2j2(i)*(12*xik0j0(i)*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) - 35*sigma1a**4) - 5*xik2j0(i)*(6*xik0j0(i)*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + 35*sigma1a**4)) - 315*xik1j1(i)*((-8*xik2j0(i)*xik3j1(i) + 7*xik2j2(i)*xik3j1(i) + 2*xik2j0(i)*xik3j3(i) + 2*xik2j2(i)*xik3j3(i))*sigma0a**2*sigma1a**2 - xik0j0(i)*(2*xik2j2(i)**2*(xik3j1(i) - 4*xik3j3(i)) + xik2j0(i)**2*(-4*xik3j1(i) + xik3j3(i)) + xik2j0(i)*xik2j2(i)*(7*xik3j1(i) + 2*xik3j3(i)) + (-4*xik3j1(i) + xik3j3(i))*sigma1a**4)) + sigma0a**2*(-60*xik2j2(i)**2*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i))*sigma0a**2 + 350*xik2j0(i)**3*sigma1a**2 - 700*xik2j2(i)**3*sigma1a**2 + 15*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i))*sigma0a**2*sigma1a**4 - 15*xik2j0(i)**2*(7*xik4j0(i)*sigma0a**2 - 5*xik4j2(i)*sigma0a**2 - 12*xik4j4(i)*sigma0a**2 + 105*xik2j2(i)*sigma1a**2) + 7*xik2j2(i)*(162*xik3j1(i)**2*sigma0a**2 + 54*xik3j1(i)*xik3j3(i)*sigma0a**2 - 108*xik3j3(i)**2*sigma0a**2 + 25*sigma1a**6) - xik2j0(i)*(567*xik3j1(i)**2*sigma0a**2 + 189*xik3j1(i)*xik3j3(i)*sigma0a**2 - 378*xik3j3(i)**2*sigma0a**2 + 10*(6*xik2j2(i)*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i))*sigma0a**2 - 210*xik2j2(i)**2*sigma1a**2 + 35*sigma1a**6))))/(1575.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(3,9)=(-1890*xik1j1(i)**3*(xik2j0(i) - 2*xik2j2(i))*(3*xik3j1(i) - 2*xik3j3(i)) + 405*xik1j1(i)**4*(7*xik4j0(i) - 20*xik4j2(i) + 8*xik4j4(i)) + 9*xik0j0(i)**2*(-42*xik2j2(i)*(3*xik3j1(i) - 2*xik3j3(i))**2 + 5*xik2j0(i)**2*(7*xik4j0(i) - 20*xik4j2(i) + 8*xik4j4(i)) + 20*xik2j2(i)**2*(7*xik4j0(i) - 20*xik4j2(i) + 8*xik4j4(i)) + xik2j0(i)*(189*xik3j1(i)**2 - 252*xik3j1(i)*xik3j3(i) + 4*(21*xik3j3(i)**2 - 5*xik2j2(i)*(7*xik4j0(i) - 20*xik4j2(i) + 8*xik4j4(i)))) - 5*(7*xik4j0(i) - 20*xik4j2(i) + 8*xik4j4(i))*sigma1a**4) - 175*xik0j0(i)*(xik2j0(i)**4 - 8*xik2j0(i)**3*xik2j2(i) + 24*xik2j0(i)**2*xik2j2(i)**2 - 32*xik2j0(i)*xik2j2(i)**3 + 16*xik2j2(i)**4 - sigma1a**8) - 630*xik1j1(i)*(3*xik3j1(i) - 2*xik3j3(i))*(-2*(xik2j0(i) - 2*xik2j2(i))*sigma0a**2*sigma1a**2 + xik0j0(i)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 + sigma1a**4)) - 3*xik1j1(i)**2*(175*xik2j0(i)**3 - 1050*xik2j0(i)**2*xik2j2(i) - 1400*xik2j2(i)**3 - 9*(21*xik0j0(i)*(3*xik3j1(i) - 2*xik3j3(i))**2 - 10*(7*xik4j0(i) - 20*xik4j2(i) + 8*xik4j4(i))*sigma0a**2*sigma1a**2) + 5*xik2j0(i)*(420*xik2j2(i)**2 - 18*xik0j0(i)*(7*xik4j0(i) - 20*xik4j2(i) + 8*xik4j4(i)) - 35*sigma1a**4) + 10*xik2j2(i)*(18*xik0j0(i)*(7*xik4j0(i) - 20*xik4j2(i) + 8*xik4j4(i)) + 35*sigma1a**4)) + sigma0a**2*(-180*xik2j2(i)**2*(7*xik4j0(i) - 20*xik4j2(i) + 8*xik4j4(i))*sigma0a**2 + 350*xik2j0(i)**3*sigma1a**2 - 2800*xik2j2(i)**3*sigma1a**2 + 45*(7*xik4j0(i) - 20*xik4j2(i) + 8*xik4j4(i))*sigma0a**2*sigma1a**4 - 15*xik2j0(i)**2*(21*xik4j0(i)*sigma0a**2 - 60*xik4j2(i)*sigma0a**2 + 24*xik4j4(i)*sigma0a**2 + 140*xik2j2(i)*sigma1a**2) + xik2j0(i)*(-1701*xik3j1(i)**2*sigma0a**2 + 2268*xik3j1(i)*xik3j3(i)*sigma0a**2 - 756*xik3j3(i)**2*sigma0a**2 + 1260*xik2j2(i)*xik4j0(i)*sigma0a**2 - 3600*xik2j2(i)*xik4j2(i)*sigma0a**2 + 1440*xik2j2(i)*xik4j4(i)*sigma0a**2 + 4200*xik2j2(i)**2*sigma1a**2 - 350*sigma1a**6) + 14*xik2j2(i)*(243*xik3j1(i)**2*sigma0a**2 - 324*xik3j1(i)*xik3j3(i)*sigma0a**2 + 108*xik3j3(i)**2*sigma0a**2 + 50*sigma1a**6)))/(1575.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(3,10)=0.0d0
vinvmat(3,11)=0.0d0
vinvmat(3,12)=0.0d0
vinvmat(4,1)=0.0d0
vinvmat(4,2)=0.0d0
vinvmat(4,3)=0.0d0
vinvmat(4,4)=sigma2a**2/15.
vinvmat(4,5)=0.0d0
vinvmat(4,6)=0.0d0
vinvmat(4,7)=0.0d0
vinvmat(4,8)=0.0d0
vinvmat(4,9)=0.0d0
vinvmat(4,10)=(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))/105.
vinvmat(4,11)=0.0d0
vinvmat(4,12)=0.0d0
vinvmat(5,1)=0.0d0
vinvmat(5,2)=0.0d0
vinvmat(5,3)=0.0d0
vinvmat(5,4)=0.0d0
vinvmat(5,5)=((9*(xik3j1(i) + xik3j3(i))**2*sigma1a**2)/(5.*(xik2j0(i)**2 + 2*xik2j0(i)*xik2j2(i) + xik2j2(i)**2 - sigma1a**4)) + sigma2a**2)/15.
vinvmat(5,6)=0.0d0
vinvmat(5,7)=0.0d0
vinvmat(5,8)=0.0d0
vinvmat(5,9)=0.0d0
vinvmat(5,10)=0.0d0
vinvmat(5,11)=(35*xik4j0(i) - 25*xik4j2(i) - 60*xik4j4(i) + (63*(xik2j0(i) + xik2j2(i))*(xik3j1(i) + xik3j3(i))**2)/(xik2j0(i)**2 + 2*xik2j0(i)*xik2j2(i) + xik2j2(i)**2 - sigma1a**4))/525.
vinvmat(5,12)=0.0d0
vinvmat(6,1)=0.0d0
vinvmat(6,2)=0.0d0
vinvmat(6,3)=0.0d0
vinvmat(6,4)=0.0d0
vinvmat(6,5)=0.0d0
vinvmat(6,6)=((9*(xik3j1(i) + xik3j3(i))**2*sigma1a**2)/(5.*(xik2j0(i)**2 + 2*xik2j0(i)*xik2j2(i) + xik2j2(i)**2 - sigma1a**4)) + sigma2a**2)/15.
vinvmat(6,7)=0.0d0
vinvmat(6,8)=0.0d0
vinvmat(6,9)=0.0d0
vinvmat(6,10)=0.0d0
vinvmat(6,11)=0.0d0
vinvmat(6,12)=(35*xik4j0(i) - 25*xik4j2(i) - 60*xik4j4(i) + (63*(xik2j0(i) + xik2j2(i))*(xik3j1(i) + xik3j3(i))**2)/(xik2j0(i)**2 + 2*xik2j0(i)*xik2j2(i) + xik2j2(i)**2 - sigma1a**4))/525.
vinvmat(7,1)=(-1890*xik1j1(i)**3*(xik2j0(i) + xik2j2(i))*(xik3j1(i) + xik3j3(i)) + 405*xik1j1(i)**4*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + 9*xik0j0(i)**2*(-42*xik2j2(i)*(xik3j1(i) + xik3j3(i))**2 + 5*xik2j0(i)**2*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + 20*xik2j2(i)**2*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + xik2j0(i)*(21*xik3j1(i)**2 + 42*xik3j1(i)*xik3j3(i) + 21*xik3j3(i)**2 - 20*xik2j2(i)*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))) - 5*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma1a**4) - 175*xik0j0(i)*(xik2j0(i)**4 - 2*xik2j0(i)**3*xik2j2(i) - 3*xik2j0(i)**2*xik2j2(i)**2 + 4*xik2j0(i)*xik2j2(i)**3 + 4*xik2j2(i)**4 - 6*xik2j0(i)*xik2j2(i)*sigma1a**4 + 3*xik2j2(i)**2*sigma1a**4 - sigma1a**8) - 630*xik1j1(i)*(xik3j1(i) + xik3j3(i))*((-2*xik2j0(i) + xik2j2(i))*sigma0a**2*sigma1a**2 + xik0j0(i)*(xik2j0(i)**2 - xik2j0(i)*xik2j2(i) - 2*xik2j2(i)**2 + sigma1a**4)) - 3*xik1j1(i)**2*(175*xik2j0(i)**3 - 350*xik2j2(i)**3 - 9*(21*xik0j0(i)*(xik3j1(i) + xik3j3(i))**2 - 10*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma0a**2*sigma1a**2) + 20*xik2j2(i)*(9*xik0j0(i)*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) - 35*sigma1a**4) - 5*xik2j0(i)*(105*xik2j2(i)**2 + 18*xik0j0(i)*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + 35*sigma1a**4)) + sigma0a**2*(-180*xik2j2(i)**2*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma0a**2 + 350*xik2j0(i)**3*sigma1a**2 + 1400*xik2j2(i)**3*sigma1a**2 + 45*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma0a**2*sigma1a**4 - 15*xik2j0(i)**2*(21*xik4j0(i)*sigma0a**2 + 30*xik4j2(i)*sigma0a**2 + 9*xik4j4(i)*sigma0a**2 + 70*xik2j2(i)*sigma1a**2) + xik2j0(i)*(-189*xik3j1(i)**2*sigma0a**2 - 378*xik3j1(i)*xik3j3(i)*sigma0a**2 - 189*xik3j3(i)**2*sigma0a**2 + 1260*xik2j2(i)*xik4j0(i)*sigma0a**2 + 1800*xik2j2(i)*xik4j2(i)*sigma0a**2 + 540*xik2j2(i)*xik4j4(i)*sigma0a**2 - 350*sigma1a**6) + 14*xik2j2(i)*(27*xik3j1(i)**2*sigma0a**2 + 54*xik3j1(i)*xik3j3(i)*sigma0a**2 + 27*xik3j3(i)**2*sigma0a**2 - 25*sigma1a**6)))/(1575.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(7,2)=(-1890*xik1j1(i)**3*(xik2j0(i) + xik2j2(i))*(xik3j1(i) + xik3j3(i)) + 135*xik1j1(i)**4*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + 3*xik0j0(i)**2*(-126*xik2j2(i)*(xik3j1(i) + xik3j3(i))**2 + 5*xik2j0(i)**2*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + 20*xik2j2(i)**2*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + xik2j0(i)*(63*xik3j1(i)**2 + 126*xik3j1(i)*xik3j3(i) + 63*xik3j3(i)**2 - 20*xik2j2(i)*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))) - 5*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma1a**4) - 175*xik0j0(i)*(xik2j0(i)**4 - 2*xik2j0(i)**3*xik2j2(i) - 3*xik2j0(i)**2*xik2j2(i)**2 + 4*xik2j0(i)*xik2j2(i)**3 + 4*xik2j2(i)**4 - 6*xik2j0(i)*xik2j2(i)*sigma1a**4 + 3*xik2j2(i)**2*sigma1a**4 - sigma1a**8) - 630*xik1j1(i)*(xik3j1(i) + xik3j3(i))*((-2*xik2j0(i) + xik2j2(i))*sigma0a**2*sigma1a**2 + xik0j0(i)*(xik2j0(i)**2 - xik2j0(i)*xik2j2(i) - 2*xik2j2(i)**2 + sigma1a**4)) - 3*xik1j1(i)**2*(175*xik2j0(i)**3 - 350*xik2j2(i)**3 + 3*(-63*xik0j0(i)*(xik3j1(i) + xik3j3(i))**2 + 10*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma0a**2*sigma1a**2) + 20*xik2j2(i)*(3*xik0j0(i)*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) - 35*sigma1a**4) - 5*xik2j0(i)*(105*xik2j2(i)**2 + 6*xik0j0(i)*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + 35*sigma1a**4)) + sigma0a**2*(-60*xik2j2(i)**2*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma0a**2 + 350*xik2j0(i)**3*sigma1a**2 + 1400*xik2j2(i)**3*sigma1a**2 + 15*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma0a**2*sigma1a**4 - 15*xik2j0(i)**2*(7*xik4j0(i)*sigma0a**2 + 10*xik4j2(i)*sigma0a**2 + 3*xik4j4(i)*sigma0a**2 + 70*xik2j2(i)*sigma1a**2) + xik2j0(i)*(-189*xik3j1(i)**2*sigma0a**2 - 378*xik3j1(i)*xik3j3(i)*sigma0a**2 - 189*xik3j3(i)**2*sigma0a**2 + 420*xik2j2(i)*xik4j0(i)*sigma0a**2 + 600*xik2j2(i)*xik4j2(i)*sigma0a**2 + 180*xik2j2(i)*xik4j4(i)*sigma0a**2 - 350*sigma1a**6) + 14*xik2j2(i)*(27*xik3j1(i)**2*sigma0a**2 + 54*xik3j1(i)*xik3j3(i)*sigma0a**2 + 27*xik3j3(i)**2*sigma0a**2 - 25*sigma1a**6)))/(1575.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(7,3)=(-945*xik1j1(i)**3*(4*xik2j0(i)*xik3j1(i) + xik2j2(i)*xik3j1(i) - xik2j0(i)*xik3j3(i) - 4*xik2j2(i)*xik3j3(i)) + 135*xik1j1(i)**4*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + 3*xik0j0(i)**2*(-126*xik2j2(i)*(3*xik3j1(i)**2 + xik3j1(i)*xik3j3(i) - 2*xik3j3(i)**2) + 5*xik2j0(i)**2*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + 20*xik2j2(i)**2*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + xik2j0(i)*(189*xik3j1(i)**2 + 63*xik3j1(i)*xik3j3(i) + 2*(-63*xik3j3(i)**2 + 10*xik2j2(i)*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i)))) + 5*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i))*sigma1a**4) - 175*xik0j0(i)*(xik2j0(i)**4 - 5*xik2j0(i)**3*xik2j2(i) + 6*xik2j0(i)**2*xik2j2(i)**2 + 4*xik2j0(i)*xik2j2(i)**3 - 8*xik2j2(i)**4 - 3*xik2j0(i)*xik2j2(i)*sigma1a**4 + 6*xik2j2(i)**2*sigma1a**4 - sigma1a**8) - 3*xik1j1(i)**2*(175*xik2j0(i)**3 - 525*xik2j0(i)**2*xik2j2(i) + 700*xik2j2(i)**3 - 3*(63*xik0j0(i)*(3*xik3j1(i)**2 + xik3j1(i)*xik3j3(i) - 2*xik3j3(i)**2) + 10*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i))*sigma0a**2*sigma1a**2) + 5*xik2j2(i)*(12*xik0j0(i)*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) - 35*sigma1a**4) - 5*xik2j0(i)*(6*xik0j0(i)*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + 35*sigma1a**4)) - 315*xik1j1(i)*((-8*xik2j0(i)*xik3j1(i) + 7*xik2j2(i)*xik3j1(i) + 2*xik2j0(i)*xik3j3(i) + 2*xik2j2(i)*xik3j3(i))*sigma0a**2*sigma1a**2 - xik0j0(i)*(2*xik2j2(i)**2*(xik3j1(i) - 4*xik3j3(i)) + xik2j0(i)**2*(-4*xik3j1(i) + xik3j3(i)) + xik2j0(i)*xik2j2(i)*(7*xik3j1(i) + 2*xik3j3(i)) + (-4*xik3j1(i) + xik3j3(i))*sigma1a**4)) + sigma0a**2*(-60*xik2j2(i)**2*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i))*sigma0a**2 + 350*xik2j0(i)**3*sigma1a**2 - 700*xik2j2(i)**3*sigma1a**2 + 15*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i))*sigma0a**2*sigma1a**4 - 15*xik2j0(i)**2*(7*xik4j0(i)*sigma0a**2 - 5*xik4j2(i)*sigma0a**2 - 12*xik4j4(i)*sigma0a**2 + 105*xik2j2(i)*sigma1a**2) + 7*xik2j2(i)*(162*xik3j1(i)**2*sigma0a**2 + 54*xik3j1(i)*xik3j3(i)*sigma0a**2 - 108*xik3j3(i)**2*sigma0a**2 + 25*sigma1a**6) - xik2j0(i)*(567*xik3j1(i)**2*sigma0a**2 + 189*xik3j1(i)*xik3j3(i)*sigma0a**2 - 378*xik3j3(i)**2*sigma0a**2 + 10*(6*xik2j2(i)*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i))*sigma0a**2 - 210*xik2j2(i)**2*sigma1a**2 + 35*sigma1a**6))))/(1575.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(7,4)=0.0d0
vinvmat(7,5)=0.0d0
vinvmat(7,6)=0.0d0
vinvmat(7,7)=(25*xik2j0(i)**4*sigma0a**2 + 100*xik2j2(i)**4*sigma0a**2 + 81*xik1j1(i)**2*xik3j1(i)**2*sigma0a**2 + 162*xik1j1(i)**2*xik3j1(i)*xik3j3(i)*sigma0a**2 + 81*xik1j1(i)**2*xik3j3(i)**2*sigma0a**2 - 200*xik0j0(i)*xik2j2(i)**3*sigma1a**2 - 270*xik1j1(i)**3*xik3j1(i)*sigma1a**2 + 27*xik0j0(i)**2*xik3j1(i)**2*sigma1a**2 - 270*xik1j1(i)**3*xik3j3(i)*sigma1a**2 + 54*xik0j0(i)**2*xik3j1(i)*xik3j3(i)*sigma1a**2 + 27*xik0j0(i)**2*xik3j3(i)**2*sigma1a**2 - 27*xik3j1(i)**2*sigma0a**4*sigma1a**2 - 54*xik3j1(i)*xik3j3(i)*sigma0a**4*sigma1a**2 - 27*xik3j3(i)**2*sigma0a**4*sigma1a**2 + 90*xik1j1(i)*xik3j1(i)*sigma0a**2*sigma1a**4 + 90*xik1j1(i)*xik3j3(i)*sigma0a**2*sigma1a**4 + 75*xik1j1(i)**2*sigma1a**6 - 25*sigma0a**2*sigma1a**8 - 50*xik2j0(i)**3*(xik2j2(i)*sigma0a**2 + xik0j0(i)*sigma1a**2) + 405*xik1j1(i)**4*sigma2a**2 - 270*xik1j1(i)**2*sigma0a**2*sigma1a**2*sigma2a**2 - 45*xik0j0(i)**2*sigma1a**4*sigma2a**2 + 45*sigma0a**4*sigma1a**4*sigma2a**2 + 10*xik0j0(i)*xik2j2(i)*(9*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma1a**2 + 5*sigma1a**6 - 54*xik1j1(i)**2*sigma2a**2) + 15*xik2j2(i)**2*(-12*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 25*xik1j1(i)**2*sigma1a**2 + 5*sigma0a**2*sigma1a**4 + 12*xik0j0(i)**2*sigma2a**2 - 12*sigma0a**4*sigma2a**2) - 15*xik2j0(i)**2*(5*xik2j2(i)**2*sigma0a**2 - 6*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 5*xik1j1(i)**2*sigma1a**2 - 10*xik0j0(i)*xik2j2(i)*sigma1a**2 + 3*(-xik0j0(i)**2 + sigma0a**4)*sigma2a**2) + 10*xik2j0(i)*(10*xik2j2(i)**3*sigma0a**2 + xik0j0(i)*(-18*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 27*xik1j1(i)**2*sigma2a**2) + 3*xik2j2(i)*(-3*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 10*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4 - 6*xik0j0(i)**2*sigma2a**2 + 6*sigma0a**4*sigma2a**2)))/(225.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(7,8)=(25*xik2j0(i)**4*sigma0a**2 + 100*xik2j2(i)**4*sigma0a**2 + 81*xik1j1(i)**2*xik3j1(i)**2*sigma0a**2 + 162*xik1j1(i)**2*xik3j1(i)*xik3j3(i)*sigma0a**2 + 81*xik1j1(i)**2*xik3j3(i)**2*sigma0a**2 - 200*xik0j0(i)*xik2j2(i)**3*sigma1a**2 - 270*xik1j1(i)**3*xik3j1(i)*sigma1a**2 + 27*xik0j0(i)**2*xik3j1(i)**2*sigma1a**2 - 270*xik1j1(i)**3*xik3j3(i)*sigma1a**2 + 54*xik0j0(i)**2*xik3j1(i)*xik3j3(i)*sigma1a**2 + 27*xik0j0(i)**2*xik3j3(i)**2*sigma1a**2 - 27*xik3j1(i)**2*sigma0a**4*sigma1a**2 - 54*xik3j1(i)*xik3j3(i)*sigma0a**4*sigma1a**2 - 27*xik3j3(i)**2*sigma0a**4*sigma1a**2 + 90*xik1j1(i)*xik3j1(i)*sigma0a**2*sigma1a**4 + 90*xik1j1(i)*xik3j3(i)*sigma0a**2*sigma1a**4 + 75*xik1j1(i)**2*sigma1a**6 - 25*sigma0a**2*sigma1a**8 - 50*xik2j0(i)**3*(xik2j2(i)*sigma0a**2 + xik0j0(i)*sigma1a**2) + 135*xik1j1(i)**4*sigma2a**2 - 90*xik1j1(i)**2*sigma0a**2*sigma1a**2*sigma2a**2 - 15*xik0j0(i)**2*sigma1a**4*sigma2a**2 + 15*sigma0a**4*sigma1a**4*sigma2a**2 + 10*xik0j0(i)*xik2j2(i)*(9*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma1a**2 + 5*sigma1a**6 - 18*xik1j1(i)**2*sigma2a**2) + 15*xik2j2(i)**2*(-12*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 25*xik1j1(i)**2*sigma1a**2 + 5*sigma0a**2*sigma1a**4 + 4*xik0j0(i)**2*sigma2a**2 - 4*sigma0a**4*sigma2a**2) - 15*xik2j0(i)**2*(5*xik2j2(i)**2*sigma0a**2 - 6*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 5*xik1j1(i)**2*sigma1a**2 - 10*xik0j0(i)*xik2j2(i)*sigma1a**2 + (-xik0j0(i)**2 + sigma0a**4)*sigma2a**2) + 10*xik2j0(i)*(10*xik2j2(i)**3*sigma0a**2 + xik0j0(i)*(-18*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 9*xik1j1(i)**2*sigma2a**2) + 3*xik2j2(i)*(-3*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 10*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4 - 2*xik0j0(i)**2*sigma2a**2 + 2*sigma0a**4*sigma2a**2)))/(225.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(7,9)=(25*xik2j0(i)**4*sigma0a**2 - 200*xik2j2(i)**4*sigma0a**2 + 243*xik1j1(i)**2*xik3j1(i)**2*sigma0a**2 + 81*xik1j1(i)**2*xik3j1(i)*xik3j3(i)*sigma0a**2 - 162*xik1j1(i)**2*xik3j3(i)**2*sigma0a**2 + 100*xik0j0(i)*xik2j2(i)**3*sigma1a**2 - 540*xik1j1(i)**3*xik3j1(i)*sigma1a**2 + 81*xik0j0(i)**2*xik3j1(i)**2*sigma1a**2 + 135*xik1j1(i)**3*xik3j3(i)*sigma1a**2 + 27*xik0j0(i)**2*xik3j1(i)*xik3j3(i)*sigma1a**2 - 54*xik0j0(i)**2*xik3j3(i)**2*sigma1a**2 - 81*xik3j1(i)**2*sigma0a**4*sigma1a**2 - 27*xik3j1(i)*xik3j3(i)*sigma0a**4*sigma1a**2 + 54*xik3j3(i)**2*sigma0a**4*sigma1a**2 + 180*xik1j1(i)*xik3j1(i)*sigma0a**2*sigma1a**4 - 45*xik1j1(i)*xik3j3(i)*sigma0a**2*sigma1a**4 + 75*xik1j1(i)**2*sigma1a**6 - 25*sigma0a**2*sigma1a**8 - 25*xik2j0(i)**3*(5*xik2j2(i)*sigma0a**2 + 2*xik0j0(i)*sigma1a**2) + 135*xik1j1(i)**4*sigma2a**2 - 90*xik1j1(i)**2*sigma0a**2*sigma1a**2*sigma2a**2 - 15*xik0j0(i)**2*sigma1a**4*sigma2a**2 + 15*sigma0a**4*sigma1a**4*sigma2a**2 - 5*xik0j0(i)*xik2j2(i)*(-9*xik1j1(i)*(7*xik3j1(i) + 2*xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 36*xik1j1(i)**2*sigma2a**2) - 30*xik2j2(i)**2*(3*xik1j1(i)*(xik3j1(i) - 4*xik3j3(i))*sigma0a**2 + 10*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4 - 2*xik0j0(i)**2*sigma2a**2 + 2*sigma0a**4*sigma2a**2) + 15*xik2j0(i)**2*(10*xik2j2(i)**2*sigma0a**2 + 3*xik1j1(i)*(4*xik3j1(i) - xik3j3(i))*sigma0a**2 - 5*xik1j1(i)**2*sigma1a**2 + 15*xik0j0(i)*xik2j2(i)*sigma1a**2 + (xik0j0(i)**2 - sigma0a**4)*sigma2a**2) + 5*xik2j0(i)*(20*xik2j2(i)**3*sigma0a**2 - 60*xik0j0(i)*xik2j2(i)**2*sigma1a**2 + 2*xik0j0(i)*(9*xik1j1(i)*(-4*xik3j1(i) + xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 9*xik1j1(i)**2*sigma2a**2) + 3*xik2j2(i)*(-3*xik1j1(i)*(7*xik3j1(i) + 2*xik3j3(i))*sigma0a**2 + 20*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4 - 4*xik0j0(i)**2*sigma2a**2 + 4*sigma0a**4*sigma2a**2)))/(225.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(7,10)=0.0d0
vinvmat(7,11)=0.0d0
vinvmat(7,12)=0.0d0
vinvmat(8,1)=(-1890*xik1j1(i)**3*(xik2j0(i) + xik2j2(i))*(xik3j1(i) + xik3j3(i)) + 135*xik1j1(i)**4*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + 3*xik0j0(i)**2*(-126*xik2j2(i)*(xik3j1(i) + xik3j3(i))**2 + 5*xik2j0(i)**2*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + 20*xik2j2(i)**2*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + xik2j0(i)*(63*xik3j1(i)**2 + 126*xik3j1(i)*xik3j3(i) + 63*xik3j3(i)**2 - 20*xik2j2(i)*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))) - 5*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma1a**4) - 175*xik0j0(i)*(xik2j0(i)**4 - 2*xik2j0(i)**3*xik2j2(i) - 3*xik2j0(i)**2*xik2j2(i)**2 + 4*xik2j0(i)*xik2j2(i)**3 + 4*xik2j2(i)**4 - 6*xik2j0(i)*xik2j2(i)*sigma1a**4 + 3*xik2j2(i)**2*sigma1a**4 - sigma1a**8) - 630*xik1j1(i)*(xik3j1(i) + xik3j3(i))*((-2*xik2j0(i) + xik2j2(i))*sigma0a**2*sigma1a**2 + xik0j0(i)*(xik2j0(i)**2 - xik2j0(i)*xik2j2(i) - 2*xik2j2(i)**2 + sigma1a**4)) - 3*xik1j1(i)**2*(175*xik2j0(i)**3 - 350*xik2j2(i)**3 + 3*(-63*xik0j0(i)*(xik3j1(i) + xik3j3(i))**2 + 10*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma0a**2*sigma1a**2) + 20*xik2j2(i)*(3*xik0j0(i)*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) - 35*sigma1a**4) - 5*xik2j0(i)*(105*xik2j2(i)**2 + 6*xik0j0(i)*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + 35*sigma1a**4)) + sigma0a**2*(-60*xik2j2(i)**2*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma0a**2 + 350*xik2j0(i)**3*sigma1a**2 + 1400*xik2j2(i)**3*sigma1a**2 + 15*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma0a**2*sigma1a**4 - 15*xik2j0(i)**2*(7*xik4j0(i)*sigma0a**2 + 10*xik4j2(i)*sigma0a**2 + 3*xik4j4(i)*sigma0a**2 + 70*xik2j2(i)*sigma1a**2) + xik2j0(i)*(-189*xik3j1(i)**2*sigma0a**2 - 378*xik3j1(i)*xik3j3(i)*sigma0a**2 - 189*xik3j3(i)**2*sigma0a**2 + 420*xik2j2(i)*xik4j0(i)*sigma0a**2 + 600*xik2j2(i)*xik4j2(i)*sigma0a**2 + 180*xik2j2(i)*xik4j4(i)*sigma0a**2 - 350*sigma1a**6) + 14*xik2j2(i)*(27*xik3j1(i)**2*sigma0a**2 + 54*xik3j1(i)*xik3j3(i)*sigma0a**2 + 27*xik3j3(i)**2*sigma0a**2 - 25*sigma1a**6)))/(1575.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(8,2)=(-1890*xik1j1(i)**3*(xik2j0(i) + xik2j2(i))*(xik3j1(i) + xik3j3(i)) + 405*xik1j1(i)**4*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + 9*xik0j0(i)**2*(-42*xik2j2(i)*(xik3j1(i) + xik3j3(i))**2 + 5*xik2j0(i)**2*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + 20*xik2j2(i)**2*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + xik2j0(i)*(21*xik3j1(i)**2 + 42*xik3j1(i)*xik3j3(i) + 21*xik3j3(i)**2 - 20*xik2j2(i)*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))) - 5*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma1a**4) - 175*xik0j0(i)*(xik2j0(i)**4 - 2*xik2j0(i)**3*xik2j2(i) - 3*xik2j0(i)**2*xik2j2(i)**2 + 4*xik2j0(i)*xik2j2(i)**3 + 4*xik2j2(i)**4 - 6*xik2j0(i)*xik2j2(i)*sigma1a**4 + 3*xik2j2(i)**2*sigma1a**4 - sigma1a**8) - 630*xik1j1(i)*(xik3j1(i) + xik3j3(i))*((-2*xik2j0(i) + xik2j2(i))*sigma0a**2*sigma1a**2 + xik0j0(i)*(xik2j0(i)**2 - xik2j0(i)*xik2j2(i) - 2*xik2j2(i)**2 + sigma1a**4)) - 3*xik1j1(i)**2*(175*xik2j0(i)**3 - 350*xik2j2(i)**3 - 9*(21*xik0j0(i)*(xik3j1(i) + xik3j3(i))**2 - 10*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma0a**2*sigma1a**2) + 20*xik2j2(i)*(9*xik0j0(i)*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) - 35*sigma1a**4) - 5*xik2j0(i)*(105*xik2j2(i)**2 + 18*xik0j0(i)*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i)) + 35*sigma1a**4)) + sigma0a**2*(-180*xik2j2(i)**2*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma0a**2 + 350*xik2j0(i)**3*sigma1a**2 + 1400*xik2j2(i)**3*sigma1a**2 + 45*(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))*sigma0a**2*sigma1a**4 - 15*xik2j0(i)**2*(21*xik4j0(i)*sigma0a**2 + 30*xik4j2(i)*sigma0a**2 + 9*xik4j4(i)*sigma0a**2 + 70*xik2j2(i)*sigma1a**2) + xik2j0(i)*(-189*xik3j1(i)**2*sigma0a**2 - 378*xik3j1(i)*xik3j3(i)*sigma0a**2 - 189*xik3j3(i)**2*sigma0a**2 + 1260*xik2j2(i)*xik4j0(i)*sigma0a**2 + 1800*xik2j2(i)*xik4j2(i)*sigma0a**2 + 540*xik2j2(i)*xik4j4(i)*sigma0a**2 - 350*sigma1a**6) + 14*xik2j2(i)*(27*xik3j1(i)**2*sigma0a**2 + 54*xik3j1(i)*xik3j3(i)*sigma0a**2 + 27*xik3j3(i)**2*sigma0a**2 - 25*sigma1a**6)))/(1575.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(8,3)=(-945*xik1j1(i)**3*(4*xik2j0(i)*xik3j1(i) + xik2j2(i)*xik3j1(i) - xik2j0(i)*xik3j3(i) - 4*xik2j2(i)*xik3j3(i)) + 135*xik1j1(i)**4*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + 3*xik0j0(i)**2*(-126*xik2j2(i)*(3*xik3j1(i)**2 + xik3j1(i)*xik3j3(i) - 2*xik3j3(i)**2) + 5*xik2j0(i)**2*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + 20*xik2j2(i)**2*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + xik2j0(i)*(189*xik3j1(i)**2 + 63*xik3j1(i)*xik3j3(i) + 2*(-63*xik3j3(i)**2 + 10*xik2j2(i)*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i)))) + 5*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i))*sigma1a**4) - 175*xik0j0(i)*(xik2j0(i)**4 - 5*xik2j0(i)**3*xik2j2(i) + 6*xik2j0(i)**2*xik2j2(i)**2 + 4*xik2j0(i)*xik2j2(i)**3 - 8*xik2j2(i)**4 - 3*xik2j0(i)*xik2j2(i)*sigma1a**4 + 6*xik2j2(i)**2*sigma1a**4 - sigma1a**8) - 3*xik1j1(i)**2*(175*xik2j0(i)**3 - 525*xik2j0(i)**2*xik2j2(i) + 700*xik2j2(i)**3 - 3*(63*xik0j0(i)*(3*xik3j1(i)**2 + xik3j1(i)*xik3j3(i) - 2*xik3j3(i)**2) + 10*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i))*sigma0a**2*sigma1a**2) + 5*xik2j2(i)*(12*xik0j0(i)*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) - 35*sigma1a**4) - 5*xik2j0(i)*(6*xik0j0(i)*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + 35*sigma1a**4)) - 315*xik1j1(i)*((-8*xik2j0(i)*xik3j1(i) + 7*xik2j2(i)*xik3j1(i) + 2*xik2j0(i)*xik3j3(i) + 2*xik2j2(i)*xik3j3(i))*sigma0a**2*sigma1a**2 - xik0j0(i)*(2*xik2j2(i)**2*(xik3j1(i) - 4*xik3j3(i)) + xik2j0(i)**2*(-4*xik3j1(i) + xik3j3(i)) + xik2j0(i)*xik2j2(i)*(7*xik3j1(i) + 2*xik3j3(i)) + (-4*xik3j1(i) + xik3j3(i))*sigma1a**4)) + sigma0a**2*(-60*xik2j2(i)**2*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i))*sigma0a**2 + 350*xik2j0(i)**3*sigma1a**2 - 700*xik2j2(i)**3*sigma1a**2 + 15*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i))*sigma0a**2*sigma1a**4 - 15*xik2j0(i)**2*(7*xik4j0(i)*sigma0a**2 - 5*xik4j2(i)*sigma0a**2 - 12*xik4j4(i)*sigma0a**2 + 105*xik2j2(i)*sigma1a**2) + 7*xik2j2(i)*(162*xik3j1(i)**2*sigma0a**2 + 54*xik3j1(i)*xik3j3(i)*sigma0a**2 - 108*xik3j3(i)**2*sigma0a**2 + 25*sigma1a**6) - xik2j0(i)*(567*xik3j1(i)**2*sigma0a**2 + 189*xik3j1(i)*xik3j3(i)*sigma0a**2 - 378*xik3j3(i)**2*sigma0a**2 + 10*(6*xik2j2(i)*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i))*sigma0a**2 - 210*xik2j2(i)**2*sigma1a**2 + 35*sigma1a**6))))/(1575.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(8,4)=0.0d0
vinvmat(8,5)=0.0d0
vinvmat(8,6)=0.0d0
vinvmat(8,7)=(25*xik2j0(i)**4*sigma0a**2 + 100*xik2j2(i)**4*sigma0a**2 + 81*xik1j1(i)**2*xik3j1(i)**2*sigma0a**2 + 162*xik1j1(i)**2*xik3j1(i)*xik3j3(i)*sigma0a**2 + 81*xik1j1(i)**2*xik3j3(i)**2*sigma0a**2 - 200*xik0j0(i)*xik2j2(i)**3*sigma1a**2 - 270*xik1j1(i)**3*xik3j1(i)*sigma1a**2 + 27*xik0j0(i)**2*xik3j1(i)**2*sigma1a**2 - 270*xik1j1(i)**3*xik3j3(i)*sigma1a**2 + 54*xik0j0(i)**2*xik3j1(i)*xik3j3(i)*sigma1a**2 + 27*xik0j0(i)**2*xik3j3(i)**2*sigma1a**2 - 27*xik3j1(i)**2*sigma0a**4*sigma1a**2 - 54*xik3j1(i)*xik3j3(i)*sigma0a**4*sigma1a**2 - 27*xik3j3(i)**2*sigma0a**4*sigma1a**2 + 90*xik1j1(i)*xik3j1(i)*sigma0a**2*sigma1a**4 + 90*xik1j1(i)*xik3j3(i)*sigma0a**2*sigma1a**4 + 75*xik1j1(i)**2*sigma1a**6 - 25*sigma0a**2*sigma1a**8 - 50*xik2j0(i)**3*(xik2j2(i)*sigma0a**2 + xik0j0(i)*sigma1a**2) + 135*xik1j1(i)**4*sigma2a**2 - 90*xik1j1(i)**2*sigma0a**2*sigma1a**2*sigma2a**2 - 15*xik0j0(i)**2*sigma1a**4*sigma2a**2 + 15*sigma0a**4*sigma1a**4*sigma2a**2 + 10*xik0j0(i)*xik2j2(i)*(9*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma1a**2 + 5*sigma1a**6 - 18*xik1j1(i)**2*sigma2a**2) + 15*xik2j2(i)**2*(-12*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 25*xik1j1(i)**2*sigma1a**2 + 5*sigma0a**2*sigma1a**4 + 4*xik0j0(i)**2*sigma2a**2 - 4*sigma0a**4*sigma2a**2) - 15*xik2j0(i)**2*(5*xik2j2(i)**2*sigma0a**2 - 6*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 5*xik1j1(i)**2*sigma1a**2 - 10*xik0j0(i)*xik2j2(i)*sigma1a**2 + (-xik0j0(i)**2 + sigma0a**4)*sigma2a**2) + 10*xik2j0(i)*(10*xik2j2(i)**3*sigma0a**2 + xik0j0(i)*(-18*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 9*xik1j1(i)**2*sigma2a**2) + 3*xik2j2(i)*(-3*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 10*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4 - 2*xik0j0(i)**2*sigma2a**2 + 2*sigma0a**4*sigma2a**2)))/(225.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(8,8)=(25*xik2j0(i)**4*sigma0a**2 + 100*xik2j2(i)**4*sigma0a**2 + 81*xik1j1(i)**2*xik3j1(i)**2*sigma0a**2 + 162*xik1j1(i)**2*xik3j1(i)*xik3j3(i)*sigma0a**2 + 81*xik1j1(i)**2*xik3j3(i)**2*sigma0a**2 - 200*xik0j0(i)*xik2j2(i)**3*sigma1a**2 - 270*xik1j1(i)**3*xik3j1(i)*sigma1a**2 + 27*xik0j0(i)**2*xik3j1(i)**2*sigma1a**2 - 270*xik1j1(i)**3*xik3j3(i)*sigma1a**2 + 54*xik0j0(i)**2*xik3j1(i)*xik3j3(i)*sigma1a**2 + 27*xik0j0(i)**2*xik3j3(i)**2*sigma1a**2 - 27*xik3j1(i)**2*sigma0a**4*sigma1a**2 - 54*xik3j1(i)*xik3j3(i)*sigma0a**4*sigma1a**2 - 27*xik3j3(i)**2*sigma0a**4*sigma1a**2 + 90*xik1j1(i)*xik3j1(i)*sigma0a**2*sigma1a**4 + 90*xik1j1(i)*xik3j3(i)*sigma0a**2*sigma1a**4 + 75*xik1j1(i)**2*sigma1a**6 - 25*sigma0a**2*sigma1a**8 - 50*xik2j0(i)**3*(xik2j2(i)*sigma0a**2 + xik0j0(i)*sigma1a**2) + 405*xik1j1(i)**4*sigma2a**2 - 270*xik1j1(i)**2*sigma0a**2*sigma1a**2*sigma2a**2 - 45*xik0j0(i)**2*sigma1a**4*sigma2a**2 + 45*sigma0a**4*sigma1a**4*sigma2a**2 + 10*xik0j0(i)*xik2j2(i)*(9*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma1a**2 + 5*sigma1a**6 - 54*xik1j1(i)**2*sigma2a**2) + 15*xik2j2(i)**2*(-12*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 25*xik1j1(i)**2*sigma1a**2 + 5*sigma0a**2*sigma1a**4 + 12*xik0j0(i)**2*sigma2a**2 - 12*sigma0a**4*sigma2a**2) - 15*xik2j0(i)**2*(5*xik2j2(i)**2*sigma0a**2 - 6*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 5*xik1j1(i)**2*sigma1a**2 - 10*xik0j0(i)*xik2j2(i)*sigma1a**2 + 3*(-xik0j0(i)**2 + sigma0a**4)*sigma2a**2) + 10*xik2j0(i)*(10*xik2j2(i)**3*sigma0a**2 + xik0j0(i)*(-18*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 27*xik1j1(i)**2*sigma2a**2) + 3*xik2j2(i)*(-3*xik1j1(i)*(xik3j1(i) + xik3j3(i))*sigma0a**2 + 10*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4 - 6*xik0j0(i)**2*sigma2a**2 + 6*sigma0a**4*sigma2a**2)))/(225.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(8,9)=(25*xik2j0(i)**4*sigma0a**2 - 200*xik2j2(i)**4*sigma0a**2 + 243*xik1j1(i)**2*xik3j1(i)**2*sigma0a**2 + 81*xik1j1(i)**2*xik3j1(i)*xik3j3(i)*sigma0a**2 - 162*xik1j1(i)**2*xik3j3(i)**2*sigma0a**2 + 100*xik0j0(i)*xik2j2(i)**3*sigma1a**2 - 540*xik1j1(i)**3*xik3j1(i)*sigma1a**2 + 81*xik0j0(i)**2*xik3j1(i)**2*sigma1a**2 + 135*xik1j1(i)**3*xik3j3(i)*sigma1a**2 + 27*xik0j0(i)**2*xik3j1(i)*xik3j3(i)*sigma1a**2 - 54*xik0j0(i)**2*xik3j3(i)**2*sigma1a**2 - 81*xik3j1(i)**2*sigma0a**4*sigma1a**2 - 27*xik3j1(i)*xik3j3(i)*sigma0a**4*sigma1a**2 + 54*xik3j3(i)**2*sigma0a**4*sigma1a**2 + 180*xik1j1(i)*xik3j1(i)*sigma0a**2*sigma1a**4 - 45*xik1j1(i)*xik3j3(i)*sigma0a**2*sigma1a**4 + 75*xik1j1(i)**2*sigma1a**6 - 25*sigma0a**2*sigma1a**8 - 25*xik2j0(i)**3*(5*xik2j2(i)*sigma0a**2 + 2*xik0j0(i)*sigma1a**2) + 135*xik1j1(i)**4*sigma2a**2 - 90*xik1j1(i)**2*sigma0a**2*sigma1a**2*sigma2a**2 - 15*xik0j0(i)**2*sigma1a**4*sigma2a**2 + 15*sigma0a**4*sigma1a**4*sigma2a**2 - 5*xik0j0(i)*xik2j2(i)*(-9*xik1j1(i)*(7*xik3j1(i) + 2*xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 36*xik1j1(i)**2*sigma2a**2) - 30*xik2j2(i)**2*(3*xik1j1(i)*(xik3j1(i) - 4*xik3j3(i))*sigma0a**2 + 10*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4 - 2*xik0j0(i)**2*sigma2a**2 + 2*sigma0a**4*sigma2a**2) + 15*xik2j0(i)**2*(10*xik2j2(i)**2*sigma0a**2 + 3*xik1j1(i)*(4*xik3j1(i) - xik3j3(i))*sigma0a**2 - 5*xik1j1(i)**2*sigma1a**2 + 15*xik0j0(i)*xik2j2(i)*sigma1a**2 + (xik0j0(i)**2 - sigma0a**4)*sigma2a**2) + 5*xik2j0(i)*(20*xik2j2(i)**3*sigma0a**2 - 60*xik0j0(i)*xik2j2(i)**2*sigma1a**2 + 2*xik0j0(i)*(9*xik1j1(i)*(-4*xik3j1(i) + xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 9*xik1j1(i)**2*sigma2a**2) + 3*xik2j2(i)*(-3*xik1j1(i)*(7*xik3j1(i) + 2*xik3j3(i))*sigma0a**2 + 20*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4 - 4*xik0j0(i)**2*sigma2a**2 + 4*sigma0a**4*sigma2a**2)))/(225.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(8,10)=0.0d0
vinvmat(8,11)=0.0d0
vinvmat(8,12)=0.0d0
vinvmat(9,1)=(-945*xik1j1(i)**3*(4*xik2j0(i)*xik3j1(i) + xik2j2(i)*xik3j1(i) - xik2j0(i)*xik3j3(i) - 4*xik2j2(i)*xik3j3(i)) + 135*xik1j1(i)**4*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + 3*xik0j0(i)**2*(-126*xik2j2(i)*(3*xik3j1(i)**2 + xik3j1(i)*xik3j3(i) - 2*xik3j3(i)**2) + 5*xik2j0(i)**2*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + 20*xik2j2(i)**2*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + xik2j0(i)*(189*xik3j1(i)**2 + 63*xik3j1(i)*xik3j3(i) + 2*(-63*xik3j3(i)**2 + 10*xik2j2(i)*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i)))) + 5*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i))*sigma1a**4) - 175*xik0j0(i)*(xik2j0(i)**4 - 5*xik2j0(i)**3*xik2j2(i) + 6*xik2j0(i)**2*xik2j2(i)**2 + 4*xik2j0(i)*xik2j2(i)**3 - 8*xik2j2(i)**4 - 3*xik2j0(i)*xik2j2(i)*sigma1a**4 + 6*xik2j2(i)**2*sigma1a**4 - sigma1a**8) - 3*xik1j1(i)**2*(175*xik2j0(i)**3 - 525*xik2j0(i)**2*xik2j2(i) + 700*xik2j2(i)**3 - 3*(63*xik0j0(i)*(3*xik3j1(i)**2 + xik3j1(i)*xik3j3(i) - 2*xik3j3(i)**2) + 10*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i))*sigma0a**2*sigma1a**2) + 5*xik2j2(i)*(12*xik0j0(i)*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) - 35*sigma1a**4) - 5*xik2j0(i)*(6*xik0j0(i)*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + 35*sigma1a**4)) - 315*xik1j1(i)*((-8*xik2j0(i)*xik3j1(i) + 7*xik2j2(i)*xik3j1(i) + 2*xik2j0(i)*xik3j3(i) + 2*xik2j2(i)*xik3j3(i))*sigma0a**2*sigma1a**2 - xik0j0(i)*(2*xik2j2(i)**2*(xik3j1(i) - 4*xik3j3(i)) + xik2j0(i)**2*(-4*xik3j1(i) + xik3j3(i)) + xik2j0(i)*xik2j2(i)*(7*xik3j1(i) + 2*xik3j3(i)) + (-4*xik3j1(i) + xik3j3(i))*sigma1a**4)) + sigma0a**2*(-60*xik2j2(i)**2*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i))*sigma0a**2 + 350*xik2j0(i)**3*sigma1a**2 - 700*xik2j2(i)**3*sigma1a**2 + 15*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i))*sigma0a**2*sigma1a**4 - 15*xik2j0(i)**2*(7*xik4j0(i)*sigma0a**2 - 5*xik4j2(i)*sigma0a**2 - 12*xik4j4(i)*sigma0a**2 + 105*xik2j2(i)*sigma1a**2) + 7*xik2j2(i)*(162*xik3j1(i)**2*sigma0a**2 + 54*xik3j1(i)*xik3j3(i)*sigma0a**2 - 108*xik3j3(i)**2*sigma0a**2 + 25*sigma1a**6) - xik2j0(i)*(567*xik3j1(i)**2*sigma0a**2 + 189*xik3j1(i)*xik3j3(i)*sigma0a**2 - 378*xik3j3(i)**2*sigma0a**2 + 10*(6*xik2j2(i)*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i))*sigma0a**2 - 210*xik2j2(i)**2*sigma1a**2 + 35*sigma1a**6))))/(1575.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(9,2)=(-945*xik1j1(i)**3*(4*xik2j0(i)*xik3j1(i) + xik2j2(i)*xik3j1(i) - xik2j0(i)*xik3j3(i) - 4*xik2j2(i)*xik3j3(i)) + 135*xik1j1(i)**4*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + 3*xik0j0(i)**2*(-126*xik2j2(i)*(3*xik3j1(i)**2 + xik3j1(i)*xik3j3(i) - 2*xik3j3(i)**2) + 5*xik2j0(i)**2*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + 20*xik2j2(i)**2*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + xik2j0(i)*(189*xik3j1(i)**2 + 63*xik3j1(i)*xik3j3(i) + 2*(-63*xik3j3(i)**2 + 10*xik2j2(i)*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i)))) + 5*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i))*sigma1a**4) - 175*xik0j0(i)*(xik2j0(i)**4 - 5*xik2j0(i)**3*xik2j2(i) + 6*xik2j0(i)**2*xik2j2(i)**2 + 4*xik2j0(i)*xik2j2(i)**3 - 8*xik2j2(i)**4 - 3*xik2j0(i)*xik2j2(i)*sigma1a**4 + 6*xik2j2(i)**2*sigma1a**4 - sigma1a**8) - 3*xik1j1(i)**2*(175*xik2j0(i)**3 - 525*xik2j0(i)**2*xik2j2(i) + 700*xik2j2(i)**3 - 3*(63*xik0j0(i)*(3*xik3j1(i)**2 + xik3j1(i)*xik3j3(i) - 2*xik3j3(i)**2) + 10*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i))*sigma0a**2*sigma1a**2) + 5*xik2j2(i)*(12*xik0j0(i)*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) - 35*sigma1a**4) - 5*xik2j0(i)*(6*xik0j0(i)*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i)) + 35*sigma1a**4)) - 315*xik1j1(i)*((-8*xik2j0(i)*xik3j1(i) + 7*xik2j2(i)*xik3j1(i) + 2*xik2j0(i)*xik3j3(i) + 2*xik2j2(i)*xik3j3(i))*sigma0a**2*sigma1a**2 - xik0j0(i)*(2*xik2j2(i)**2*(xik3j1(i) - 4*xik3j3(i)) + xik2j0(i)**2*(-4*xik3j1(i) + xik3j3(i)) + xik2j0(i)*xik2j2(i)*(7*xik3j1(i) + 2*xik3j3(i)) + (-4*xik3j1(i) + xik3j3(i))*sigma1a**4)) + sigma0a**2*(-60*xik2j2(i)**2*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i))*sigma0a**2 + 350*xik2j0(i)**3*sigma1a**2 - 700*xik2j2(i)**3*sigma1a**2 + 15*(7*xik4j0(i) - 5*xik4j2(i) - 12*xik4j4(i))*sigma0a**2*sigma1a**4 - 15*xik2j0(i)**2*(7*xik4j0(i)*sigma0a**2 - 5*xik4j2(i)*sigma0a**2 - 12*xik4j4(i)*sigma0a**2 + 105*xik2j2(i)*sigma1a**2) + 7*xik2j2(i)*(162*xik3j1(i)**2*sigma0a**2 + 54*xik3j1(i)*xik3j3(i)*sigma0a**2 - 108*xik3j3(i)**2*sigma0a**2 + 25*sigma1a**6) - xik2j0(i)*(567*xik3j1(i)**2*sigma0a**2 + 189*xik3j1(i)*xik3j3(i)*sigma0a**2 - 378*xik3j3(i)**2*sigma0a**2 + 10*(6*xik2j2(i)*(-7*xik4j0(i) + 5*xik4j2(i) + 12*xik4j4(i))*sigma0a**2 - 210*xik2j2(i)**2*sigma1a**2 + 35*sigma1a**6))))/(1575.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(9,3)=(-1890*xik1j1(i)**3*(xik2j0(i) - 2*xik2j2(i))*(3*xik3j1(i) - 2*xik3j3(i)) + 405*xik1j1(i)**4*(7*xik4j0(i) - 20*xik4j2(i) + 8*xik4j4(i)) + 9*xik0j0(i)**2*(-42*xik2j2(i)*(3*xik3j1(i) - 2*xik3j3(i))**2 + 5*xik2j0(i)**2*(7*xik4j0(i) - 20*xik4j2(i) + 8*xik4j4(i)) + 20*xik2j2(i)**2*(7*xik4j0(i) - 20*xik4j2(i) + 8*xik4j4(i)) + xik2j0(i)*(189*xik3j1(i)**2 - 252*xik3j1(i)*xik3j3(i) + 4*(21*xik3j3(i)**2 - 5*xik2j2(i)*(7*xik4j0(i) - 20*xik4j2(i) + 8*xik4j4(i)))) - 5*(7*xik4j0(i) - 20*xik4j2(i) + 8*xik4j4(i))*sigma1a**4) - 175*xik0j0(i)*(xik2j0(i)**4 - 8*xik2j0(i)**3*xik2j2(i) + 24*xik2j0(i)**2*xik2j2(i)**2 - 32*xik2j0(i)*xik2j2(i)**3 + 16*xik2j2(i)**4 - sigma1a**8) - 630*xik1j1(i)*(3*xik3j1(i) - 2*xik3j3(i))*(-2*(xik2j0(i) - 2*xik2j2(i))*sigma0a**2*sigma1a**2 + xik0j0(i)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 + sigma1a**4)) - 3*xik1j1(i)**2*(175*xik2j0(i)**3 - 1050*xik2j0(i)**2*xik2j2(i) - 1400*xik2j2(i)**3 - 9*(21*xik0j0(i)*(3*xik3j1(i) - 2*xik3j3(i))**2 - 10*(7*xik4j0(i) - 20*xik4j2(i) + 8*xik4j4(i))*sigma0a**2*sigma1a**2) + 5*xik2j0(i)*(420*xik2j2(i)**2 - 18*xik0j0(i)*(7*xik4j0(i) - 20*xik4j2(i) + 8*xik4j4(i)) - 35*sigma1a**4) + 10*xik2j2(i)*(18*xik0j0(i)*(7*xik4j0(i) - 20*xik4j2(i) + 8*xik4j4(i)) + 35*sigma1a**4)) + sigma0a**2*(-180*xik2j2(i)**2*(7*xik4j0(i) - 20*xik4j2(i) + 8*xik4j4(i))*sigma0a**2 + 350*xik2j0(i)**3*sigma1a**2 - 2800*xik2j2(i)**3*sigma1a**2 + 45*(7*xik4j0(i) - 20*xik4j2(i) + 8*xik4j4(i))*sigma0a**2*sigma1a**4 - 15*xik2j0(i)**2*(21*xik4j0(i)*sigma0a**2 - 60*xik4j2(i)*sigma0a**2 + 24*xik4j4(i)*sigma0a**2 + 140*xik2j2(i)*sigma1a**2) + xik2j0(i)*(-1701*xik3j1(i)**2*sigma0a**2 + 2268*xik3j1(i)*xik3j3(i)*sigma0a**2 - 756*xik3j3(i)**2*sigma0a**2 + 1260*xik2j2(i)*xik4j0(i)*sigma0a**2 - 3600*xik2j2(i)*xik4j2(i)*sigma0a**2 + 1440*xik2j2(i)*xik4j4(i)*sigma0a**2 + 4200*xik2j2(i)**2*sigma1a**2 - 350*sigma1a**6) + 14*xik2j2(i)*(243*xik3j1(i)**2*sigma0a**2 - 324*xik3j1(i)*xik3j3(i)*sigma0a**2 + 108*xik3j3(i)**2*sigma0a**2 + 50*sigma1a**6)))/(1575.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(9,4)=0.0d0
vinvmat(9,5)=0.0d0
vinvmat(9,6)=0.0d0
vinvmat(9,7)=(25*xik2j0(i)**4*sigma0a**2 - 200*xik2j2(i)**4*sigma0a**2 + 243*xik1j1(i)**2*xik3j1(i)**2*sigma0a**2 + 81*xik1j1(i)**2*xik3j1(i)*xik3j3(i)*sigma0a**2 - 162*xik1j1(i)**2*xik3j3(i)**2*sigma0a**2 + 100*xik0j0(i)*xik2j2(i)**3*sigma1a**2 - 540*xik1j1(i)**3*xik3j1(i)*sigma1a**2 + 81*xik0j0(i)**2*xik3j1(i)**2*sigma1a**2 + 135*xik1j1(i)**3*xik3j3(i)*sigma1a**2 + 27*xik0j0(i)**2*xik3j1(i)*xik3j3(i)*sigma1a**2 - 54*xik0j0(i)**2*xik3j3(i)**2*sigma1a**2 - 81*xik3j1(i)**2*sigma0a**4*sigma1a**2 - 27*xik3j1(i)*xik3j3(i)*sigma0a**4*sigma1a**2 + 54*xik3j3(i)**2*sigma0a**4*sigma1a**2 + 180*xik1j1(i)*xik3j1(i)*sigma0a**2*sigma1a**4 - 45*xik1j1(i)*xik3j3(i)*sigma0a**2*sigma1a**4 + 75*xik1j1(i)**2*sigma1a**6 - 25*sigma0a**2*sigma1a**8 - 25*xik2j0(i)**3*(5*xik2j2(i)*sigma0a**2 + 2*xik0j0(i)*sigma1a**2) + 135*xik1j1(i)**4*sigma2a**2 - 90*xik1j1(i)**2*sigma0a**2*sigma1a**2*sigma2a**2 - 15*xik0j0(i)**2*sigma1a**4*sigma2a**2 + 15*sigma0a**4*sigma1a**4*sigma2a**2 - 5*xik0j0(i)*xik2j2(i)*(-9*xik1j1(i)*(7*xik3j1(i) + 2*xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 36*xik1j1(i)**2*sigma2a**2) - 30*xik2j2(i)**2*(3*xik1j1(i)*(xik3j1(i) - 4*xik3j3(i))*sigma0a**2 + 10*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4 - 2*xik0j0(i)**2*sigma2a**2 + 2*sigma0a**4*sigma2a**2) + 15*xik2j0(i)**2*(10*xik2j2(i)**2*sigma0a**2 + 3*xik1j1(i)*(4*xik3j1(i) - xik3j3(i))*sigma0a**2 - 5*xik1j1(i)**2*sigma1a**2 + 15*xik0j0(i)*xik2j2(i)*sigma1a**2 + (xik0j0(i)**2 - sigma0a**4)*sigma2a**2) + 5*xik2j0(i)*(20*xik2j2(i)**3*sigma0a**2 - 60*xik0j0(i)*xik2j2(i)**2*sigma1a**2 + 2*xik0j0(i)*(9*xik1j1(i)*(-4*xik3j1(i) + xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 9*xik1j1(i)**2*sigma2a**2) + 3*xik2j2(i)*(-3*xik1j1(i)*(7*xik3j1(i) + 2*xik3j3(i))*sigma0a**2 + 20*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4 - 4*xik0j0(i)**2*sigma2a**2 + 4*sigma0a**4*sigma2a**2)))/(225.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(9,8)=(25*xik2j0(i)**4*sigma0a**2 - 200*xik2j2(i)**4*sigma0a**2 + 243*xik1j1(i)**2*xik3j1(i)**2*sigma0a**2 + 81*xik1j1(i)**2*xik3j1(i)*xik3j3(i)*sigma0a**2 - 162*xik1j1(i)**2*xik3j3(i)**2*sigma0a**2 + 100*xik0j0(i)*xik2j2(i)**3*sigma1a**2 - 540*xik1j1(i)**3*xik3j1(i)*sigma1a**2 + 81*xik0j0(i)**2*xik3j1(i)**2*sigma1a**2 + 135*xik1j1(i)**3*xik3j3(i)*sigma1a**2 + 27*xik0j0(i)**2*xik3j1(i)*xik3j3(i)*sigma1a**2 - 54*xik0j0(i)**2*xik3j3(i)**2*sigma1a**2 - 81*xik3j1(i)**2*sigma0a**4*sigma1a**2 - 27*xik3j1(i)*xik3j3(i)*sigma0a**4*sigma1a**2 + 54*xik3j3(i)**2*sigma0a**4*sigma1a**2 + 180*xik1j1(i)*xik3j1(i)*sigma0a**2*sigma1a**4 - 45*xik1j1(i)*xik3j3(i)*sigma0a**2*sigma1a**4 + 75*xik1j1(i)**2*sigma1a**6 - 25*sigma0a**2*sigma1a**8 - 25*xik2j0(i)**3*(5*xik2j2(i)*sigma0a**2 + 2*xik0j0(i)*sigma1a**2) + 135*xik1j1(i)**4*sigma2a**2 - 90*xik1j1(i)**2*sigma0a**2*sigma1a**2*sigma2a**2 - 15*xik0j0(i)**2*sigma1a**4*sigma2a**2 + 15*sigma0a**4*sigma1a**4*sigma2a**2 - 5*xik0j0(i)*xik2j2(i)*(-9*xik1j1(i)*(7*xik3j1(i) + 2*xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 36*xik1j1(i)**2*sigma2a**2) - 30*xik2j2(i)**2*(3*xik1j1(i)*(xik3j1(i) - 4*xik3j3(i))*sigma0a**2 + 10*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4 - 2*xik0j0(i)**2*sigma2a**2 + 2*sigma0a**4*sigma2a**2) + 15*xik2j0(i)**2*(10*xik2j2(i)**2*sigma0a**2 + 3*xik1j1(i)*(4*xik3j1(i) - xik3j3(i))*sigma0a**2 - 5*xik1j1(i)**2*sigma1a**2 + 15*xik0j0(i)*xik2j2(i)*sigma1a**2 + (xik0j0(i)**2 - sigma0a**4)*sigma2a**2) + 5*xik2j0(i)*(20*xik2j2(i)**3*sigma0a**2 - 60*xik0j0(i)*xik2j2(i)**2*sigma1a**2 + 2*xik0j0(i)*(9*xik1j1(i)*(-4*xik3j1(i) + xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 9*xik1j1(i)**2*sigma2a**2) + 3*xik2j2(i)*(-3*xik1j1(i)*(7*xik3j1(i) + 2*xik3j3(i))*sigma0a**2 + 20*xik1j1(i)**2*sigma1a**2 - 5*sigma0a**2*sigma1a**4 - 4*xik0j0(i)**2*sigma2a**2 + 4*sigma0a**4*sigma2a**2)))/(225.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(9,9)=(25*xik2j0(i)**4*sigma0a**2 + 400*xik2j2(i)**4*sigma0a**2 + 729*xik1j1(i)**2*xik3j1(i)**2*sigma0a**2 - 972*xik1j1(i)**2*xik3j1(i)*xik3j3(i)*sigma0a**2 + 324*xik1j1(i)**2*xik3j3(i)**2*sigma0a**2 + 400*xik0j0(i)*xik2j2(i)**3*sigma1a**2 - 810*xik1j1(i)**3*xik3j1(i)*sigma1a**2 + 243*xik0j0(i)**2*xik3j1(i)**2*sigma1a**2 + 540*xik1j1(i)**3*xik3j3(i)*sigma1a**2 - 324*xik0j0(i)**2*xik3j1(i)*xik3j3(i)*sigma1a**2 + 108*xik0j0(i)**2*xik3j3(i)**2*sigma1a**2 - 243*xik3j1(i)**2*sigma0a**4*sigma1a**2 + 324*xik3j1(i)*xik3j3(i)*sigma0a**4*sigma1a**2 - 108*xik3j3(i)**2*sigma0a**4*sigma1a**2 + 270*xik1j1(i)*xik3j1(i)*sigma0a**2*sigma1a**4 - 180*xik1j1(i)*xik3j3(i)*sigma0a**2*sigma1a**4 + 75*xik1j1(i)**2*sigma1a**6 - 25*sigma0a**2*sigma1a**8 - 50*xik2j0(i)**3*(4*xik2j2(i)*sigma0a**2 + xik0j0(i)*sigma1a**2) + 405*xik1j1(i)**4*sigma2a**2 - 270*xik1j1(i)**2*sigma0a**2*sigma1a**2*sigma2a**2 - 45*xik0j0(i)**2*sigma1a**4*sigma2a**2 + 45*sigma0a**4*sigma1a**4*sigma2a**2 - 20*xik0j0(i)*xik2j2(i)*(18*xik1j1(i)*(-3*xik3j1(i) + 2*xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 27*xik1j1(i)**2*sigma2a**2) + 15*xik2j0(i)**2*(40*xik2j2(i)**2*sigma0a**2 + 6*xik1j1(i)*(3*xik3j1(i) - 2*xik3j3(i))*sigma0a**2 - 5*xik1j1(i)**2*sigma1a**2 + 20*xik0j0(i)*xik2j2(i)*sigma1a**2 + 3*(xik0j0(i)**2 - sigma0a**4)*sigma2a**2) - 60*xik2j2(i)**2*(-6*xik1j1(i)*(3*xik3j1(i) - 2*xik3j3(i))*sigma0a**2 + 5*xik1j1(i)**2*sigma1a**2 + 3*(-xik0j0(i)**2 + sigma0a**4)*sigma2a**2) - 10*xik2j0(i)*(80*xik2j2(i)**3*sigma0a**2 + 60*xik0j0(i)*xik2j2(i)**2*sigma1a**2 - xik0j0(i)*(18*xik1j1(i)*(-3*xik3j1(i) + 2*xik3j3(i))*sigma1a**2 + 5*sigma1a**6 + 27*xik1j1(i)**2*sigma2a**2) - 6*xik2j2(i)*(-6*xik1j1(i)*(3*xik3j1(i) - 2*xik3j3(i))*sigma0a**2 + 5*xik1j1(i)**2*sigma1a**2 + 3*(-xik0j0(i)**2 + sigma0a**4)*sigma2a**2)))/(225.*(9*xik1j1(i)**4 + 6*xik1j1(i)**2*(xik0j0(i)*(xik2j0(i) - 2*xik2j2(i)) - sigma0a**2*sigma1a**2) + (xik0j0(i)**2 - sigma0a**4)*(xik2j0(i)**2 - 4*xik2j0(i)*xik2j2(i) + 4*xik2j2(i)**2 - sigma1a**4)))
vinvmat(9,10)=0.0d0
vinvmat(9,11)=0.0d0
vinvmat(9,12)=0.0d0
vinvmat(10,1)=0.0d0
vinvmat(10,2)=0.0d0
vinvmat(10,3)=0.0d0
vinvmat(10,4)=(7*xik4j0(i) + 10*xik4j2(i) + 3*xik4j4(i))/105.
vinvmat(10,5)=0.0d0
vinvmat(10,6)=0.0d0
vinvmat(10,7)=0.0d0
vinvmat(10,8)=0.0d0
vinvmat(10,9)=0.0d0
vinvmat(10,10)=sigma2a**2/15.
vinvmat(10,11)=0.0d0
vinvmat(10,12)=0.0d0
vinvmat(11,1)=0.0d0
vinvmat(11,2)=0.0d0
vinvmat(11,3)=0.0d0
vinvmat(11,4)=0.0d0
vinvmat(11,5)=(35*xik4j0(i) - 25*xik4j2(i) - 60*xik4j4(i) + (63*(xik2j0(i) + xik2j2(i))*(xik3j1(i) + xik3j3(i))**2)/(xik2j0(i)**2 + 2*xik2j0(i)*xik2j2(i) + xik2j2(i)**2 - sigma1a**4))/525.
vinvmat(11,6)=0.0d0
vinvmat(11,7)=0.0d0
vinvmat(11,8)=0.0d0
vinvmat(11,9)=0.0d0
vinvmat(11,10)=0.0d0
vinvmat(11,11)=((9*(xik3j1(i) + xik3j3(i))**2*sigma1a**2)/(5.*(xik2j0(i)**2 + 2*xik2j0(i)*xik2j2(i) + xik2j2(i)**2 - sigma1a**4)) + sigma2a**2)/15.
vinvmat(11,12)=0.0d0
vinvmat(12,1)=0.0d0
vinvmat(12,2)=0.0d0
vinvmat(12,3)=0.0d0
vinvmat(12,4)=0.0d0
vinvmat(12,5)=0.0d0
vinvmat(12,6)=(35*xik4j0(i) - 25*xik4j2(i) - 60*xik4j4(i) + (63*(xik2j0(i) + xik2j2(i))*(xik3j1(i) + xik3j3(i))**2)/(xik2j0(i)**2 + 2*xik2j0(i)*xik2j2(i) + xik2j2(i)**2 - sigma1a**4))/525.
vinvmat(12,7)=0.0d0
vinvmat(12,8)=0.0d0
vinvmat(12,9)=0.0d0
vinvmat(12,10)=0.0d0
vinvmat(12,11)=0.0d0
vinvmat(12,12)=((9*(xik3j1(i) + xik3j3(i))**2*sigma1a**2)/(5.*(xik2j0(i)**2 + 2*xik2j0(i)*xik2j2(i) + xik2j2(i)**2 - sigma1a**4)) + sigma2a**2)/15.


return
end subroutine setmat

!=======================================================================================
subroutine normalrandvec(w)
implicit none
real(8), dimension(12,1) :: w
real(4) :: normalrand
real(4) :: x1,x2,a1
integer :: i

do i=1,6
	x1=0.0
	do while(x1==0.0)
		!print*,'zero'
		call random_number(x1)
	enddo
	call random_number(x2)
	!print*,x1,x2,dble(sqrt(-2.0*log(x1))*cos(2.0*pi*x2))
	w(2*i-1,1)=dble(sqrt(-2.0*log(x1))*cos(2.0*pi*x2))
	w(2*i,1)=dble(sqrt(-2.0*log(x1))*sin(2.0*pi*x2))
enddo

return
end subroutine normalrandvec

!=======================================================================================
function normalrand()
implicit none
real(4) :: x1,x2,a1
real(4) :: normalrand

call random_number(x1)
do while(x1==0.0)
	print*,'zero'
	call random_number(x1)
enddo
call random_number(x2)
normalrand=dble(sqrt(-2.0*log(x1))*cos(2.0*pi*x2))

return
end function normalrand

SUBROUTINE jacobi(a,d,v,nrot)
IMPLICIT NONE
INTEGER(4), INTENT(OUT) :: nrot
REAL(4), DIMENSION(:), INTENT(OUT) :: d
REAL(4), DIMENSION(:,:), INTENT(INOUT) :: a
REAL(4), DIMENSION(:,:), INTENT(OUT) :: v
INTEGER(4) :: i,ip,iq,n
REAL(4) :: c,g,h,s,sm,t,tau,theta,tresh
REAL(4), DIMENSION(size(d)) :: b,z
call unit_matrix(v(:,:))
b(:)=get_diag(a(:,:))

d(:)=b(:)
z(:)=0.0
nrot=0
n=3
do i=1,50
	
   sm=sum(abs(a),mask=upper_triangle(n,n))
!print*,i,sm
   if (sm == 0.0) RETURN
   tresh=merge(0.2*sm/real(n)**2.0,0.0, i < 4 )
   do ip=1,n-1
      do iq=ip+1,n
         g=100.0*abs(a(ip,iq))
         if ((i > 4) .and. (abs(d(ip))+g == abs(d(ip))) &
            .and. (abs(d(iq))+g == abs(d(iq)))) then
            a(ip,iq)=0.0
         else if (abs(a(ip,iq)) > tresh) then
            h=d(iq)-d(ip)
            if (abs(h)+g == abs(h)) then
               t=a(ip,iq)/h
            else
               theta=0.5*h/a(ip,iq)
               t=1.0/(abs(theta)+sqrt(1.0+theta**2))
               if (theta < 0.0) t=-t
            end if
            c=1.0/sqrt(1+t**2)
            s=t*c
            tau=s/(1.0+c)
            h=t*a(ip,iq)
            z(ip)=z(ip)-h
            z(iq)=z(iq)+h
            d(ip)=d(ip)-h
            d(iq)=d(iq)+h
            a(ip,iq)=0.0
            call jrotate(a(1:ip-1,ip),a(1:ip-1,iq),s,tau)
            call jrotate(a(ip,ip+1:iq-1),a(ip+1:iq-1,iq),s,tau)
            call jrotate(a(ip,iq+1:n),a(iq,iq+1:n),s,tau)
            call jrotate(v(:,ip),v(:,iq),s,tau)
            nrot=nrot+1
         end if
      end do
   end do
   b(:)=b(:)+z(:)
   d(:)=b(:)
   z(:)=0.0
end do

END SUBROUTINE jacobi

!BL
SUBROUTINE jrotate(a1,a2,s,tau)
REAL(4), DIMENSION(:), INTENT(INOUT) :: a1,a2
REAL(4), DIMENSION(size(a1)) :: wk1
real(4) :: s,tau
wk1(:)=a1(:)
a1(:)=a1(:)-s*(a2(:)+a1(:)*tau)
a2(:)=a2(:)+s*(wk1(:)-a2(:)*tau)
END SUBROUTINE jrotate


FUNCTION get_diag(mat)
REAL(4), DIMENSION(:,:), INTENT(IN) :: mat
REAL(4), DIMENSION(size(mat,1)) :: get_diag
INTEGER(4) :: j
do j=1,size(mat,1)
   get_diag(j)=mat(j,j)
end do
END FUNCTION get_diag

SUBROUTINE unit_matrix(mat)
REAL(4), DIMENSION(:,:), INTENT(OUT) :: mat
INTEGER(4) :: i,n
n=min(size(mat,1),size(mat,2))
mat(:,:)=0.0
do i=1,n
   mat(i,i)=1.0
end do
END SUBROUTINE unit_matrix
!BL
FUNCTION upper_triangle(j,k,extra)
INTEGER(4), INTENT(IN) :: j,k
INTEGER(4), OPTIONAL, INTENT(IN) :: extra
LOGICAL, DIMENSION(j,k) :: upper_triangle
INTEGER(4) :: n
n=0
if (present(extra)) n=extra
upper_triangle=(outerdiff_i(arth_i(1,1,j),arth_i(1,1,k)) <n)
END FUNCTION upper_triangle

FUNCTION arth_i(first,increment,n)
INTEGER(4), INTENT(IN) :: first,increment,n
INTEGER(4), DIMENSION(n) :: arth_i
INTEGER(4) :: k,k2,temp
if (n > 0) arth_i(1)=first
if (n <= NPAR_ARTH) then
   do k=2,n
      arth_i(k)=arth_i(k-1)+increment
   end do
else
   do k=2,NPAR2_ARTH
      arth_i(k)=arth_i(k-1)+increment
   end do
   temp=increment*NPAR2_ARTH
   k=NPAR2_ARTH
   do
      if (k >= n) exit
      k2=k+k
      arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
      temp=temp+temp
      k=k2
   end do
end if
END FUNCTION arth_i

FUNCTION outerdiff_i(a,b)
INTEGER(4), DIMENSION(:), INTENT(IN) :: a,b
INTEGER(4), DIMENSION(size(a),size(b)) :: outerdiff_i
outerdiff_i = spread(a,dim=2,ncopies=size(b)) - &
   spread(b,dim=1,ncopies=size(a))
END FUNCTION outerdiff_i


end program pksampler


  
