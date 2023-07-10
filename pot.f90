!!!
! 7 July 2023: Potential expansion in Legendre Polinomias
! JCP, 139, 144305 (2013)
! For any question ask: cristina.sanz@uam.es
! Be careful and be sure that the values of k1, l1, k2,l2 and n
! correspond to the values used during your fitting
! E0 has to be changed to the value corresponding to long range of your system.
!!!

subroutine vpot(br,the,v)
! program potential
  implicit none
 !in and out variables
  real(kind=8), intent(in) :: br,the
  real(kind=8), intent(out) :: v
!  real(kind=8) :: sr, br,the
!  real(kind=8) :: v
 !parameter for potential expansion, short range summations
  integer, dimension(4), parameter :: k1=(/0,1,2,3/)
  integer, dimension(3), parameter :: l1=(/0,2,4/)        
 !parameter for potential expansion, long rage summations
  integer, dimension(1), parameter :: k2=(/0/)
  integer, dimension(2), parameter :: l2=(/0,2/)
  integer, dimension(3), parameter :: n=(/6,8,10/)

 !other parameters
  integer, parameter :: nvex1=size(k1) !non-linear parameters first summation
  integer, parameter :: nvex2=1 !non-linear parameter second summation
  integer, parameter :: nk1=size(k1),nl1=size(l1),nk2=size(k2),nl2=size(l2),nn=size(n)
  real(kind=8), parameter :: E0=-528.266863429459d0

  !other variables needed in subroutine
  integer, parameter :: nsum=nk1*nl1*3
  real(kind=8), dimension(nsum) :: xx
  real(kind=8), dimension(nk1,nl1) :: f,g,h
  real(kind=8), dimension(nvex1), save :: vex1
  real(kind=8), dimension(nvex2), save :: vex2
  real(kind=8) :: q !second expansion
  real(kind=8) :: sr,cthet,auxsr,auxbr,auxth, polleg
  real(kind=8), parameter :: pi=dacos(-1.d0)
  real(kind=8), dimension(nsum), save :: sm
  real(kind=8), dimension(nk2,nl2,nn), save :: cnkl !non linear coefficients of second term
  logical, save :: first=.true.
  integer :: i,j,k,l,nnn,ixx
  interface
   real(kind=8) function pleg(l,x)
     implicit none
     integer, intent(in) :: l
     real(kind=8), intent(in) :: x
   end function pleg

  real(kind=8) function fn(n,d,x)
    implicit none
    integer, intent(in) :: n
    real(kind=8), intent(in) :: d,x
  end function fn
  
  real(kind=8) function fact(nf)
    implicit none
    integer, intent(in) :: nf
  end function fact
  
  end interface

  !sr is used but it is not need at all, rigid rotor
  sr=1.4d0 !For the fitting only considering R and theta sr is not really used

  if(first)then
     open(13,file="coeffs-lr.txt",status="old",action="read")
     do k=1,nk2
        do l=1,nl2
           do nnn=1,nn
              read(13,*)cnkl(k,l,nnn)
!              print*,"Cnkl",k,l,nnn,cnkl(k,l,nnn)
           enddo
        enddo
     enddo

     open(14,file="coeffs-sr.txt",status="old",action="read")
     do i=1,nvex1
        read(14,*)vex1(i)
!        print*,vex1(i)
     enddo
     do i=1,nvex2
        read(14,*)vex2(i)
!        print*,vex2(i)
     enddo
     ixx=1
     do i=1,nk1
        do j=1,nl1
           read(14,*)sm(ixx)
           ixx=ixx+1
        enddo
        do j=1,nl1
           read(14,*)sm(ixx)
           ixx=ixx+1
        enddo
        do j=1,nl1
           read(14,*)sm(ixx)
           ixx=ixx+1
        enddo
     enddo
     first=.false.
  endif
!!!!!!!!!!!!!!!!!!!!
!Second summation term
!!!!!!!!!!!!!!!!!!!!
  
  cthet=dcos(the*pi/180.d0)
  q=0.d0
  do i=1,nk2
     auxsr=sr**k2(i)
    !caution: when more than 1 non-linear parameter in the second term
     do j=1,nl2
        polleg=pleg(l2(j),cthet)
        do l=1,nn
           q=q+auxsr * fn(n(l),vex2(nvex2),br)*cnkl(i,j,l)*polleg*(br**(-n(l)))
        enddo
     enddo
  enddo
  
!!!!!!!!!!!!!!!!!!!!
!First summation term
!!!!!!!!!!!!!!!!!!!!

  v=0.d0
  cthet=dcos(the*pi/180.d0)
  do i=1,nk1
     auxsr = sr**k1(i)
     auxbr = dexp(-vex1(i)*br)     
     do j=1,nl1
        f(i,j)= auxsr * auxbr * pleg(l1(j),cthet)
        g(i,j)= br * f(i,j)
        h(i,j)= br * g(i,j)
     enddo
  enddo
  ixx=1
  do i=1,nk1  
   do j=1,nl1
      xx(ixx)=f(i,j)
      ixx=ixx+1
   enddo 
   do j=1,nl1
      xx(ixx)=g(i,j)
      ixx=ixx+1
   enddo 
   do j=1,nl1
      xx(ixx)=h(i,j)
      ixx=ixx+1
   enddo 
   xx(ixx-1)=xx(ixx-1)-q
  enddo

  v = 0.0d0
  do k=1,nsum
   v = v + sm(k)*xx(k)
  enddo
  v=v+E0

  return
end subroutine vpot


  !---------- P0
  
 real(kind=8) function pleg(l,x)
   implicit none
   integer, intent(in) :: l
   real(kind=8), intent(in) :: x

  select case(l)
     case(0)
        pleg=1.d0
     case(1)
        pleg=x
     case(2)
        pleg = 0.5d0*(3.d0*x*x - 1.d0)
     case(3)
        pleg = 0.5d0*(5.0d0*x*x*x - 3.d0*x)
     case(4)
        pleg = 0.125d0*(35.0d0*x*x*x*x - 30.d0*x*x + 3.d0)
     case(5)
        pleg = 0.125d0*(63.0d0*x*x*x*x*x - 70.d0*x*x*x + 15.d0*x)
     case(6)
        pleg = 0.0625d0*(231.0d0*x*x*x*x*x*x - 315.d0*x*x*x*x + &
       &       105.d0*x*x - 5.d0)
     case(7)
        pleg = 0.0625d0*(429.0d0*x*x*x*x*x*x*x - 693.d0*x*x*x*x*x + &
      &       315.d0*x*x*x - 35.d0*x)  
     end select
  return
 end function pleg

  !---------- Tang-Toennies damping function, K.T.Tang, J. P. Toennies, JCP. 80, 3726 (1984)

  real(kind=8) function fn(n,d,x)
  implicit none
  integer, intent(in) :: n
  real(kind=8), intent(in) :: d,x
  real(kind=8) :: aux
  integer :: i
  interface
  real(kind=8) function fact(nf)
    implicit none
    integer, intent(in) :: nf
  end function fact
  end interface

  aux=1.d0 !!First term is always 1, factorial of 0
  do i=1,n
     aux=aux+((d*x)**i)/fact(i)
  enddo   
  fn=1.d0-dexp(-d*x)*aux
  
  return
  end function fn

  real(kind=8) function fact(nf)
  implicit none
  integer, intent(in) :: nf
  integer :: i

  if (nf < 0) error stop 'factorial is singular for negative integers'
  fact = 1.0d0
  do i = 2, nf
    fact = fact * dble(i)
  enddo
end function fact
