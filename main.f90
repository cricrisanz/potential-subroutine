program testvpot
  implicit none
  real(kind=8) :: br,sr,t,v,dbr,dt,e
  integer, parameter :: nbr=100,nt=10,np=3500
  real(kind=8), parameter :: bri=4.0625d0,brf=40.d0,ti=0.0001d0,tf=90.d0
  integer :: ibr,it,ip
  
  dbr=(brf-bri)/dble(nbr-1)
  dt=(tf-ti)/dble(nt-1)

  open(8,file='fit.inp',status='old',action='read')
  t=ti
!  do it = 1,nt
!     br=bri
!     do ibr =1,nbr
  do ip=1,np
      read(8,*)br,sr,t,e
      call vpot(br,t,v)
      write(*,'(2(f8.4,2x),2(f18.12,1x),f14.8)') br,t,v,e,(v-e)*219474.63
      !      br=br+dbr
      !     enddo
!     t=t+dt
  enddo

  stop
end program testvpot
  
  
  
