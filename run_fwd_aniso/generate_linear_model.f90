program main
  implicit none

  integer :: n1, n2, i1, i2
  real, dimension(:,:), allocatable :: vp, vpfinal
  real, dimension(:), allocatable :: vpnew
  real :: s
  
  n1 = 151
  n2 = 461

  allocate(vp(n1, n2))
  allocate(vpfinal(n1, n2))
  allocate(vpnew(n1))

  open(10,file='vp_marm_true',action='read',access='direct',recl=4*n1*n2,status='old')
  read(10,rec=1) vp(:,:)
  close(10)

  
  do i1=1,n1
     s = 0
     do i2=1,n2
        s = s + vp(i1,i2)
     enddo
     vpnew(i1) = s/n2
  enddo

  do i2=1,n2
     do i1=1,n1
        if(i1<4) then
           vpfinal(i1,i2) = vp(i1,i2)
        else
           vpfinal(i1,i2) = vpnew(4) + (vpnew(n1)-vpnew(4) )*(i1-4)/(n1-4)
        endif
     enddo
  enddo
  
  open(55,file='vp_linear',access='direct',status='replace',recl=4*n1*n2)
  write(55,rec=1) vpfinal(:,:)
  close(55)



  deallocate(vp)
  deallocate(vpnew)
  deallocate(vpfinal)
end program main
