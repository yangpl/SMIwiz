!!$ deduce the density from velocity by gardner's law
program gardner
 implicit none
 
 integer::n1,n2,i1,i2
 real ::h,velwater,depthwater
 real,dimension(:,:),allocatable::x,x1
 character(len=80):: fvp,frho

 read(*,*) fvp,frho
 read(*,*) n1,n2 !n1 (fast) 
 read(*,*) h,depthwater,velwater

 allocate (x(n1,n2))
 allocate (x1(n1,n2))

 open(10,file=fvp,recl=n1*n2*4,access='direct')
 read(10,rec=1) x
 close(10)

 
 do i2=1,n2
    do i1=1,n1
       if ( (i1-1)*h<=depthwater) then
          x1(i1,i2)=velwater
       else
          x1(i1,i2)=1741*(0.001*x(i1,i2))**0.25
       endif
    enddo
 enddo

 open(11,file=frho,recl=n1*n2*4,access='direct')
 write(11,rec=1) x1 
 close(11)

 deallocate(x)
 deallocate(x1)
end program gardner
