!!$ deduce the density from velocity by gardner's law
!========================================================
program gardner
 implicit none
 
 integer::n1,n2,n3,i1,i2,i3
 real ::h,depthwater,rhowater
 real,dimension(:,:,:),allocatable::x,x1
 character(len=80):: fvp,frho

 read(*,*) fvp,frho
 read(*,*) n1,n2,n3
 read(*,*) h,depthwater,rhowater

 allocate (x(n1,n2,n3))
 allocate (x1(n1,n2,n3))

 open(10,file=fvp,recl=n1*n2*n3*4,access='direct')
 read(10,rec=1) x
 close(10)

 do i3=1,n3 
    do i2=1,n2
       do i1=1,n1
          if ( (i1-1)*h<=depthwater) then
             x1(i1,i2,i3)=rhowater
          else
             x1(i1,i2,i3)=1741*(0.001*x(i1,i2,i3))**0.25
          endif
       enddo
    enddo
 enddo

 open(11,file=frho,recl=n1*n2*n3*4,access='direct')
 write(11,rec=1) x1 
 close(11)

 deallocate(x)
 deallocate(x1)
end program gardner
