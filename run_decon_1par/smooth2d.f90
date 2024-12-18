!=================================================================
module triangle
  implicit none

  integer::np,nq
  real,dimension(:),allocatable::bb1,bb2,pp,qq,yy

contains
  subroutine triangle_init(nbox,nd)
    integer::nbox,nd

    if(nbox < 1) stop 'nbox<1'
    if(nbox > nd) stop 'nbox>nx'
    np=nbox+nd-1
    nq=nbox+np-1
    allocate(pp(np))
    allocate(qq(nq))
    allocate(bb1(np))
    allocate(bb2(nq))
	allocate(yy(nd))
  end subroutine triangle_init

  subroutine triangle_close()
    deallocate(pp)
    deallocate(qq)
    deallocate(bb1)
    deallocate(bb2)
    deallocate(yy)
  end subroutine triangle_close

  subroutine boxconv_lop(nbox, nx, xx, yy, bb)
    integer,intent(in):: nx,nbox
    real,dimension(nx),intent(in)::xx
    real,dimension(nx+nbox-1),intent(out)::yy,bb

    integer::i,ny

    ny=nx+nbox-1
    do i= 1, ny
       bb(i) = 0.
    enddo
    bb(1) = xx(1)
    do i= 2, nx
       bb(i) = bb(i-1) + xx(i)  ! make B(Z) = X(Z)/(1-Z)
    enddo
    do i= nx+1, ny
       bb(i) = bb(i-1)
    enddo
    do i= 1, nbox
       yy(i) = bb(i)
    enddo
    do i= nbox+1, ny
       yy(i) = bb(i) - bb(i-nbox) ! make Y(Z) = B(Z)*(1-Z**nbox)
    enddo
    do i= 1, ny
       yy(i) = yy(i) / nbox
    enddo
  end subroutine boxconv_lop

  subroutine triangle_lop(nbox,nd,xx)
    integer, intent(in)::nbox,nd
    real,dimension(nd),intent(inout)::xx

    integer::i

    call boxconv_lop(nbox,nd,xx,pp,bb1)
    call boxconv_lop(nbox,np,pp,qq,bb2)
    do i=1,nd
       yy(i)=qq(i+nbox-1)
    enddo
    do i=1,nbox-1
       yy(i)=yy(i)+qq(nbox-i) !fold back near end
    enddo
    do i=1,nbox-1
       yy(nd-i+1)=yy(nd-i+1)+qq(nd+(nbox-1)+i) !fold back far end
    enddo
    xx(:)=yy(:)
  end subroutine triangle_lop
end module triangle

!!$===========================================================================
program smooth2d
  use triangle
  implicit none

  integer::n1,n2,r1,r2, repeat
  real,dimension(:,:),allocatable::mod0
  real,dimension(:),allocatable:: bathy
  real :: val
  integer::i1,i2,i
  character(len=80)::inputfile,outputfile

  read(*,*) inputfile,outputfile
  read(*,*) n1,n2,r1,r2,repeat !size of data, smoothing radius,repeating times

  
  allocate(mod0(n1,n2))
  allocate(bathy(n2))
  
  open(10,file=inputfile,access='direct',recl=4*n1*n2,status='old')
  read(10,rec=1) mod0
  close(10)

  val = mod0(1,1)
  do i=1,repeat
     !----------------------------------------
     call triangle_init(r1,n1)
     do i2=1,n2
        call triangle_lop(r1,n1,mod0(:,i2))
     enddo
     call triangle_close()

     !-----------------------------------------
     call triangle_init(r2,n2)
     do i1=1,n1
        call triangle_lop(r2,n2,mod0(i1,:))
     enddo
     call triangle_close()
  enddo

  open(10,file=outputfile,access='direct',recl=4*n1*n2,status='replace')
  write(10,rec=1) mod0
  close(10)


  deallocate(mod0)
end program smooth2d
