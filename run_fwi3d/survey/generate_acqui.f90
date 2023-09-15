!===========================================================================
!create the acquisition file 
!Author: Pengliang Yang
!E-mail: ypl.2100@gmail.com
!==========================================================================
program createacqui
  implicit none

  integer :: isrc
  real :: xsrc, ysrc, zsrc
  real :: xrec, yrec, zrec
  integer :: iys, ixs
  integer :: iyr, ixr
  integer :: nxs = 20
  integer :: nys = 20
  integer :: nxr = 100
  integer :: nyr = 100
  real :: dxr = 100
  real :: dyr = 100
  real :: oxr = 50
  real :: oyr = 50
  real :: dxs = 500
  real :: dys = 500
  real :: oxs = 250
  real :: oys = 250
  
  zsrc = 5
  zrec = 10
  
  open(10,file='acqui.txt',status='replace')
  write(10,*) 'z x y azimuth dip src/rec(0/1)'
  do iys = 1,nys
     ysrc = oys + (iys-1)*dys
     do ixs = 1,nxs
        xsrc = oxs + (ixs-1)*dxs

        write(10,*) zsrc,xsrc,ysrc,0,0,0 !source code=0
        do iyr = 1,nyr
           yrec = oyr + (iyr-1)*dyr
           do ixr = 1,nxr
              xrec = oxr + (ixr-1)*dxr

              write(10,*) zrec,xrec,yrec,0,0,1 !receiver code=1
           enddo
        enddo
     enddo
  enddo
  close(10)

  open(10,file='sources.txt',status='replace')
  do iys = 1,nys
     ysrc = oys + (iys-1)*dys
     do ixs = 1,nxs
        xsrc = oxs + (ixs-1)*dxs
        
        isrc = ixs + nxs*(iys-1)

        write(10,*) zsrc,xsrc,ysrc,0,0,0,isrc !source code=0
     enddo
  enddo
  close(10)

  open(10,file='receivers.txt',status='replace')
  do iyr = 1,nyr
     yrec = oyr + (iyr-1)*dyr
     do ixr = 1,nxr
        xrec = oxr + (ixr-1)*dxr

        write(10,*) zrec,xrec,yrec,0,0,1 !receiver code=1
     enddo
  enddo
  close(10)
  
end program createacqui
