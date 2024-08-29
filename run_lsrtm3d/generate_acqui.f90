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
  integer :: nxs = 16
  integer :: nys = 16
  integer :: nxr = 200
  integer :: nyr = 200
  real :: dxr = 25
  real :: dyr = 25
  real :: oxr = 10
  real :: oyr = 10
  real :: dxs = 300
  real :: dys = 300
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
