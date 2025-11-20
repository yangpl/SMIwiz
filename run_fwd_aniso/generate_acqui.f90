!===========================================================================
!create the acquisition file for toyxd_time
!Pengliang Yang,
!E-mail: pengliang.yang@univ-grenoble-alpes.fr
!==========================================================================
program createacqui
  implicit none

  real::xmin,xmax,zmin,zmax
  real::zsou,xsou,zrec,xrec
  real::xs1,xs2,dxs
  real::xr1,xr2,dxr
  real::zr1,zr2,dzr
  !character(len=80)::filename
  integer::ns !number of sources
  integer::nr, nrz !number of receivers per source asigned by offset
  integer::is,ir

  zmin=0.
  zmax=3000.
  xmin=0.
  xmax=9200
  zsou=5
  zrec=20
  xs1=100.
  xs2=9100
  dxs=380
  xr1=100.
  xr2=9100.
  dxr=20.

  zr1=zrec + 50
  zr2=zmax - 50
  dzr=20.
  
  ns=int((xs2-xs1)/dxs) + 1
  nr=int((xr2-xr1)/dxr) + 1
  nrz=int((zr2-zr1)/dxr) + 1
  print *, 'ns=', ns
  print *, 'nr=', nr
  print *, 'nrz=', nrz
  
  
  open(10, file='acqui.txt', status='replace')
  write(10,*) 'z     x    y     azimuth    dip    src/rec(0/1)'
  do is=1,ns
     xsou=xs1+(is-1)*dxs
     write(10,*) zsou,xsou,0,0,0,0 !source code=0

     do ir=1,nr
        xrec=xr1+(ir-1)*dxr
        write(10,*) zrec,xrec,0,0,0,1 !receiver code=1
     enddo

     !two vertical lines of receivers
     !================================
     ! xrec = xr1
     ! do ir=1,nrz
     !    zrec=zr1+(ir-1)*dzr
     !    write(10,*) zrec,xrec,0,0,0,1 !receiver code=1
     ! enddo

     ! xrec = xr2
     ! do ir=1,nrz
     !    zrec=zr1+(ir-1)*dzr
     !    write(10,*) zrec,xrec,0,0,0,1 !receiver code=1
     ! enddo
     
  enddo
  close(10)
end program createacqui
