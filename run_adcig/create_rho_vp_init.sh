gfortran smooth2d.f90
./a.out<<EOF
vp.rsf@ vp_true
251 767 4 4 1
12. 40.
EOF
./a.out<<EOF
vp.rsf@ vp_init
251 767 4 4 4
12. 40.
EOF

./a.out<<EOF
rho.rsf@ rho_true
251 767 4 4 1
12. 40.
EOF
./a.out<<EOF
rho.rsf@ rho_init
251 767 4 4 4
12. 40.
EOF

# suaddhead< vp.rsf@ ns=251 |suk1k2filter d1=12 d2=12 k1=0,0.05 k2=0,0.05 amps1=1,0 amps2=1,0 |sustrip > vp_true
# suaddhead< rho.rsf@ ns=251 |suk1k2filter d1=12 d2=12 k1=0,0.05 k2=0,0.05 amps1=1,0 amps2=1,0 |sustrip > rho_true

# suaddhead< vp.rsf@ ns=251 |suk1k2filter d1=12 d2=12 k1=0,0.05 k2=0,0.03 amps1=1,0 amps2=1,0 |sustrip > vp_init
# suaddhead< rho.rsf@ ns=251 |suk1k2filter d1=12 d2=12 k1=0,0.05 k2=0,0.03 amps1=1,0 amps2=1,0 |sustrip > rho_init
 
