gfortran smooth2d.f90
./a.out<<EOF
vp_marm_true vp_marm_init
151 461 5 5 8
EOF

gfortran gardner.f90
./a.out<<EOF
vp_marm_true rho_true
151 461 1
20. 0. 1000.
EOF

./a.out<<EOF
vp_marm_init rho_init
151 461 1
20. 0. 1000.
EOF
