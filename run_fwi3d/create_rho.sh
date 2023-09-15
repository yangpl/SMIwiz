gfortran gardner.f90 -o gardner
./gardner<<EOF
vp rho
61 201 201
50. 0. 1000.
EOF

./gardner<<EOF
vp_init rho_init
61 201 201
50. 0. 1000.
EOF
