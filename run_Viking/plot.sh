cmap='wrgb=0,0,1 grgb=1,1,0 brgb=1,0,0'
psimage < vp_init n1=201 d1=20 n2=1568 d2=20 d1num=1000 d2num=5000 label1="Z (m)" label2="X (m)" title='(a) Vp' legend=1 lstyle=vertright lheight=4 width=6. height=4. $cmap> vp_init.ps

psimage < rho_init n1=201 d1=20 n2=1568 d2=20 d1num=1000 d2num=5000 label1="Z (m)" label2="X (m)" title='(b) rho' legend=1 lstyle=vertright lheight=4 width=6. height=4. $cmap> rho_init.ps

psmerge in=vp_init.ps translate=0,0 in=rho_init.ps translate=8,0 > vp_rho.ps
epstopdf vp_rho.ps

#suflip < dat_0001.su flip=2 |sustrip  > dat_0001
sustrip < dat_0001.su > dat_0001
psimage < dat_0001 n1=1500 d1=0.004 n2=120 d2=25 label1='Time (s)' label2='Distance (m)' title='(a)' perc=99> dat_0001.ps

pswigb < stf n1=2500 d1=0.0024 label1='Time (s)' title='(b)'> stf.ps

psmerge in=dat_0001.ps translate=0,0 in=stf.ps translate=7,0 > dat_wlt.ps
epstopdf dat_wlt.ps

rm *.ps
