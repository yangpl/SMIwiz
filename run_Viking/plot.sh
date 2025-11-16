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


dd if=param_final_rtm of=rtm_ip skip=1 bs=1260672
dd if=param_final of=lsrtm_ip skip=1 bs=1260672

psimage < rtm_ip n1=201 d1=20 n2=1568 d2=20 d1num=1000 d2num=5000 perc=99 label1="Z (m)" label2="X (m)" title='(a)' width=6. height=4.> rtm_ip.ps
psimage < lsrtm_ip n1=201 d1=20 n2=1568 d2=20 d1num=1000 d2num=5000 perc=99 label1="Z (m)" label2="X (m)" title='(b)' width=6. height=4.> lsrtm_ip.ps

psimage < param_final_psf n1=201 d1=20 n2=1568 d2=20 d1num=1000 d2num=5000 perc=99 label1="Z (m)" label2="X (m)" title='(c)' width=6. height=4.> psf_ip.ps
psimage < param_final_fft n1=201 d1=20 n2=1568 d2=20 d1num=1000 d2num=5000 perc=99 label1="Z (m)" label2="X (m)" title='(d)' width=6. height=4.> fft_ip.ps

psmerge in=rtm_ip.ps translate=0,0 in=lsrtm_ip.ps translate=8,0 in=psf_ip.ps translate=0,-5.5 in=fft_ip.ps translate=8,-5.5 > viking_rtm_lsrtm.ps
epstopdf viking_rtm_lsrtm.ps


rm *.ps
