psimage < wave_up.bin n1=141 n2=141 title="(a) Up-going" height=6. width=6.> wave_up.eps
psimage < wave_down.bin n1=141 n2=141 title="(b) Down-going" height=6. width=6. > wave_down.eps

psmerge in=wave_up.eps translate=0,0 in=wave_down.eps translate=7,0 > updown.eps

epstopdf updown.eps
rm *.eps
