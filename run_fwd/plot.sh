psimage < snapshot.bin n1=141 d1=5 d2=5 label1='Z (m)' label2='X (m)' height=6. width=6. perc=90 >snapshot.ps

epstopdf snapshot.ps
rm *.ps
