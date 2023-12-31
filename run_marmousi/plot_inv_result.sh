#colormap
colormap1='wrgb=1.0,0,0 grgb=1.0,1.0,1.0 brgb=0,0,1.0' 
colormap2='wrgb=0,0,1 grgb=1,1,0 brgb=1,0,0' #used by Ludovic
colormap3='bhls=0.666666,.5,1  ghls=.333333,.5,1 whls=0,.5,1 bps=24 ' #used by Wei Zhou
colormap4='bhls=.2,.5,1  ghls=0,.5,1 whls=.8,.5,1 bps=24' #used by Wei Zhou

option="width=7 height=4 label1=Z(km) label2=X(km) d1num=1 d2num=3 legend=1 lheight=4 lwidth=0.3 lstyle=vertright "


dd if=param_final of=vp_inv bs=278444 count=1 skip=0  
dd if=param_final of=rho_inv bs=278444 count=1 skip=1

psimage < vp_inv n1=151 n2=461 d1=0.02 d2=0.02  title="(a)" $option $colormap2 lbeg=1500 lend=5500 >vp_inv.ps
psimage < rho_inv n1=151 n2=461 d1=0.02 d2=0.02 title="(b)" $option $colormap2 lbeg=1000 lend=2700 > rho_inv.ps

psmerge in=vp_inv.ps translate=0,0 in=rho_inv.ps translate=9,0  > marm_inv.eps

rm *.ps
epstopdf marm_inv.eps

