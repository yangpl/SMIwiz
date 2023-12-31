#colormap
colormap1='wrgb=1.0,0,0 grgb=1.0,1.0,1.0 brgb=0,0,1.0' 
colormap2='wrgb=0,0,1 grgb=1,1,0 brgb=1,0,0' #used by Ludovic
colormap3='bhls=0.666666,.5,1  ghls=.333333,.5,1 whls=0,.5,1 bps=24 ' #used by Wei Zhou
colormap4='bhls=.2,.5,1  ghls=0,.5,1 whls=.8,.5,1 bps=24' #used by Wei Zhou

option="width=7 height=4 label1=Z(km) label2=X(km) d1num=1 d2num=3 legend=1 lheight=4 lwidth=0.3 lstyle=vertright "



psimage < vp_marm_true n1=151 n2=461 d1=0.02 d2=0.02  title="(a)" $option $colormap2 lbeg=1500 lend=5500 >vp_true.ps
psimage < rho_true n1=151 n2=461 d1=0.02 d2=0.02 title="(b)" $option $colormap2 lbeg=1000 lend=2700 > rho_true.ps

psimage < vp_marm_init n1=151 n2=461 d1=0.02 d2=0.02 title="(c)" $option $colormap2 lbeg=1500 lend=5500 >vp_init.ps
psimage < rho_init n1=151 n2=461 d1=0.02 d2=0.02 title="(d)" $option $colormap2 lbeg=1000 lend=2700> rho_init.ps

psmerge in=vp_true.ps translate=0,0 in=rho_true.ps translate=9,0 \
	in=vp_init.ps translate=0,-5.5 in=rho_init.ps translate=9,-5.5 > marmousi.eps

rm *.ps
epstopdf marmousi.eps

