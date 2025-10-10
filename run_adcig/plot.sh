#colormap
colormap1='wrgb=1.0,0,0 grgb=1.0,1.0,1.0 brgb=0,0,1.0' 
colormap2='wrgb=0,0,1 grgb=1,1,0 brgb=1,0,0' #used by Ludovic
colormap3='bhls=0.666666,.5,1  ghls=.333333,.5,1 whls=0,.5,1 bps=24 ' #used by Wei Zhou
colormap4='bhls=.2,.5,1  ghls=0,.5,1 whls=.8,.5,1 bps=24' #used by Wei Zhou

option="width=7 height=4 label1=Z(km) label2=X(km) d1num=1 d2num=3 legend=1 lheight=4 lwidth=0.3 lstyle=vertright "

farith < vp_true in2=vp_init op=div | farith op=log > rvp
farith < rho_true in2=rho_init op=div | farith op=log > rrho

psimage < vp_init n1=251 n2=767 d1=0.012 d2=0.012 title="(a)" $option $colormap2 lbeg=1500 lend=5500 >vp_init.ps
psimage < rho_init n1=251 n2=767 d1=0.012 d2=0.012 title="(b)" $option $colormap2 lbeg=1000 lend=2700 > rho_init.ps

psimage < rvp n1=251 n2=767 d1=0.012 d2=0.012 title="(c)" $option legend=1 lheight=4 lwidth=0.3 lstyle=vertright lbeg=-0.2 lend=0.2 >rvp.ps
psimage < rrho n1=251 n2=767 d1=0.012 d2=0.012 title="(d)" $option legend=1 lheight=4 lwidth=0.3 lstyle=vertright lbeg=-0.2 lend=0.2 > rrho.ps
psmerge in=vp_init.ps translate=0,0 in=rho_init.ps translate=9,0 \
	in=rvp.ps translate=0,-5.5 in=rrho.ps translate=9,-5.5  > vp_rho_smooth.eps

psimage < image_xcorr n1=251 n2=767 d1=0.012 d2=0.012 perc=99 title="(a)" $option legend=0 >image1.ps
psimage < image_normalized_xcorr n1=251 n2=767 d1=0.012 d2=0.012 perc=99 title="(b)" $option legend=0> image2.ps
psmerge in=image1.ps translate=0,0 in=image2.ps translate=9,0  > rtm.eps

rm *.ps
epstopdf vp_rho_smooth.eps
epstopdf rtm.eps

