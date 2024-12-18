#colormap
colormap1='wrgb=1.0,0,0 grgb=1.0,1.0,1.0 brgb=0,0,1.0' 
colormap2='wrgb=0,0,1 grgb=1,1,0 brgb=1,0,0' #used by Ludovic
colormap3='bhls=0.666666,.5,1  ghls=.333333,.5,1 whls=0,.5,1 bps=24 ' #used by Wei Zhou
colormap4='bhls=.2,.5,1  ghls=0,.5,1 whls=.8,.5,1 bps=24' #used by Wei Zhou

option="width=7 height=4 label1=Z(km) label2=X(km) d1num=1 d2num=3 legend=1 lheight=4 lwidth=0.3 lstyle=vertright "

psimage < decon_psf_1par n1=251 n2=767 d1=0.012 d2=0.012 title="(a) PSF-based mig-decon" $option lbeg=-0.2 lend=0.2 >decon_psf.eps
psimage < decon_fft_1par n1=251 n2=767 d1=0.012 d2=0.012 title="(b) FFT-based mig-decon" $option lbeg=-0.2 lend=0.2 >decon_fft.eps


psmerge in=decon_psf.eps translate=0,0 in=decon_fft.eps translate=9,0  > decon_1par.eps

epstopdf decon_1par.eps
rm *.eps
