#get Viking Graben data from the web
wget https://s3.amazonaws.com/open.source.geoscience/open_data/Mobil_Avo_Viking_Graben_Line_12/seismic.segy

#clean up all non-segy data
rm *.su *.bin *.txt

#convert segy to su
segyread tape=seismic.segy hfile=header.txt bfile=binary.bin | segyclean> seismic.su

#kill traces with amplitude=0


#3D (scales with 1/r) to 2D (scales with 1/sqrt(r)) conversion, multiply by sqrt(r) (equivalent to multply by sqrt(t))
sugain < seismic.su tpow=0.5 > seismic2d.su

#select shots every 500 m, surange < seismic_converted.su will give you an idea on the range of sx
suwind < seismic2d.su key=sx s=3237 j=250 > seismic_subsampled.su

#we use reciprocotiy to switch source and receiver locations:
#get headers associated with sx,sy,selev,gx,gy,gelev and print in txt file as 6 columns with order: gx,gy,gelev,sx,sy,-selev
#we expect gy=0,sy=0 corresponding to the 2nd and the 5th column, gelev=-10, selev=6
# sugethw <seismic_clean.su output=geom key=sx,sy,selev,gx,gy,gelev | awk '{print $4,$5,$6,$1,$2,-$3}'> sorted_header.txt
# a2b < sorted_header.txt n1=6 >sorted_header.bin
# sushw <seismic_clean.su infile=sorted_header.bin key=sx,sy,selev,gx,gy,gelev >seismic_corrected.su

#split data into shots using keyword=fldr/ep, named as shot_xxxx.su (not necessarily ensure we have shot_0001.su, it may be killed out)
#fldr = Field Record Number (usually shot number), ep = Energy Source Point number (alternative shot number)
susplit < seismic_subsampled.su key=fldr numlength=4 stem=shot middle=_

#We also need to know the frequency range
#suspecfx < shot_0003.su | suximage perc=99
#fmin ~= 7Hz; fmax ~= 60Hz

k=1
for i in `ls shot_*.su`; do
    echo $i
    j=`printf "%04d" $k`
    #minimum phase bandpass filtering, re-numbering from 0001 to the last
    subfilt <$i fstoplo=5 fpasslo=10 fpasshi=15 fstophi=20 zerophase=0 >dat_$j.su 
    k=$[$k + 1] #important to have space around +
done
rm shot_*.su
