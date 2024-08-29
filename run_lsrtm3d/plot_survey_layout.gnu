set terminal postscript eps enhanced color font 'Helvetica,10' size 7.5cm,7.5cm
set output 'land_survey.eps'

#set terminal pngcairo size 800,800
# set xtics font "Helvetica,15"
# set ytics font "Helvetica,15"
# set key font "Helvetica,15"
# set title font 'Helvetica,20'
#set output '| display png:-'


#set multiplot layout 2,1
set grid back
set title "Acquisition geometry"
set xlabel 'X(m)'
set ylabel 'Y(m)'


plot "receivers.txt" using 2:3 with points ps 1 pt 1 title 'Receivers', \
"sources.txt" using 2:3 with points ps 1 pt 7 title 'Transmitters' 


