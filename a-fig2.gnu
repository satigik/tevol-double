# 
reset
set encoding utf8
# set key left top
# set log z
# set format z "10^{%L}"
set format z "%h"
set view 60.0000, 30.000, 1, 1
# set log zcb
# set format cb "%0.g"
# set view map
# unset surface
# set view map scale 1
# set pm3d map interpolate 0,0
# set palette color
# set palette model RGB
# set pal defined (1 '#00008f', 8 '#0000ff', 24 '#00ffff', 40 '#ffff00', 56 '#ff0000', 64 '#800000')
# set xrange [0:16]
# set style line 11 lc rgb '#000000' lt 1
# set border lw 1.5
# set tics right
# set cblabel "Fe(u)" rotate by 0
set xrange [-0.6:0.6]
set yrange [-0.6:0.6] 
set xlabel "q_x"
set ylabel "q_z"
#set zlabel " {/Symbol e}@^{L}_{q}({/Symbol t})" offset -5

# set yrange [1.0e-16:0.5]
# set yrange [2.e-6:5]
# set grid lw 1.0
# set xtics "-16",4,"16" right
# set ytics "-16",4,"16" right
# set cbtics left # "1e-16",1e-1,"1"
# set size square 1,1
set contour
# set cntrparam bspline #cubic
# set cntrparam levels auto 6
# set cntrlabel format "10^{%L}"
# set cntrlabel onecolor
# set key at -18,18,1e-4 title 'Levels'
#set hidden3d
# set zrange [1e-6:1e-4]
# splot "IL0.wt" t "{/Symbol t}=0" w l lc rgb "light-red"
# set term post fontfile "/usr/share/texmf/fonts/type1/public/cm-super/sfrm1000.pfb" "SFRM1000"
# set terminal post enhanced color lw 0.9 round dl 1 fontscale 0.9 font 'SFRM1000' 24 size 30cm,20cm
# set out "IL0.eps"
# replot

# set zlabel " {/Symbol e}@^{L}_{q}({/Symbol t})"
set zlabel "  P@^{L}_{q}" offset -3 # rotate by 0
splot "BremL.wt" t "Electrostatic bremsstrahlung" w l lc rgb "royalblue"
set terminal pngcairo size 1600,1400 enhanced color fontscale 0.9 font 'Times-New-Roman,24'
set out "BremL.png"
replot

set zlabel "  {/Symbol g}@^{L_{coll}}_{q}" offset -2 # rotate by 0
splot "GcollL.wt" t "Collisional damping" w l lc rgb "orange-red"
set terminal pngcairo size 1600,1400 enhanced color fontscale 0.9 font 'Times-New-Roman,24'
set out "GcollL.png"
replot
pause -1 

# splot "IL_0100.wt" t "{/Symbol t}=100" w l lc rgb "light-red"
# set term post fontfile "/usr/share/texmf/fonts/type1/public/cm-super/sfrm1000.pfb" "SFRM1000"
# set terminal post enhanced color lw 0.9 round dl 1 fontscale 0.9 font 'SFRM1000' 24 size 30cm,20cm
# set out "IL_0100.eps"
# replot

# splot "IL_0200.wt" t "{/Symbol t}=200" w l lc rgb "light-red"
# set term post fontfile "/usr/share/texmf/fonts/type1/public/cm-super/sfrm1000.pfb" "SFRM1000"
# set terminal post enhanced color lw 0.9 round dl 1 fontscale 0.9 font 'SFRM1000' 24 size 30cm,20cm
# set out "IL_0200.eps"
# replot

# splot "IL_0300.wt" t "{/Symbol t}=300" w l lc rgb "light-red"
# set term post fontfile "/usr/share/texmf/fonts/type1/public/cm-super/sfrm1000.pfb" "SFRM1000"
# set terminal post enhanced color lw 0.9 round dl 1 fontscale 0.9 font 'SFRM1000' 24 size 30cm,20cm
# set out "IL_0300.eps"
# replot

# splot "IL_0400.wt" t "{/Symbol t}=400" w l lc rgb "light-red"
# set term post fontfile "/usr/share/texmf/fonts/type1/public/cm-super/sfrm1000.pfb" "SFRM1000"
# set terminal post enhanced color lw 0.9 round dl 1 fontscale 0.9 font 'SFRM1000' 24 size 30cm,20cm
# set out "IL_0400.eps"
# replot

# splot "IL_0500.wt" t "{/Symbol t}=500" w l lc rgb "light-red"
# set term post fontfile "/usr/share/texmf/fonts/type1/public/cm-super/sfrm1000.pfb" "SFRM1000"
# set terminal post enhanced color lw 0.9 round dl 1 fontscale 0.9 font 'SFRM1000' 24 size 30cm,20cm
# set out "IL_0500.eps"
# replot

# splot "IL_0600.wt" t "{/Symbol t}=600" w l lc rgb "light-red"
# set term post fontfile "/usr/share/texmf/fonts/type1/public/cm-super/sfrm1000.pfb" "SFRM1000"
# set terminal post enhanced color lw 0.9 round dl 1 fontscale 0.9 font 'SFRM1000' 24 size 30cm,20cm
# set out "IL_0600.eps"
# replot

# splot "IL_0700.wt" t "{/Symbol t}=700" w l lc rgb "light-red"
# set term post fontfile "/usr/share/texmf/fonts/type1/public/cm-super/sfrm1000.pfb" "SFRM1000"
# set terminal post enhanced color lw 0.9 round dl 1 fontscale 0.9 font 'SFRM1000' 24 size 30cm,20cm
# set out "IL_0700.eps"
# replot

# splot "IL_0800.wt" t "{/Symbol t}=800" w l lc rgb "light-red"
# set term post fontfile "/usr/share/texmf/fonts/type1/public/cm-super/sfrm1000.pfb" "SFRM1000"
# set terminal post enhanced color lw 0.9 round dl 1 fontscale 0.9 font 'SFRM1000' 24 size 30cm,20cm
# set out "IL_0800.eps"
# replot

# splot "IL_0900.wt" t "{/Symbol t}=900" w l lc rgb "light-red"
# set term post fontfile "/usr/share/texmf/fonts/type1/public/cm-super/sfrm1000.pfb" "SFRM1000"
# set terminal post enhanced color lw 0.9 round dl 1 fontscale 0.9 font 'SFRM1000' 24 size 30cm,20cm
# set out "IL_0900.eps"
# replot

# splot "IL_1000.wt" t "{/Symbol t}=1000" w l lc rgb "light-red"
# set term post fontfile "/usr/share/texmf/fonts/type1/public/cm-super/sfrm1000.pfb" "SFRM1000"
# set terminal post enhanced color lw 0.9 round dl 1 fontscale 0.9 font 'SFRM1000' 24 size 30cm,20cm
# set out "IL_1000.eps"
# replot

# pause -1

# # set title "Time evolution in {/Symbol t}=0100 with collisions and no collisions"
# # splot "../g5em3/IL_0100.wt" t "{/Symbol t}=0100, no coll" w l lc rgb "royalblue"
#  splot "IL_0100.wt" t "{/Symbol t}=0100" w l lc rgb "light-red"
#  set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
#  set out "IL_0100.png"
#  replot


# # set title "Time evolution in {/Symbol t}=2000 with collisions and no collisions"
# # splot "../g5em3/IL_2000.wt" t "{/Symbol t}=2000, no coll" w l lc rgb "royalblue"
#  splot "IL_2000.wt" t "{/Symbol t}=2000" w l lc rgb "light-red"
#  set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
#  set out "IL_2000.png"
#  replot

# # set title "Time evolution in {/Symbol t}=3000 with collisions and no collisions"
# # splot "../g5em3/IL_3000.wt" t "{/Symbol t}=3000, no coll" w l lc rgb "royalblue"
#  splot "IL_3000.wt" t "{/Symbol t}=3000" w l lc rgb "light-red"
#  set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
#  set out "IL_3000.png"
#  replot

# # set title "Time evolution in {/Symbol t}=4000 with collisions and no collisions"
# # splot "../g5em3/IL_4000.wt" t "{/Symbol t}=4000, no coll" w l lc rgb "royalblue"
#  splot "IL_4000.wt" t "{/Symbol t}=4000" w l lc rgb "light-red"
#  set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
#  set out "IL_4000.png"
#  replot

# # set title "Time evolution in {/Symbol t}=5000 with collisions and no collisions"
# # splot "../g5em3/IL_5000.wt" t "{/Symbol t}=5000, no coll" w l lc rgb "royalblue"
#  splot "IL_5000.wt" t "{/Symbol t}=5000" w l lc rgb "light-red"
#  set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
#  set out "IL_5000.png"
#  replot

# # set title "Time evolution in {/Symbol t}=6000 with collisions and no collisions"
# # splot "../g5em3/IL_6000.wt" t "{/Symbol t}=6000, no coll" w l lc rgb "royalblue"
#  splot "IL_6000.wt" t "{/Symbol t}=6000" w l lc rgb "light-red"
#  set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
#  set out "IL_6000.png"
#  replot

# # set title "Time evolution in {/Symbol t}=7000 with collisions and no collisions"
# # splot "../g5em3/IL_7000.wt" t "{/Symbol t}=7000, no coll" w l lc rgb "royalblue"
#  splot "IL_7000.wt" t "{/Symbol t}=7000" w l lc rgb "light-red"
#  set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
#  set out "IL_7000.png"
#  replot

# # set title "Time evolution in {/Symbol t}=8000 with collisions and no collisions"
# # splot "../g5em3/IL_8000.wt" t "{/Symbol t}=8000, no coll" w l lc rgb "royalblue"
#  splot "IL_8000.wt" t "{/Symbol t}=8000" w l lc rgb "light-red"
#  set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
#  set out "IL_8000.png"
#  replot

# # set title "Time evolution in {/Symbol t}=9000 with collisions and no collisions"
# # splot "../g5em3/IL_9000.wt" t "{/Symbol t}=9000, no coll" w l lc rgb "royalblue"
#  splot "IL_9000.wt" t "{/Symbol t}=9000" w l lc rgb "light-red"
#  set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
#  set out "IL_9000.png"
#  replot

# # set title "Time evolution in {/Symbol t}=10000 with collisions and no collisions"
# # splot "../g5em3/IL_10000.wt" t "{/Symbol t}=10000, no coll" w l lc rgb "royalblue"
#  splot "IL_10000.wt" t "{/Symbol t}=10000" w l lc rgb "light-red"
#  set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
#  set out "IL_10000.png"
#  replot

# # set title "Time evolution in {/Symbol t}=11000 with collisions and no collisions"
# # splot "../g5em3/IL_11000.wt" t "{/Symbol t}=11000, no coll" w l lc rgb "royalblue"
#  splot "IL_11000.wt" t "{/Symbol t}=11000" w l lc rgb "light-red"
#  set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
#  set out "IL_11000.png"
#  replot

# # set title "Time evolution in {/Symbol t}=12000 with collisions and no collisions"
# # splot "../g5em3/IL_12000.wt" t "{/Symbol t}=12000, no coll" w l lc rgb "royalblue"
#  splot "IL_12000.wt" t "{/Symbol t}=12000" w l lc rgb "light-red"
#  set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
#  set out "IL_12000.png"
#  replot

# # set title "Time evolution in {/Symbol t}=13000 with collisions and no collisions"
# # splot "../g5em3/IL_13000.wt" t "{/Symbol t}=13000, no coll" w l lc rgb "royalblue"
#  splot "IL_13000.wt" t "{/Symbol t}=13000" w l lc rgb "light-red"
#  set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
#  set out "IL_13000.png"
#  replot

# # set title "Time evolution in {/Symbol t}=14000 with collisions and no collisions"
# # splot "../g5em3/IL_14000.wt" t "{/Symbol t}=14000, no coll" w l lc rgb "royalblue"
#  splot "IL_14000.wt" t "{/Symbol t}=14000" w l lc rgb "light-red"
#  set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
#  set out "IL_14000.png"
#  replot

# # set title "Time evolution in {/Symbol t}=15000 with collisions and no collisions"
# # splot "../g5em3/IL_15000.wt" t "{/Symbol t}=15000, no coll" w l lc rgb "royalblue"
#  splot "IL_15000.wt" t "{/Symbol t}=15000" w l lc rgb "light-red"
#  set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
#  set out "IL_15000.png"
#  replot

# # set title "Time evolution in {/Symbol t}=16000 with collisions and no collisions"
# # splot "../g5em3/IL_16000.wt" t "{/Symbol t}=16000, no coll" w l lc rgb "royalblue"
# # replot "IL_16000.wt" t "{/Symbol t}=16000, coll" w l lc rgb "light-red"
# # set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
# # set out "comp_IL_16000.png"
# # replot

# # set title "Time evolution in {/Symbol t}=17000 with collisions and no collisions"
# # splot "../g5em3/IL_17000.wt" t "{/Symbol t}=17000, no coll" w l lc rgb "royalblue"
# # replot "IL_17000.wt" t "{/Symbol t}=17000, coll" w l lc rgb "light-red"
# # set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
# # set out "comp_IL_17000.png"
# # replot

# # set title "Time evolution in {/Symbol t}=18000 with collisions and no collisions"
# # splot "../g5em3/IL_18000.wt" t "{/Symbol t}=18000, no coll" w l lc rgb "royalblue"
# # replot "IL_18000.wt" t "{/Symbol t}=18000, coll" w l lc rgb "light-red"
# # set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
# # set out "comp_IL_18000.png"
# # replot

# # set title "Time evolution in {/Symbol t}=19000 with collisions and no collisions"
# # splot "../g5em3/IL_19000.wt" t "{/Symbol t}=19000, no coll" w l lc rgb "royalblue"
# # replot "IL_19000.wt" t "{/Symbol t}=19000, coll" w l lc rgb "light-red"
# # set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
# # set out "comp_IL_19000.png"
# # replot

# # set title "Time evolution in {/Symbol t}=20000 with collisions and no collisions"
# # splot "../g5em3/IL_20000.wt" t "{/Symbol t}=20000, no coll" w l lc rgb "royalblue"
#  splot "IL_20000.wt" t "{/Symbol t}=20000" w l lc rgb "light-red"
#  set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
#  set out "IL_20000.png"
#  replot

# # set title "Time evolution in {/Symbol t}=21000 with collisions and no collisions"
# # splot "../reb/IL_21000.wt" t "{/Symbol t}=21000 , no coll" w l lc rgb "royalblue"
# # replot "IL_21000.wt" t "{/Symbol t}=21000, coll" w l lc rgb "light-red"
# # set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
# # set out "comp_IL_21000.png"
# # replot

# # set title "Time evolution in {/Symbol t}=22000 with collisions and no collisions"
# # splot "../reb/IL_22000.wt" t "{/Symbol t}=22000, no coll" w l lc rgb "royalblue"
# # replot "IL_22000.wt" t "{/Symbol t}=22000, coll" w l lc rgb "light-red"
# # set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
# # set out "comp_IL_22000.png"
# # replot

# # set title "Time evolution in {/Symbol t}=23000 with collisions and no collisions"
# # splot "../reb/IL_23000.wt" t "{/Symbol t}=23000, no coll" w l lc rgb "royalblue"
# # replot "IL_23000.wt" t "{/Symbol t}=23000, coll" w l lc rgb "light-red"
# # set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
# # set out "comp_IL_23000.png"
# # replot

# # set title "Time evolution in {/Symbol t}=24000 with collisions and no collisions"
# # splot "../reb/IL_24000.wt" t "{/Symbol t}=24000, no coll" w l lc rgb "royalblue"
# # replot "IL_24000.wt" t "{/Symbol t}=24000, coll" w l lc rgb "light-red"
# # set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
# # set out "comp_IL_24000.png"
# # replot

# # set title "Time evolution in {/Symbol t}=25000 with collisions and no collisions"
# # splot "../reb/IL_25000.wt" t "{/Symbol t}=25000, no coll" w l lc rgb "royalblue"
#  splot "IL_25000.wt" t "{/Symbol t}=25000" w l lc rgb "light-red"
#  set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
#  set out "IL_25000.png"
#  replot

# # set title "Time evolution in {/Symbol t}=26000 with collisions and no collisions"
# # splot "../reb/IL_26000.wt" t "{/Symbol t}=26000, no coll" w l lc rgb "royalblue"
# # replot "IL_26000.wt" t "{/Symbol t}=26000, coll" w l lc rgb "light-red"
# # set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
# # set out "comp_IL_26000.png"
# # replot

# # set title "Time evolution in {/Symbol t}=27000 with collisions and no collisions"
# # splot "../reb/IL_27000.wt" t "{/Symbol t}=27000, no coll" w l lc rgb "royalblue"
# # replot "IL_27000.wt" t "{/Symbol t}=27000, coll" w l lc rgb "light-red"
# # set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
# # set out "comp_IL_27000.png"
# # replot

# # set title "Time evolution in {/Symbol t}=28000 with collisions and no collisions"
# # splot "../reb/IL_28000.wt" t "{/Symbol t}=28000, no coll" w l lc rgb "royalblue"
# # replot "IL_28000.wt" t "{/Symbol t}=28000, coll" w l lc rgb "light-red"
# # set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
# # set out "comp_IL_28000.png"
# # replot

# # set title "Time evolution in {/Symbol t}=29000 with collisions and no collisions"
# # splot "../reb/IL_29000.wt" t "{/Symbol t}=29000, no coll" w l lc rgb "royalblue"
# # replot "IL_29000.wt" t "{/Symbol t}=29000, coll" w l lc rgb "light-red"
# # set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
# # set out "comp_IL_29000.png"
# # replot

# # set title "Time evolution in {/Symbol t}=30000 with collisions and no collisions"
# # splot "../reb/IL_30000.wt" t "{/Symbol t}=30000, no coll" w l lc rgb "royalblue"
#  splot "IL_30000.wt" t "{/Symbol t}=30000" w l lc rgb "light-red"
#  set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
#  set out "IL_30000.png"
#  replot

# set title "Time evolution in {/Symbol t}=31000 with collisions and no collisions"
# splot "../reb/IL_31000.wt" t "{/Symbol t}=31000, no coll" w l lc rgb "royalblue"
# replot "IL_31000.wt" t "{/Symbol t}=31000, coll" w l lc rgb "light-red"
# set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
# set out "comp_IL_31000.png"
# replot

# set title "Time evolution in {/Symbol t}=32000 with collisions and no collisions"
# splot "../reb/IL_32000.wt" t "{/Symbol t}=32000, no coll" w l lc rgb "royalblue"
# replot "IL_32000.wt" t "{/Symbol t}=32000, coll" w l lc rgb "light-red"
# set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
# set out "comp_IL_32000.png"
# replot

# set title "Time evolution in {/Symbol t}=33000 with collisions and no collisions"
# splot "../reb/IL_33000.wt" t "{/Symbol t}=33000, no coll" w l lc rgb "royalblue"
# replot "IL_33000.wt" t "{/Symbol t}=33000, coll" w l lc rgb "light-red"
# set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
# set out "comp_IL_33000.png"
# replot

# set title "Time evolution in {/Symbol t}=34000 with collisions and no collisions"
# splot "../reb/IL_34000.wt" t "{/Symbol t}=34000, no coll" w l lc rgb "royalblue"
# replot "IL_34000.wt" t "{/Symbol t}=34000, coll" w l lc rgb "light-red"
# set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
# set out "comp_IL_34000.png"
# replot

# set title "Time evolution in {/Symbol t}=35000 with collisions and no collisions"
# splot "../reb/IL_35000.wt" t "{/Symbol t}=35000, no coll" w l lc rgb "royalblue"
# replot "IL_35000.wt" t "{/Symbol t}=35000, coll" w l lc rgb "light-red"
# set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
# set out "comp_IL_35000.png"
# replot

# set title "Time evolution in {/Symbol t}=36000 with collisions and no collisions"
# splot "../reb/IL_36000.wt" t "{/Symbol t}=36000, no coll" w l lc rgb "royalblue"
# replot "IL_36000.wt" t "{/Symbol t}=36000, coll" w l lc rgb "light-red"
# set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
# set out "comp_IL_36000.png"
# replot

# set title "Time evolution in {/Symbol t}=37000 with collisions and no collisions"
# splot "../reb/IL_37000.wt" t "{/Symbol t}=37000, no coll" w l lc rgb "royalblue"
# replot "IL_37000.wt" t "{/Symbol t}=37000, coll" w l lc rgb "light-red"
# set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
# set out "comp_IL_37000.png"
# replot

# set title "Time evolution in {/Symbol t}=38000 with collisions and no collisions"
# splot "../reb/IL_38000.wt" t "{/Symbol t}=38000, no coll" w l lc rgb "royalblue"
# replot "IL_38000.wt" t "{/Symbol t}=38000, coll" w l lc rgb "light-red"
# set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
# set out "comp_IL_38000.png"
# replot

# set title "Time evolution in {/Symbol t}=39000 with collisions and no collisions
# splot "../reb/IL_39000.wt" t "{/Symbol t}=39000, no coll" w l lc rgb "royalblue"
# replot "IL_39000.wt" t "{/Symbol t}=39000, coll" w l lc rgb "light-red"
# set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
# set out "comp_IL_39000.png"
# replot

# # set title "Time evolution in {/Symbol t}=40000 with collisions and no collisions"
# # splot "../reb/IL_40000.wt" t "{/Symbol t}=40000, no coll" w l lc rgb "royalblue"
#  splot "IL_40000.wt" t "{/Symbol t}=40000" w l lc rgb "light-red"
#  set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
#  set out "IL_40000.png"
#  replot

# # set title "Time evolution in {/Symbol t}=40000 with collisions and no collisions"
# # splot "../reb/IL_40000.wt" t "{/Symbol t}=40000, no coll" w l lc rgb "royalblue"
#  splot "IL_50000.wt" t "{/Symbol t}=50000" w l lc rgb "light-red"
#  set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
#  set out "IL_50000.png"
#  replot

##### 1D #####
# set title "Evolution of L waves at {/Symbol t}=7000,{/Symbol t}=10000 and {/Symbol t}=13000 - 1D "
# plot "../time-evol-ne/g5em3/IL1D_13000.wt" t "{/Symbol t}=13000" w l lc rgb "royalblue" lw 1.5 dt 2
# replot "../time-evol-ne/g5em3/IL1D_10000.wt" t "{/Symbol t}=10000" w l lc rgb "orange-red" lw 1.5 dt 1
# replot "../time-evol-ne/g5em3/IL1D_7000.wt" t "{/Symbol t}=7000" w l lc rgb "dark-yellow" lw 1.5 dt 4
# set terminal pngcairo size 1600,1400 enhanced color font 'Times-New-Roman,24'
# set out "../time-evol-ne/g5em3/Lw-evol1D.png"
# replot

# pause -1