reset
set termoption enhanced
#save_encoding = GPVAL_ENCODING
set encoding utf8
set log y  ## no log y for collisional damping
# set log x # for BremL and GcollL
set format y "%h"
#set format y "%0.e"
#set format y "10^{%L}"
#set border  lw 0.8 
# set log z
# set format z "%0.e"
# set view 79.0000, 342.000, 1, 1
# set contour
# set cntrparam bspline #cubicset tics 
# set cntrparam levels auto 6
# set cntrlabel format "%0.e"
# set yrange [3e-6:1]
# set yrange [1e-14:1]
# #set title "Time evolution of Langmuir spectrum"
# set key bottom right
# set key at 10,-160 # GcollL
# set key at 10,1.1e-4
# set key at 10,0.011
# set key noautotitle
# set label "(b)" at 0.08,2e-7
# set xrange [0.1:10.0] ## L waves range
# set xrange [0.:0.6]
#set yrange [0.:2.e-5]
# set yrange [3.e-6:3.e-5]
# set yrange [0:6.0e-5]
# set xrange [0:10]
# set yrange [0:0.012] ## BremL range
# set yrange [-2.2e-3:1.2e-2]
# set ytics "0",2e-3,"6e-3" right
# set yrange [-180:0] # GcollL range (new)
# set yrange [-210:0] ## GcollL range (old)
# set yrange [2.5e-5:1.2e-4] ## BremL/GcollL range
# set ytics "-5e-5",2.e-5,"-8e-6"
# set xlabel "log(q)" # offset 0,0.5
set xlabel "q"
# set ylabel "  {/Symbol g}@^{L_{coll}}_{q}" rotate by 0 offset 1
# set ylabel "  P@^{L}_{q}" rotate by 0 offset 0.5
# set ylabel "-{P@^{L}_{q}/{/Symbol g}@^{L_{coll}}_{q}}" rotate by 0 offset 3,0
set ylabel "{/Symbol e}@^{L}_{q}(0)" rotate by 0 offset 2.8 ## L waves
# set title "Collisional damping"
# set title "Electrostatic bremsstrahlung in the new approximation"
# set title "-P@^{L}_{q}/{/Symbol g}@^{L_{coll}}_{q}"
# set title "Eletrostatic bremsstrahlung"
# set title "Initial Langmuir spectrum (log q)"
# set title "Initial Langmuir spectrum (log {/Symbol e}@^{L}_{q}(0))"
#set xrange [0:0.6] ## L waves range
# set yrange [3.e-8:2.e-5]
set key spacing 1.5
# plot "GcollL1D.wt" w l lw 2 lc 7
# plot "GcollL1D.wt" w l t "{/Symbol g}@^{L_{coll}}_{q}  " lw 3 lc 7
# plot "BrGcL1D.wt" w l t "{/*0.9 -P@^{L}_{q}/{/Symbol g}@^{L_{coll}}_{q}}    " lw 2 lc 4
# plot "BrGcL1D.wt" w l lw 2 lc 4
# plot "BremL1D.wt" w l lw 2 lc "blue"
# plot "BremL1D.wt" w l t "{/*0.95 Electrostatic bremsstrahlung of L waves}" lw 4 lc "blue"
# plot "GcollL1D.wt" u 1:3 w l t "{/Symbol g}@^{L_{coll}}_{q}    * Geff" lw 3 lc 2
# replot "BremL1D.wt" w l t "P@^{L}_{q}" lw 3 lc 3
# plot "IL01D.wt" t "Only spontaneus and induced emissions" w l lw 2
# replot "IL1D.wt" t "With P@^{L}_{q} and {{/Symbol g}@^{L_{coll}}_{q}}    included" w l lw 2 lc 4
plot "IL01D.wt" t "When IL0 calculated in double precision is read" w l lw 2
replot "bk-IL01D.wt" t "When IL0 is calculated in single precision" w l lw 2 lc 3
# plot "IL01D.wt" w l t "{Initial state of L waves}" lw 3
# replot "IL1D_0100.wt" w l t "{/*0.7{/Symbol t}=0100}" lw 3
# replot "IL1D_0200.wt" w l t "{/*0.7{/Symbol t}=0200}" lw 3
# replot "IL1D_0300.wt" w l t "{/*0.7{/Symbol t}=0300}" lw 3
# replot "IL1D_0400.wt" w l t "{/*0.7{/Symbol t}=0400}" lw 3
# replot "IL1D_0500.wt" w l t "{/*0.7{/Symbol t}=0500}" lw 3
# replot "IL1D_2500.wt" w l t "{/*0.7{/Symbol t}=2500}" lw 4
# replot "IL1D_5000.wt" w l t "{/*0.7{/Symbol t}=5000}" lc "orange-red" lw 4
# replot "IL1D_10000.wt" w l t "{/*0.7{/Symbol t}=10000}" lw 4
# #replot "IL1D_20000.wt" w l t "{/*0.7{/Symbol t}=20000}" lw 4
# replot "../reb/IL1D_30000.wt" w l t "{/*0.7{/Symbol t}=30000}" lc "black" lw 4
#replot "../br-coll/IS1D_50000.wt" w l t "{/Symbol t}=50000" lc "black" lw 1.5 dt "-.---"
# set terminal pngcairo size 1280,896 enhanced rounded color lw 1.5 dl 2 fontscale 1.0 font 'TeX Gyre Pagella, 16'
# set terminal pngcairo size 1000,700 enhanced rounded color lw 1.5 dl 2 fontscale 2.1 font 'TeX Gyre Pagella, 10'
set term post fontfile "/usr/share/texmf/fonts/type1/public/cm-super/sfrm1000.pfb" "SFRM1000"
set terminal post enhanced color lw 2 round dl 1 fontscale 0.9 font 'SFRM1000' 24 size 30cm,25cm
set out "IL0-comparison-log.eps"
replot
pause -1

# set xlabel "u"
# set ylabel "    Fe(u)" rotate by 0 offset -01,0.5
# set xrange [-16:16]
# set yrange [1e-14:1]
# set key spacing 0.7
# plot "Fe01D.wt" w l t "{/*0.7 {/Symbol t}=0}" lw 3
# # replot "Fe1D_0100.wt" w l t "{/*0.6{/Symbol t}=0100}" lw 4
# # replot "Fe1D_0500.wt" w l t "{/*0.6{/Symbol t}=0500}" lw 4
# # replot "Fe1D_2500.wt" w l t "{/*0.6{/Symbol t}=2500}" lw 4
# replot "Fe1D_5000.wt" w l t "{/*0.7{/Symbol t}=5000}" lw 3
# replot "Fe1D_10000.wt" w l t "{/*0.7{/Symbol t}=10000}" lw 3
# replot "Fe1D_15000.wt" w l t "{/*0.7{/Symbol t}=15000}" lw 3
# replot "Fe1D_20000.wt" w l t "{/*0.7{/Symbol t}=20000}" lc "orange-red" lw 3
# replot "../reb/Fe1D_30000.wt" w l t "{/*0.7{/Symbol t}=30000}" lw 3
# #replot "../br-coll/IS1D_50000.wt" w l t "{/Symbol t}=50000" lc "black" lw 1.5 dt "-.---"
# # set terminal pngcairo size 1280,896 enhanced rounded color lw 1.5 dl 2 fontscale 1.0 font 'TeX Gyre Pagella, 16'
# # set terminal pngcairo size 1000,700 enhanced rounded color lw 1.5 dl 2 fontscale 2.1 font 'TeX Gyre Pagella, 10'
# set term post fontfile "/usr/share/texmf/fonts/type1/public/cm-super/sfrm1000.pfb" "SFRM1000"
# set terminal post enhanced color lw 2 dl 1 fontscale 1.2 font 'SFRM1000' 28 size 30cm,20cm
# set out "../../../brem_draft/PRL_MS/Fe1D_tevol-0-30k.eps"
# replot
# pause -1


# set yrange [1e-14:1]
# set label "(b)" at -10.5, 0.2
# plot "Fe01D.wt" w l t "{/Symbol t}=0" lc "red" lw 2.5 dt 1
# replot "Fe1D_5000.wt" w l t "{/Symbol t}=5000" lc "green" lw 2.5 dt "-"
# replot "Fe1D_10000.wt" w l t "{/Symbol t}=10000" lc "blue" lw 2.5 dt "..."
# replot "Fe1D_20000.wt" w l t "{/Symbol t}=20000" lc "magenta" lw 2.5 dt  "-."
# replot "../reb/Fe1D_30000.wt" w l t "{/Symbol t}=30000" lc "black" lw 2.5 dt "-.."
# # replot "../reb/Fe1D_45000.wt" w l t "{/Symbol t}=45000" lc "cyan" lw 2.5 dt "...-" 
# # replot "Fe1D_20000.wt" w l t "{/Symbol t}=20000" lc "cyan" lw 2.5 dt "...-"
# # replot  w l t "{/Symbol t}=30000" lc "yellow" lw 2.5 dt "--."
# # replot "../reb/Fe1D_40000.wt" w l t "{/Symbol t}=40000" lc "orange-red" lw 2.5 dt "--.."
# # set terminal pngcairo size 1600,1200 enhanced rounded color lw 2.0 dl 2 fontscale 1.0 font 'roman,24'
# #set terminal pngcairo size 1280,896 enhanced rounded color lw 1.5 dl 2 fontscale 2.1 font 'TeX Gyre Pagella, 10'
# #set terminal pngcairo size 1000,700 enhanced rounded color lw 1.5 dl 2 fontscale 2.0 font 'TeX Gyre Pagella, 10'
# set term post fontfile "/usr/share/texmf/fonts/type1/public/cm-super/sfrm1000.pfb" "SFRM1000"
# set terminal post enhanced rounded color lw 2 dl 1 fontscale 2.1 font 'SFRM1000' 12
# set out "../../../brem_draft/Fe1D_tevol-5k-30k.eps"
# replot

# pause -1