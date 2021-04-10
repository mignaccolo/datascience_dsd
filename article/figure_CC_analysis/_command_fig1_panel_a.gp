reset

# adjust range
set xr [-1:167]
set yr [-1:1]

# adjust tics
set ytics format ""
set ytics 0.25
set noxtics

# set grid on y-axis
set grid y front

# create color boxes corresponding to different instrument
set obj rect from 1, graph 0 to 33, graph 1 back fc "#DCDCDC" fs solid 0.3
set obj rect from 34, graph 0 to 69, graph 1 back fc "#A9A9A9" fs solid 0.3
set obj rect from 70, graph 0 to 101, graph 1 back fc "#AFEEEE" fs solid 0.3
set obj rect from 102, graph 0 to 166, graph 1 back fc "#D6D6FA" fs solid 0.3

# create legend
set label '(N,{/Symbol m})' at 150,-0.6 left
set arrow nohead from 135,-0.6 to 147.5,-0.6 lw 4 lc "#FF0000"
set label '(N,{/Symbol s})' at 150,-0.75 left
set arrow nohead from 135,-0.75 to 147.5,-0.75 lw 4 lc "#0000FF"
set label '(N,{/Symbol g})' at 150,-0.9 left
set arrow nohead from 135,-0.9 to 147.5,-0.9 lw 4 lc "#00CD00"
set label '(N,{/Symbol k})' at 115,-0.6 left
set arrow nohead from 100,-0.6 to 112.5,-0.6 lw 4 lc "#00BFFF"
set label '(N,{/Symbol h})' at 115,-0.75 left
set arrow nohead from 100,-0.75 to 112.5,-0.75 lw 4 lc "#FF00FF"

# put tics on y-axis
set label '-1' at -2,-1 right
set label '-0.75' at -2,-0.75 right
set label '-0.5' at -2,-0.5 right
set label '-0.25' at -2,-0.25 right
set label '0' at -2,0 right
set label '0.75' at -2,0.75 right
set label '0.5' at -2,0.5 right
set label '0.25' at -2,0.25 right
set label '1' at -2,1 right

#put panel label
set label '(a)' at 157,0.925

# plot data
p "crosscorr_groundpdf_matrixformat_allsites_Nvsmu" u 1:2 title "" w l lt 1 lw 2 lc "#FF0000","crosscorr_groundpdf_matrixformat_allsites_Nvssigma" u 1:2 title "" w l lt 1 lw 2 lc "#0000FF","crosscorr_groundpdf_matrixformat_allsites_Nvsgamma" u 1:2 title "" w l lt 1 lw 2 lc "#00CD00","crosscorr_groundpdf_matrixformat_allsites_Nvskappa" u 1:2 title "" w l lt 1 lw 2 lc "#00BFFF","crosscorr_groundpdf_matrixformat_allsites_Nvseta" u 1:2 title "" w l lt 1 lw 2 lc "#FF00FF",0 title "" w l lt 1 lc "#000000" lw 3



reset