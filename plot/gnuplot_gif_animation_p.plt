if(exist('n')==0 || n<0) n = n0
time =  sprintf('./result/m=%03d.png',n)
#ttl =  sprintf('t=%d{/Symbol D}t',n)

#set title ttl
set output time
set xr[0:40]
set yr[0:15]
#set cbrange[0.:0.01]
set size ratio -1
set pm3d interpolate 10,10
unset key
str = sprintf("../data/data_pm3d_%08d", n)
splot str with pm3d


n = n + dn
if (n < n1) reread
undefine n


system("convert -delay 10 -loop 0 ./*.png ./animation.gif")