set cblabel "p"
set palette defined ( -1 '#000030', 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')
#set palette defined ( 0 0 0 0, 0.1667 0 0 1, 0.5 0 1 0,\
     0.8333 1 0 0, 1 1 1 1 )
#set term gif animate
 #set terminal postscript eps size 600,480 dashed enhanced color font 'Roman,10'
set terminal pngcairo font "Times,14" enhanced size 800,450


set pm3d 
set pm3d map
set xlabel "{x/r_0}"
set ylabel "{z/r_0}"
#set autoscale z
#

n0 = 1
n1 = 500
dn = 1
load 'gnuplot_gif_animation_p.plt'
