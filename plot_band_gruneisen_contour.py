#Plot phonon dispersion: 
#created by Prasad Matukumilli matukumilli@gmail.com
#This script read bands.gnuplot files (produced from phonopy band.yaml with 
#BAND_CONNECTION = .TRUE., and BAND_POINTS = 101)
#that are calculated separately at three different uc volumes (orig, plus, minus)
#NOTE: band.yaml and bands.gnuplot files (for 3:orig,plus,minus) 
#must be generated separately for each high-symmetry path
 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
plt.rcParams.update({'font.size': 30})
plt.rcParams.update({'legend.handlelength': 0.5})
plt.rcParams['lines.linewidth'] = 10

#conv factor
thz2cminv=33.356
conv=thz2cminv

#set yaxis limits
ymax=2*conv
ymin=-0.1*conv

#set mode Gruneisen parameter bounds (colorbar limits)
cmin=-10
cmax=20

#fucntion to draw line plot with colourred according to a parameter (here it is mode Gruneisem)
def plot_colourline(x,y,c):
    c = cm.jet((c-cmin)/(cmax-cmin))
    #c = cm.jet((c-np.min(c))/(np.max(c)-np.min(c)))
    ax = plt.gca()
    for i in np.arange(len(x)-1):
        ax.plot([x[i],x[i+1]], [y[i],y[i+1]], c=c[i])
    #c = plt.colorbar(fig)
    return

#set figure sizes
fig = plt.figure(figsize=(16,10))
ax1 = plt.subplot2grid((1, 4), (0, 0), colspan=3)

#set row numbers that covers the data in band.gnuplot file. (To plot the whole dispersion i.e, all high symmetry paths, set b= max(rows)
a=0
b=1616
####################High-symmetry path Gamma-X
pbands = np.loadtxt('../plus/1-GX/bands.gnuplot')
px=pbands[a:b,0]
py=pbands[a:b,1]

mbands = np.loadtxt('../minus/1-GX/bands.gnuplot')
mx=mbands[a:b,0]
my=mbands[a:b,1]

bands = np.loadtxt('../orig/1-GX/bands.gnuplot')
x=bands[a:b,0]
y=bands[a:b:,1]

x3 = np.ma.masked_where(x > 0.0296, x)
px3 = np.ma.masked_where(px > 0.029, px)
mx3 = np.ma.masked_where(mx > 0.0302, mx)

#print(bands[0:110,0])
qpts1 = [0.00000000, 0.02963920, 0.05423330, 0.08387250, 0.10846660, 0.17039090]#, 0.20003010, 0.22462420, 0.25426340, 0.27885750]		#BiSeI
qptsText = ['$\Gamma$', 'X', 'S', 'Y', '$\Gamma$', 'Z']#, 'U', 'R', 'T', 'Z']
qptsp = [0.00000000, 0.02905800, 0.05316990, 0.08222790, 0.10633980, 0.16704990, 0.19610790, 0.22021980, 0.24927780, 0.27338970]
qptsm = [0.00000000, 0.03024410, 0.05534010, 0.08558420, 0.11068020, 0.17386830, 0.20411230, 0.22920830, 0.25945240, 0.28454840]
g= 0.02963920

####################High-symmetry path X-S
pbands = np.loadtxt('../plus/2-XS/bands.gnuplot')
px2=pbands[a:b,0]
py2=pbands[a:b,1]

mbands = np.loadtxt('../minus/2-XS/bands.gnuplot')
mx2=mbands[a:b,0]
my2=mbands[a:b,1]

bands = np.loadtxt('../orig/2-XS/bands.gnuplot')
x2=bands[a:b,0]
y2=bands[a:b:,1]

x32 = np.ma.masked_where(x2 > 0.0245, x2)
px32 = np.ma.masked_where(px2 > 0.024, px2)
mx32 = np.ma.masked_where(mx2 > 0.025, mx2)
gx = 0.02459410
pgx = 0.02411190
mgx = 0.02509600 

####################High-symmetry path S-Y
pbands = np.loadtxt('../plus/3-SY/bands.gnuplot')
px3=pbands[a:b,0]
py3=pbands[a:b,1]

mbands = np.loadtxt('../minus/3-SY/bands.gnuplot')
mx3=mbands[a:b,0]
my3=mbands[a:b,1]

bands = np.loadtxt('../orig/3-SY/bands.gnuplot')
x3=bands[a:b,0]
y3=bands[a:b:,1]

x33 = np.ma.masked_where(x3 > 0.0296, x3)
px33 = np.ma.masked_where(px3 > 0.029, px3)
mx33 = np.ma.masked_where(mx3 > 0.0302, mx3)
gxs = 0.02963920
pgxs = 0.02905800
mgxs = 0.03024410 

####################High-symmetry path Y-Gamma
pbands = np.loadtxt('../plus/4-YG/bands.gnuplot')
px4=pbands[a:b,0]
py4=pbands[a:b,1]

mbands = np.loadtxt('../minus/4-YG/bands.gnuplot')
mx4=mbands[a:b,0]
my4=mbands[a:b,1]

bands = np.loadtxt('../orig/4-YG/bands.gnuplot')
x4=bands[a:b,0]
y4=bands[a:b:,1]

x34 = np.ma.masked_where(x4 > 0.0245, x4)
px34 = np.ma.masked_where(px4 > 0.0241, px4)
mx34 = np.ma.masked_where(mx4 > 0.02509, mx4)
gxsy = 0.02459410
pgxsy = 0.02411190
mgxsy = 0.02509600 

####################High-symmetry path Gamma-Z
pbands = np.loadtxt('../plus/5-GZ/bands.gnuplot')
px5=pbands[a:b,0]
py5=pbands[a:b,1]

mbands = np.loadtxt('../minus/5-GZ/bands.gnuplot')
mx5=mbands[a:b,0]
my5=mbands[a:b,1]

bands = np.loadtxt('../orig/5-GZ/bands.gnuplot')
x5=bands[a:b,0]
y5=bands[a:b:,1]

x35 = np.ma.masked_where(x5 > 0.0619, x5)
px35 = np.ma.masked_where(px5 > 0.0607, px5)
mx35 = np.ma.masked_where(mx5 > 0.0631, mx5)
gxsyg = 0.06192430
pgxsyg = 0.06071010
mgxsyg = 0.06318810 

a1, b1, a2, b2, a3, b3, a4, b4, a5, b5, a6, b6, a7, b7, a8, b8, a9, b9, a10, b10, a11, b11, a12, b12, a13, b13, a14, b14, a15, b15, a16, b16 = 0, 101, 101, 202, 202, 303, 303, 404, 404, 505, 505, 606, 606, 707, 707, 808, 808, 909, 909, 1010, 1010, 1111, 1111, 1212, 1212, 1313, 1313, 1414, 1414, 1515, 1515, 1616

nb = 16 	# number of branches to plot
#######################################################################################
##1
####################High-symmetry path Gamma-X
numa = (py-my)/y
deno = (2938.6792 - 2606.3328)/2769.1828   #(vol_p - vol_m)/vol_o
gruu = -1*numa/deno
gru = y*conv
#ax1  = fig.add_subplot(111)

ykey=0
for i in range(nb):
    i = i+1
    xkey= ykey #str("a"+str(i))
    ykey= 101 + xkey #str("b"+str(i))
    plot_colourline(x3[xkey:ykey], gru[xkey:ykey], gruu[xkey:ykey])

##2
####################High-symmetry path X-S
numa = (py2-my2)/y2
gruu2 = -1*numa/deno
gru2 = y2*conv

ykey=0
for i in range(nb):
    i = i+1
    xkey= ykey #str("a"+str(i))
    ykey= 101 + xkey #str("b"+str(i))
    plot_colourline(g+x32[xkey:ykey], gru2[xkey:ykey], gruu2[xkey:ykey])

##3
####################High-symmetry path S-Y
numa = (py3-my3)/y3
gruu3 = -1*numa/deno
gru3 = y3*conv

ykey=0
for i in range(nb):
    i = i+1
    xkey= ykey #str("a"+str(i))
    ykey= 101 + xkey #str("b"+str(i))
    plot_colourline(g+gx+x33[xkey:ykey], gru3[xkey:ykey], gruu3[xkey:ykey])

##4
####################High-symmetry path Y-Gamma
gru41=(py4[a4:b4]-my4[a2:b2])/y4[a4:b4]/deno	#LA
#pra=ax1.plot(g+gx+gxs+x34[a1:b1], gru41, marker='*', linestyle='-',  linewidth=2, color='blue')
gru42=(py4[a1:b1]-my4[a3:b3])/y4[a2:b2]/deno	#TA
#pra=ax1.plot(g+gx+gxs+x34[a1:b1], gru42, marker='*', linestyle='-',  linewidth=2, color='green')
gru43=(py4[a3:b3]-my4[a1:b1])/y4[a3:b3]/deno	#TA1
#pra=ax1.plot(g+gx+gxs+x34[a1:b1], gru43, marker='*', linestyle='-',  linewidth=2, color='red')
gru44=(py4[a2:b2]-my4[a4:b4])/y4[a1:1]/deno	#TO
#pra=ax1.plot(g+gx+gxs+x34[a1:b1], gru44, marker='*', linestyle='-',  linewidth=2, color='orange')

numa = (py4-my4)/y4
gruu4 = -1*numa/deno
gru4 = y4*conv

ykey=0
for i in range(nb):
    i = i+1
    xkey= ykey #str("a"+str(i))
    ykey= 101 + xkey #str("b"+str(i))
    plot_colourline(g+gx+gxs+x34[xkey:ykey], gru4[xkey:ykey], gruu4[xkey:ykey])

# sometimes the mode correspondence across volumes has to be confirmed through inspection of eigen mode displacemnts 
plot_colourline(g+gx+gxs+x34[a1:b1], gru4[a1:b1], gru44)#gruu4[a1:b1])
plot_colourline(g+gx+gxs+x34[a2:b2-1], gru4[a2:b2-1], gru42[0:100])# gruu4[a2:b2-1])
plot_colourline(g+gx+gxs+x34[a3:b3], gru4[a3:b3], gru43)#gruu4[a3:b3])
plot_colourline(g+gx+gxs+x34[a4:b4], gru4[a4:b4], gru41)#u4[a4:b4])

##5
####################High-symmetry path Gamma-Z
gru51=(py5[a2:b2]-my5[a1:b1])/y5[a1:b1]/deno	#
#pra=ax1.plot(g+gx+gxs+x34[a1:b1], gru41, marker='*', linestyle='-',  linewidth=2, color='blue')
gru52=(py5[a1:b1]-my5[a2:b2])/y5[a2:b2]/deno	#
#pra=ax1.plot(g+gx+gxs+x34[a1:b1], gru42, marker='*', linestyle='-',  linewidth=2, color='green')
gru53=(py5[a3:b3]-my5[a3:b3])/y5[a3:b3]/deno	#
#pra=ax1.plot(g+gx+gxs+x34[a1:b1], gru43, marker='*', linestyle='-',  linewidth=2, color='red')
gru54=(py5[a4:b4]-my5[a4:b4])/y5[a4:b4]/deno	#TO
#pra=ax1.plot(g+gx+gxs+x34[a1:b1], gru44, marker='*', linestyle='-',  linewidth=2, color='orange')

numa = (py5-my5)/y5
gruu5 = -1*numa/deno
gru5 = y5*conv

ykey=0
for i in range(nb):
    i = i+1
    xkey= ykey #str("a"+str(i))
    ykey= 101 + xkey #str("b"+str(i))
    plot_colourline(g+gx+gxs+gxsy+x35[xkey:ykey], gru5[xkey:ykey], gruu5[xkey:ykey])

# sometimes the mode correspondence across volumes has to be confirmed through inspection of eigen mode displacemnts 
plot_colourline(g+gx+gxs+gxsy+x35[a1:b1], gru5[a1:b1], gru51)#gruu4[a1:b1])
plot_colourline(g+gx+gxs+gxsy+x35[a2:b2], gru5[a2:b2], gru52[0:100])# gruu4[a2:b2-1])
plot_colourline(g+gx+gxs+gxsy+x35[a3+1:b3], gru5[a3+1:b3], gru53[1:101])#gruu4[a3:b3])
plot_colourline(g+gx+gxs+gxsy+x35[a4:b4], gru5[a4:b4], gru54)#u4[a4:b4])


ax1.legend(frameon=False, loc='upper right')
ax1.set_xlim([0,max(qpts1)])
ax1.set_ylim([0,60])
ax1.set_xticks(qpts1)
ax1.set_xticklabels(qptsText)
plt.rcParams['lines.linewidth'] = 2
for xc in qpts1:
    ax1.axvline(x=xc, color='grey', linestyle='-')
    
ax1.axhline(y=0, color='grey', linestyle='-')

ax1.set_ylabel('Frequency (cm$^{-1}$)')

plt.tight_layout()
#plt.margins(0.2)
##############################COLOR BAR
plt.subplots_adjust(wspace=0.05, hspace=0)
hh=plt.scatter(g+x32[a1:b1], gru2[a1:b1], gruu2[a1:b1], cmap='cm.jet')  #choose a colorbar cmap
cbar=plt.colorbar(hh)	#, label='$\gamma$', )
plt.set_cmap(cm.jet)
plt.clim(cmin,cmax)
cbar.ax.set_ylabel('$\gamma$', rotation=0)

##############################Save Plot
plt.savefig('ABC-Grubands.contour.png', dpi=300)
plt.show()
