from itertools import product
import numpy as np
import argparse as arp
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import FuncAnimation
import matplotlib

## Parser

parser = arp.ArgumentParser()
parser.add_argument("diff_sol")
args = parser.parse_args()

read_data = np.genfromtxt(args.diff_sol, delimiter=",")

## Masses and radii

M1 = read_data[0,0] ; R1 = read_data[0,3]
M2 = read_data[0,1] ; R2 = read_data[0,4]
M3 = read_data[0,2] ; R3 = read_data[0,5]

x_cm = read_data[0,6] ; y_cm = read_data[0,7]

## X,Y Positions

X1 = read_data[1:,0] ; Y1 = read_data[1:,1]
X2 = read_data[1:,2] ; Y2 = read_data[1:,3]
X3 = read_data[1:,4] ; Y3 = read_data[1:,5]


S3 = np.sqrt( (X1-X2)**2 + (Y1-Y2)**2 )
S2 = np.sqrt( (X1-X3)**2 + (Y1-Y3)**2 )
S1 = np.sqrt( (X3-X2)**2 + (Y3-Y2)**2 )

## X,Y Momenta

KX1 = read_data[1:,6] ; KY1 = read_data[1:,7]
KX2 = read_data[1:,8] ; KY2 = read_data[1:,9]
KX3 = read_data[1:,10] ; KY3 = read_data[1:,11]

# Time

t = read_data[1:,12]

# Energy plots

K1 = 0.5/M1*KX1**2 + 0.5/M1*KY1**2
K2 = 0.5/M2*KX2**2 + 0.5/M2*KY2**2
K3 = 0.5/M3*KX3**2 + 0.5/M3*KY3**2

G = 39.43279791722677
U= -G* ( M1*M2/S3 + M1*M3/S2 + M2*M3/S1)

# Preview plot

fig = plt.figure(figsize = [14,6])
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax1.scatter([X1[0],X2[0],X3[0]],[Y1[0],Y2[0],Y3[0]], marker = "x", s = 10, color = ["#FF5733","#C44C33","#873828"], label = "Starting points")

ax1.plot(X1,Y1,"--", alpha = 0.6, lw = 1.2, color = "#FF5733", label = "$m_1 = {} \odot$".format(M1))
ax1.plot(X2,Y2,"--", alpha = 0.6, lw = 1.2, color = "#C44C33", label = "$m_2 = {} \odot$".format(M2))
ax1.plot(X3,Y3,"--", alpha = 0.6, lw = 1.2, color = "#873828", label = "$m_3 = {} \odot$".format(M3))

ax1.set_title("Animation preview", fontsize = 13)
ax1.set_xlabel("X [AU]", fontsize = 10)
ax1.set_ylabel("Y [AU]", fontsize = 10)

ax1.legend()


ax2.plot(t,K1, color = "k", label = "Kinetic energy")
ax2.plot(t,K2, color = "k")
ax2.plot(t,K3, color = "k")
ax2.plot(t, U, label = "Potential energy")
ax2.plot(t,K1+K2+K3+U, color = "r", label = "Total energy")

ax2.set_xlabel("Time")
ax2.set_ylabel("Energy")

ax2.set_title("Energy plots")

ax2.legend()

plt.show()

if input("Desea guardar la animacion? (y/n): ") == "y":

    matplotlib.use("Agg")

    ## Setting the limits of the graph

    xmax = max( np.concatenate([X1,X2,X3]) ) + max(R1,R2,R3)
    xmin = min( np.concatenate([X1,X2,X3]) ) - max(R1,R2,R3)

    ymax = max( np.concatenate([Y1,Y2,Y3]) ) + max(R1,R2,R3)
    ymin = min( np.concatenate([Y1,Y2,Y3]) ) - max(R1,R2,R3)

    z = 0.3

    if xmax - xmin > ymax - ymin:
    
        mid_point = (ymax + ymin)/2
        cm , cM = xmin - z*(xmax-xmin) , xmax + z*(xmax-xmin)
        X = True
    
    else:
    
        mid_point = (xmax + xmin)/2
        cm , cM = ymin - z*(ymax-ymin) ,  ymax + z*(ymax-ymin)
        X = False
    
    ## Plotting

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    if X:
        ax1.set_xlim([cm, cM] )
        ax1.set_ylim([cm + mid_point, cM + mid_point] )
    
    else:
        ax1.set_xlim([cm + mid_point, cM + mid_point] )
        ax1.set_ylim([cm, cM] )
    
    ax1.get_xaxis().set_ticks([])
    ax1.get_yaxis().set_ticks([])

    M = ax1.transData.get_matrix()

    # Animating

    line1, = ax1.plot([], [], lw = 1.2, color = "#FF5733", alpha = 0.2)
    line2, = ax1.plot([], [], lw = 1.2, color = "#C44C33", alpha = 0.2)
    line3, = ax1.plot([], [], lw = 1.2, color = "#873828", alpha = 0.2)

    pts1, = ax1.plot([], [], "o" , ls= " ", ms = 1.33*M[0][0]*R1, color= "#FF5733", label = "$m_1 = {} \odot$".format(M1))
    pts2, = ax1.plot([], [], "o" , ls= " ", ms = 1.33*M[0][0]*R2, color= "#C44C33", label = "$m_2 = {} \odot$".format(M2))
    pts3, = ax1.plot([], [], "o" , ls= " ", ms = 1.33*M[0][0]*R3, color= "#873828", label = "$m_3 = {} \odot$".format(M3))
    
    info = ax1.text(0.02,0.95," ",transform = ax1.transAxes,fontsize = 8)
    ax1.set_title("Three body problem animation", fontsize = 10)
    
    legend = ax1.legend(loc = "upper right", fontsize = 8)
    for i in [0,1,2] : legend.legendHandles[i]._legmarker.set_markersize(6)

    
    def init():
    
        for l in [line1,line2,line3,pts1,pts2,pts3]:
            l.set_data([],[])
            
        return line1,line2,line3,pts1,pts2,pts3

    def animate(i):
    
        line1.set_data(X1[:i],Y1[:i])
        pts1.set_data(X1[i], Y1[i])
    
        line2.set_data(X2[:i],Y2[:i])
        pts2.set_data(X2[i], Y2[i])
    
        line3.set_data(X3[:i],Y3[:i])
        pts3.set_data(X3[i], Y3[i])
        
        info.set_text("{} years".format(t[i]))
    
        return line1,line2,line3,pts1,pts2,pts3,info
    
    fig.canvas.draw()


    ani = FuncAnimation( fig, animate, init_func = init, frames = len(X1) , interval = 1, blit = True, repeat = False)

    # Writer = animation.writers['ffmpeg']
    
    ## Change fps if want a more speeded animation
    # writer = Writer(fps=90, metadata=dict(artist='Me'), bitrate=1800)
    # ani.save('ani.mp4',writer = writer)
    ani.save('ani.gif', writer='imagemagick', fps=120)

else: 
    pass



