from itertools import product
import numpy as np
import argparse as arp
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import FuncAnimation

## Parser

parser = arp.ArgumentParser()
parser.add_argument("diff_sol")
args = parser.parse_args()

read_data = np.genfromtxt(args.diff_sol, delimiter=",")

R1 = read_data[0,3]
R2 = read_data[0,4]
R3 = read_data[0,5]

X1 = read_data[1:,0] ; Y1 = read_data[1:,1]
X2 = read_data[1:,2] ; Y2 = read_data[1:,3]
X3 = read_data[1:,4] ; Y3 = read_data[1:,5]

t = read_data[:,6]

#plt.scatter( [ X1[0],X2[0],X3[0] ] , [Y1[0],Y2[0],Y3[0]] , color = "k")

plt.plot(X1,Y1)
plt.plot(X2,Y2)
plt.plot(X3,Y3)

plt.show()

if input("Desea guardar la animacion? (y/n): ") == "y":

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
    

    M = ax1.transData.get_matrix()

    # Animating

    line1, = ax1.plot([], [], lw = 1.2, color = "gray", alpha = 0.2)
    line2, = ax1.plot([], [], lw = 1.2, color = "gray", alpha = 0.2)
    line3, = ax1.plot([], [], lw = 1.2, color = "gray", alpha = 0.2)

    pts1, = ax1.plot([], [], "o" , ls= " ", ms = 1.33*M[0][0]*R1, color= "k")
    pts2, = ax1.plot([], [], "o" , ls= " ", ms = 1.33*M[0][0]*R2, color= "k")
    pts3, = ax1.plot([], [], "o" , ls= " ", ms = 1.33*M[0][0]*R3, color= "k")

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
    
        return line1,line2,line3,pts1,pts2,pts3
    
    fig.canvas.draw()


    ani = FuncAnimation( fig, animate, init_func = init, frames = len(X1) , interval = 1, blit = True, repeat = False)

    Writer = animation.writers['ffmpeg']
    
    ## Change fps if want a more speeded animation
    writer = Writer(fps=90, metadata=dict(artist='Me'), bitrate=1800)
    ani.save('ani.mp4',writer = writer)

else: 
    pass




