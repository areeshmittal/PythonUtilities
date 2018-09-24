__author__ = 'Areesh Mittal'

'''Produces Time tables. I don't know the math behind it. I just know that it 
produces very nice looking figures.

Watch the video https://www.youtube.com/watch?v=qhbuKbxJsk8
for more details and algorithm used here.

See time_table_examples.pdf in the github folder for a few examples

Simply run the file to produce the image'''

import matplotlib.pyplot as plt
import scipy

def plot_time_table(Npoints,Fac,ax = None,lw = 1):
    
    thetas = scipy.arange(Npoints) * (2*scipy.pi) / Npoints
    points = scipy.array([scipy.cos(thetas),scipy.sin(thetas)])
    
    if ax is None:
        fig,ax = plt.subplots(1,1,figsize = (10,10))
    ax.set_aspect('equal')
    ax.axis('off')
    
    #draw boundary
    ax.plot(points[0,:],points[1,:],'k')
    ax.plot(points[0,[-1,0]],points[1,[-1,0]],'k')

    for idx1 in range(Npoints):#draw connecting lines
        idx2 = (Fac*idx1) % Npoints
    #    ax.plot(points[0,[idx1,idx2]],points[1,[idx1,idx2]],'k',lw = lw)
        ax.plot(*points[:,[idx1,idx2]],'k',lw = lw)
        
Npoints = 400
Fac = 73

plot_time_table(Npoints,Fac) 