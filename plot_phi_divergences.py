__author__ = 'Areesh Mittal'

import scipy
import matplotlib.pyplot as plt
import warnings
from scipy.spatial import ConvexHull

projection_matrix = scipy.array([[1,.5,0], [0,3**.5/2,0]])
uniform = scipy.ones(3)

def project(point):
    '''project point in 3d prob simplex to 2d triangle'''
    return projection_matrix @ point

def plot_grid_lines(ax):
    style, alpha = 'k--', 0.3
    
    def draw_line_between(p1,p2):#project p1 and p2 and connect by a line
        p1,p2 = project(p1),project(p2)
        ax.plot([p1[0],p2[0]],[p1[1],p2[1]], style, alpha = alpha)
        
    for i in range(1,10):
        prob = i/10
        draw_line_between([prob, 0, 1-prob],[prob, 1-prob, 0])
        draw_line_between([0, prob, 1-prob],[1-prob, prob, 0])    
        draw_line_between([1-prob, 0, prob],[0, 1-prob, prob])
            
def draw_prob_simplex(ax,draw_grid = True):
    '''removes axis lines and plots boundary of probability simplex'''
    ax.set_aspect('equal')
    ax.axis('off')
    boundary = [[0,.5,1,0],[0,3**.5/2,0,0]]
    ax.plot(*boundary,'k')

    ax.text(-.1,-.05,'(0,0,1)',fontsize = 20)
    ax.text(.9, -.05,'(1,0,0)',fontsize = 20)
    ax.text(.4,  .9, '(0,1,0)',fontsize = 20)
    
    if draw_grid == True:
        plot_grid_lines(ax)
        
def plot_phi_div(phi, rho, q = uniform, npoints = 50000, draw_grid = True):
    '''plot the region {p \in R^3 | D_phi(p,q) = sum(q_i phi(p_i/q_i) <= rho,
       sum(p) == 1 , p >= 0} by projecting it in 2 dimensions
    
       phi: must be a vectorized function; phi(1) = 0. Convexity desired, but
            not verified; unexpected things can happen if phi is not convex;
            eg for KL divergence, phi = lambda x: x*scipy.log(x)
            
       rho: real positive number; radius of the ball
       
       The baseline distribution 'q' must be a 1-dimensional list or array with
       3 positive elements. Elements of q are normalized so that they sum to 1;
       
       npoints: number of points initially in the simplex. After filtering
       for distance, number of points is smaller
       
       draw_grid: boolean; indicates whether to draw a triangular grid'''
    
    if not scipy.isclose(phi(1),0):
        warnings.warn('phi(1) is not equal to 0')
    
    #normalizing exponential rvs provide uniform samples from simplex
    p = -scipy.log(scipy.rand(npoints,3))
    p = p/p.sum(axis = 1,keepdims = True)
    
    q = scipy.array(q)# if q is a list
    q = q/q.sum() # if q doesn't sum to 1
    
    dist = (q * phi(p/q)).sum(axis = 1)
    p = p[dist <= rho,:]
    
    fig,ax = plt.subplots(1,1,figsize = (7,7))
    draw_prob_simplex(ax,draw_grid = draw_grid)
    
    points = p @ projection_matrix.T
#    ax.scatter(points[:,0],points[:,1],s = 1)
    
    hull = ConvexHull(points)
    ax.plot(points[hull.vertices,0],points[hull.vertices,1],'k')
    ax.plot(points[hull.vertices[[-1,0]],0],points[hull.vertices[[-1,0]],1],'k')
    
    proj_q = project(q)
    ax.scatter(proj_q[0],proj_q[1],c = 'r')
    return fig,ax

###### SOME POPULAR PHI-DIVERGENCES ###############
def        kl(x): return x*scipy.log(x)
def    rev_kl(x): return -scipy.log(x)
def     j_div(x): return (x-1)*scipy.log(x)
def variation(x): return abs(x-1)
def hellinger(x): return (scipy.sqrt(x)-1)**2
###################################################

fig,ax = plot_phi_div(kl, rho =  0.15,
                      q = [.2,.2,.6], npoints = 50000, draw_grid = True)
#fig.show()