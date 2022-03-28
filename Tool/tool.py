import chart_studio.plotly as py
import plotly.graph_objs as go
import plotly.tools as tls
import plotly.express as px
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import emcee
import corner
import sys
import scipy as sp
import math as mt
from numpy import linalg as LA
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

def set_position(x,y,z,x0,y0,z0,alpha,beta):
    # transform -> rotate
    rotbeta1 = [np.cos(beta), 0, -np.sin(beta)]
    rotbeta2 = [0,1,0]
    rotbeta3 = [np.sin(beta), 0, np.cos(beta)]
    rotbeta  = np.array([rotbeta1,rotbeta2,rotbeta3])
    rotalpha1 = [1,0,0]
    rotalpha2 = [0,np.cos(alpha),np.sin(alpha)]
    rotalpha3 = [0,-np.sin(alpha),np.cos(alpha)]
    rotalpha  = np.array([rotalpha1,rotalpha2,rotalpha3]) # y-axis rotate
    data = np.array(([x-x0,y-y0,z-z0]))
    rot  = np.dot(rotalpha,rotbeta)
    pos  = np.dot(rot,data)
    return pos
def set_position_make_data(x,y,z,x0,y0,z0,alpha,beta):
    # rotate -> offset
    rotbeta1 = [np.cos(beta), 0, -np.sin(beta)]
    rotbeta2 = [0,1,0]
    rotbeta3 = [np.sin(beta), 0, np.cos(beta)]
    rotbeta  = np.array([rotbeta1,rotbeta2,rotbeta3])
    rotalpha1 = [1,0,0]
    rotalpha2 = [0,np.cos(alpha),np.sin(alpha)]
    rotalpha3 = [0,-np.sin(alpha),np.cos(alpha)]
    rotalpha  = np.array([rotalpha1,rotalpha2,rotalpha3]) # y-axis rotate
    rot  = np.dot(rotalpha,rotbeta)
    data = np.array([x,y,z])
    data = np.dot(rot,data)
    pos  = np.array([data[0]-x0,data[1]-y0,data[2]-z0]) 
    return pos
def set_position2(x,y,z,x0,y0,z0,vector):
    # transform -> rotater
    data = np.array(([x-x0,y-y0,z-z0]))
    rot  = rotM(vector)
    pos  = np.dot(rot,data)
    return pos
def rotM(p):
    # Three dimension rotation (radian)
    px = p[0]
    py = p[1]
    pz = p[2]
    Rx = np.array([[1, 0, 0],
                   [0, np.cos(px), np.sin(px)],
                   [0, -np.sin(px), np.cos(px)]])
    Ry = np.array([[np.cos(py), 0, -np.sin(py)],
                   [0, 1, 0],
                   [np.sin(py), 0, np.cos(py)]])
    Rz = np.array([[np.cos(pz), np.sin(pz), 0],
                   [-np.sin(pz), np.cos(pz), 0],
                   [0, 0, 1]])
    R = Rz.dot(Ry).dot(Rx)
    return R
def Parabola(x,y,z,f):
    return np.abs(np.abs(z) - 0.25*(x*x+y*y)/f) #return z_obs - z_model
def Parabola_(x,y,z,f):
    return (np.abs(z) - 0.25*(x*x+y*y)/f) #return z_obs - z_model
def Parabola2(x,y,f):
    return (x*x+y*y)/4/f #return x^2+y^2=4fz
def Cylinder(x,y,r):
    return r*r-x*x-y*y #return model - predict
def res_cylin(x,y,r):
    return r-np.sqrt(x*x+y*y) #return model - predict
def Cylinder2(x,y):
    return x*x+y*y #return model - predict
def plane_res(data,a,b,c):
    "ax+by+cz+d=0"
    "a/d*x+b/d*y+c/d*z+1=0"
    x,y,z = data
    residuals = a*x+b*y+c*z+1
    return residuals
def paraBolEqn(data,x0,y0,z0,f,alpha,beta): 
    x,y,z = data
    pos = set_position(x,y,z,x0,y0,z0,alpha,beta)
    return Parabola(pos[0],pos[1],pos[2],f)
def cylinderEqn(data,x0,y0,z0,r,alpha,beta):
    x,y,z = data
    pos = set_position(x,y,z,x0,y0,z0,alpha,beta)
    return Cylinder(pos[0],pos[1],r)
def distance_pl(oshift, tilt, point):
    AP=(oshift-point)
    PA=-AP
    e=tilt/np.sqrt(np.dot(tilt,tilt.T))
    d=np.sqrt(np.dot(AP,AP.T)-(np.dot(PA,e)**2))
    return d
def rms_distance(oshift, tilt, point):
    d = np.zeros(len(point))
    for i in range(len(point)):
        d[i] = distance_pl(oshift, tilt, point[i])
        print(i, distance_pl(oshift, tilt, point[i]))
    ave = np.sum(d)/len(d)
    sigma = np.sqrt(np.sum((d-ave)**2)/(len(d)-1))
    rms = np.sqrt(np.sum(d*d)/(len(d)))
    print ("mean of distance: ",ave)
    print ("sigma of distance: ",sigma)
def tangent_angle(u: np.ndarray, v: np.ndarray):
    i = np.inner(u, v)
    n = LA.norm(u) * LA.norm(v)
    c = i / n
    return np.rad2deg(np.arccos(np.clip(c, -1.0, 1.0)))
def rms_angle(a,b):
    angle = np.zeros(len(a))
    angle = angle*3600 #[arcsec]
    for i in range(len(angle)):
        angle[i] = tangent_angle(a[i],b[i])
        print(i, angle[i], " [arcsec]")
    ave = np.sum(angle)/len(angle)
    sigma = np.sqrt(np.sum((angle-ave)**2)/(len(angle)-1))
    #rms = np.sqrt(np.sum(angle*angle)/(len(angle)))
    print ("mean of angle: {:.3e}  [as]".format(ave))
    print ("sigma of angle: ",sigma)
def nor_vec(oshift, tilt, point):
    # Derive normal vector from the point to the line
    #vec = np.zeros(3)
    #s = -1*(tilt[0]*(oshift[0]-point[0])+tilt[1]*(oshift[1]-point[1])+tilt[2]*(oshift[2]-point[2]))/(tilt[0]**2+tilt[1]**2+tilt[2]**2)
    #vec[0] = oshift[0]-point[0] + s*tilt[0]
    #vec[1] = oshift[1]-point[1] + s*tilt[1]
    #vec[2] = oshift[2]-point[2] + s*tilt[2]
    #return np.array([vec[0],vec[1],vec[2]])
    vec = np.zeros(3)
    x0,y0,z0 = oshift #
    vx,vy,vz = tilt
    x1,y1,z1 = point 
    denomi = vx*(x1-x0)+vy*(y1-y0)+vz*(z1-z0)
    v2 = (vx*vx+vy*vy+vz*vz)
    t =  denomi/ v2
    vec[0] = x0+vx*t-x1
    vec[1] = y0+vy*t-y1
    vec[2] = z0+vz*t-z1
    return np.array([vec[0],vec[1],vec[2]])
def Rot_Matrix(A,B):
    # returning rotation matrix form A
    # A -> B
    # define the unit vector for each vector
    a = A / np.linalg.norm(A)
    b = B / np.linalg.norm(B)
    v = np.cross(a,b)
    c = np.dot(a,b)
    I = np.eye(3)
    vx = np.array([ [0,-v[2],v[1]], \
                    [v[2],0,-v[0]], \
                    [-v[1],v[0],0]
                  ])
    R = I + vx + np.dot(vx,vx)*(1/(1+c))
    return R
def define_plane_3points(point):
    #determine plane equation from 3 points
    import numpy as np
    p1,p2,p3 = point
    v1 = p1 - p3
    v2 = p2 - p3
    cp = np.cross(v1, v2)
    cp = cp/np.linalg.norm(cp)
    d = np.dot(cp, p3)
    a,b,c = cp
    print('The equation is {0}x + {1}y + {2}z = {3}'.format(a, b, c, d))
    return a,b,c,d
def intersection_plane_line(planeNormal, planePoint, rayDirection, rayPoint, epsilon=1e-6):
    ndotu = planeNormal.dot(rayDirection)
    if abs(ndotu) < epsilon:
        raise RuntimeError("no intersection or line is within plane")
    w = rayPoint - planePoint
    si = -planeNormal.dot(w) / ndotu
    Psi = w + si * rayDirection + planePoint
    #print(Psi)
    return Psi
def solve_3d_affine(p1, p2, p3, p4, s1, s2, s3, s4):
    x = np.transpose(np.matrix([p1,p2,p3,p4]))
    y = np.transpose(np.matrix([s1,s2,s3,s4]))
    # add ones on the bottom of x and y
    x = np.vstack((x,[1,1,1,1]))
    y = np.vstack((y,[1,1,1,1]))
    # solve for A2
    A2 = y * x.I
    #print(A2)
    # return function that takes input x and transforms it
    # don't need to return the 4th row as it is 
    return lambda x: (A2*np.vstack((np.matrix(x).reshape(3,1),1)))[0:3,:]
class estimate_ABvector():
    def __init__(self,norm_vector1,norm_point1,ref1,ref2):
        self.nvec = norm_vector1 
        self.n_p  = norm_point1   
        self.ref1 = ref1
        self.ref2 = ref2
        nvec = norm_vector1
        n_p  = norm_point1
        self.trans_function = solve_3d_affine(ref1[0],ref1[1],ref1[2],ref1[3],ref2[0],ref2[1],ref2[2],ref2[3])
        self.insec_ref1     = intersection_plane_line(nvec,ref1[0],nvec,n_p)
        self.insec_ref1to2  = np.array(self.trans_function(self.insec_ref1).T).ravel()
        ref1_vec = ref1[3]- ref1[0]
        ref2_vec = ref2[3] - ref2[0]
        self.rot_ref12 = Rot_Matrix(ref1_vec,ref2_vec)
        self.nvec_pos2 = np.dot(self.rot_ref12,nvec)
    def AB_info_on_pos2(self):
        nvec          = self.nvec
        nvec_pos2     = self.nvec_pos2
        insec_ref1    = self.insec_ref1
        insec_ref1to2 = self.insec_ref1to2
        print("Optical Axis Vector on Pos1       : ", nvec)
        print("Optical Axis Intersection on Pos1 : ", insec_ref1)
        print("Optical Axis Vector on Pos2       : ", nvec_pos2)
        print("Optical Axis Intersection on Pos2 : ", insec_ref1to2)
    def convert_position(self,pos):
        #print("A position on Pos1         :", pos)
        #print("Converted position on Pos2 :", self.trans_function(pos).T.ravel().flatten())
        return self.trans_function(pos).T.ravel().flatten()
def get_center_of_mass(array):
    print('center_of_mass: ',np.average(array, axis=0))
    return np.average(array, axis=0)
def grep_refs_mean(df,ref_num,obs_num):
    ref = np.zeros(3*ref_num).reshape(ref_num,3)
    for i in range(ref_num): 
        ref[i] = np.mean((df.x[i*5:i*5+5],df.y[i*5:i*5+5],df.z[i*5:i*5+5]),axis=1)
    return ref
def calc_R(ref_origin,ref_move):
    mean_origin = np.mean(ref_origin,axis=0)
    mean_move   = np.mean(ref_move,axis=0)
    sub_origin  = ref_origin - mean_origin 
    sub_move    = ref_move   - mean_move  
    R, rmsd = align.rotation_matrix(sub_move,sub_origin)
    print("RMS:" ,rmsd)
    return R
def my_convert(p_origin,p_move,data,R):
    mean_origin = np.mean(p_origin,axis=0)
    mean_move   = np.mean(p_move,axis=0)
    #print(data.dtype)
    #print(np.dtype(mean_move))
    #print("CM-diff:",mean_origin- mean_move)
    sub_data    = data - mean_move
    points = np.dot(R,sub_data.T).T + mean_origin
    return points 
