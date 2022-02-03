"""
Some Jlab functions trabslated to python

"""

import numpy as np
import cmath
import math
from numpy.lib.scimath import sqrt as csqrt
import gsw

"""
# test data

import numpy as np
from numpy import loadtxt
import matplotlib.pyplot as plt

a = loadtxt('timelatlon.txt')
lat = a[:, 1]
lon = a[:, 2]
time = a[:, 0]

lato = lat[0]
lono = lon[0]

"""

# Radius of the Earth (km)
R = 6371

def rot(x):
    """
    Complex-valued rotation
    rot(x) = exp(j*x)

    """
    # dimension manipulation
    x = np.asarray(x)
    if x.size == 1:
        x = x[np.newaxis]

    # put x between -pi and pi 
    x = x % (2*np.pi) # remainder
    x[x>np.pi] = x[x>np.pi] - 2 * np.pi

    x = np.asarray(np.exp(x * 1j))

    # cases
    x[x==np.pi/2] = 1j

    x[x==-np.pi/2] = -1j

    x[np.logical_or(x==np.pi, x==-np.pi)] = -1

    x[x==0] = 1
    x[x==2*np.pi] = 1
    x[x==-2*np.pi] = 1

    return x



def sind(alpha):
    return np.sin(np.deg2rad(alpha))

def cosd(alpha):
    return np.cos(np.deg2rad(alpha))

def jrad2deg(*alpha):
    """
    Convert radians to degrees.
    Output angles are in the range [-180,180]
    """
    c = 2 * np.pi / 360.

    theta = []
    for var in alpha:
        theta_i = np.angle(np.exp(var * 1j)) / c
        theta.append(theta_i)

    theta = np.asarray(theta)
    theta[theta==180] = -180

    return theta


def jdeg2rad(*alpha):
    """
    Convert degrees to radians.
    Output angles are in the range [-pi,pi]
    """
    c = 2 * np.pi / 360.

    theta = []
    for var in alpha:
        theta_i = np.angle(np.exp(var * 1j * c))
        theta.append(theta_i)

    theta = np.asarray(theta)
    theta[theta==np.pi] = -np.pi

    return theta    

def imlog(x):
    """
    Imaginary part of log

    """
    return np.imag(np.log(x))


def latlon2xy_tangent(lat, lon, lato, lono):
    """
    returns x, y with units of km
    Converts lat, lon with units of degrees into displacements (x, y) in a 
    plane tangent to the earth at the point (lato, lono)
    """

    lon = lon - lono

    coslat=cosd(lat)
    coslon=cosd(lon)
    coslato=cosd(lato)
    sinlato=sind(lato)
    sinlat=sind(lat)

    # dot product
    dotp = coslat * coslon * coslato + sinlato * sinlat

    cond = dotp > 0

    # initialise arrays; 
    # change format so this can work if lat lon are size 1
    x = np.ones(len(lat))
    y = np.ones(len(lat))
    x[:] = y[:] = np.nan

    if len(dotp) > 0:
        x[cond] = R * coslat[cond] * sind(lon[cond])
        y[cond] = ((-1) * R * coslat[cond] * sinlato * coslon[cond]
                        + R * coslato * sinlat[cond])

    return x, y



def xy2latlon_sphere(x, y, lato, lono):
    """
    Convert (x,y) position with units of km, specifying a location in a plane 
    tangent to the earth at the point (LATO,LONO), into 
    latitude and longitude locations

    """
    # convert lato, lono to rad
    lato, lono = jdeg2rad(lato, lono)

    # distance in tangent plane
    r = csqrt(x**2 + y**2)

    x = np.asarray(x)
    y = np.asarray(y)

    x[r>=R] = np.nan
    y[r>=R] = np.nan

    r1 = R - csqrt(R**2 - r**2)

    # r1 = tangential distance from tangent plane to surface of earth
    #  choosing smaller root, corresponding to near side of earth
    
    #Now choose an xyz coordinate system with x=east, z= north
    
    #tangent point = point on sphere at which plane is tangent
    #contact point = point on sphere directly underneath (x,y) point in plane

    R1 = csqrt((R-r1)**2 + y**2)
    #R1 = distance from center of earth to contact point 
    # projected onto the xz plane, i.e., looking down the y-axis
    
    gamma=np.arcsin(y/R1)
    #gamma = angle (rad) spanned between contact point and tangent point
    #projected onto the xz plane, i.e., looking down the y-axis
    
    phi = lato + gamma
    #gamma = angle spanned between contact point and x-axis
    #projected onto the xz plane, i.e., looking down the y-axis
    
    xo = R1 * np.cos(phi)
    zo = R1 * np.sin(phi)

    yo = csqrt(R**2 - xo**2 - zo**2)

    if (x<0).size>1:
        yo[x<0] = -yo[x<0]


    R1 = csqrt(abs(xo)**2 + abs(yo)**2 + abs(zo)**2)
    x1 = xo/R1
    y1 = yo/R1
    z1 = zo/R1
    
    phi = np.arcsin(z1)
    th = imlog(x1 + 1j*y1)
    
    [lat,lon] = jrad2deg(phi, th)

    lon=jrad2deg(np.angle(rot(jdeg2rad(lon)+lono)))

    return lat, lon




#------------------------------------------------------------------------


def latlon2uv_forward(time, lat, lon):
    """
    time: float, in seconds

    Forward velocity (m/s)
    """
    dt = np.diff(time)

    # shift values forward by one (first value goes to the end)
    lat2 = lat[1:]
    lat2 = np.hstack((lat2, lat[0]))

    lon2 = lon[1:]
    lon2 = np.hstack((lon2, lon[0]))

    # distance between successive points (metres)
    dr = gsw.distance(lon, lat)

    xx = sind(lon2 - lon) * cosd(lat2)
    yy = (cosd(lat) * sind(lat2) 
         - sind(lat) * cosd(lat2) * cosd(lon2-lon))
    gamma = np.arctan2(yy, xx)

    # match size of arrays to lat/lon/time by repeating last point
    dt = np.hstack((dt.flatten(), dt[-1]))
    dr = np.hstack((dr.flatten(), dr[-1]))
    gamma[-1] = gamma[-2]

    return dr, dt, gamma

def latlon2uv_celloop_one(time, lat, lon, str):

    if str=='forward':
        dr, dt, gamma = latlon2uv_forward(time, lat, lon)
    elif str == 'backward':
        dr, dt, gamma = latlon2uv_forward(np.flipud(time), np.flipud(lat), np.flipud(lon))
        dr = np.flipud(dr)
        dt = np.flipud(dt)
        gamma = np.flipud(gamma)

    u = dr / dt * np.cos(gamma)
    v = dr / dt * np.sin(gamma)
    cond = np.logical_or(np.isnan(u), np.isnan(v))
    u[cond] = np.nan
    v[cond] = np.nan

    return u, v

def latlon2uv(time, lat, lon):

    u_forw, v_forw = latlon2uv_celloop_one(time, lat, lon, 'forward')
    u_back, v_back = latlon2uv_celloop_one(time, lat, lon, 'backward')

    u = np.zeros(len(lat))
    u[0] = u_forw[0]
    u[-1] = u_back[-1]
    u[1:-1] = 0.5*(u_forw[1:-1] + u_back[1:-1])

    v = np.zeros(len(lat))
    v[0] = v_forw[0]
    v[-1] = v_back[-1]
    v[1:-1] = 0.5*(v_forw[1:-1] + v_back[1:-1])

    return u, v

def latlon2uv_forward_mine(time, lat, lon):

    u_forw, v_forw = latlon2uv_celloop_one(time, lat, lon, 'forward')

#    u = np.zeros(len(lat))
#    u[0] = u_forw[0]
#    u[-1] = u_back[-1]
#    u[1:-1] = 0.5*(u_forw[1:-1] + u_back[1:-1])

#    v = np.zeros(len(lat))
#    v[0] = v_forw[0]
#    v[-1] = v_back[-1]
#    v[1:-1] = 0.5*(v_forw[1:-1] + v_back[1:-1])

    return u_forw, v_forw


