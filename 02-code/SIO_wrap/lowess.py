# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -- - 
# LOWESS : Locally WEighted Scatter plot Smoothing estimator with variable 
# bandwidth
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -- - 


import numpy as np
from scipy import linalg
from geopy import distance, Point

from .jlab_python import latlon2xy_tangent, xy2latlon_sphere

def kernelB(t):
    #the kernel is zero outside of the normalized badwidth 1 by construction
    K = np.asarray((1 - abs(t)**2)**2)
    K[abs(t)>=1] = 0
    return K

def kernelHH(t, h):
    #the kernel is zero outside of the normalized bandwidth 1 by construction
    K = np.asarray(1 - abs(t/h)**3)
    K[K<0] = 0
    K = (1/h)  * (70/81)*K**3
    return K 

def LatLonLocalWess(ti, loni, lati, p, b):
    """
    Lowess smoother: Robust locally weighted regression.
    The lowess function fits a nonparametric regression curve to a scatterplot.
    The arrays x and y contain an equal number of elements; each pair
    (x[i], y[i]) defines a data point in the scatterplot. The function returns
    the estimated (smooth) values of y.

    """
    # don't compute CI just yet - to be added later
    # input is already sorted by time, no need to add that for now
    #b = 2 # variable bandwidth
    #p = 1 # polynomial order
    #ti = time
    #lati, loni = lat, lon
    N = 3 # number of iterations
    pp = p + 2

    # data are already sorted by time
    n1 = ti.shape[0] - 1

    d = np.ones(ti.shape)

    res = np.zeros(ti.shape)
    lat_e = np.zeros(ti.shape)
    lon_e = np.zeros(ti.shape)
    res[:] = np.nan
    lat_e[:] = np.nan
    lon_e[:] = np.nan

    for m in range(1, N+1): 
        for k in range(len(ti)): 
            #print("iteration: %s, timestep %s" % (m, k+1))
            c = 0
            b2 = b

            # initialize q to no data point
            # initialize condition of matrix
            g = 0
            q = []

            while (len(q) < pp or g<1) :
                h = 0.5 * (ti[min(n1, k+b2)] - ti[max(0, k-b2)]) + 0*10/(24*60)
                w = kernelHH(ti - ti[k], h)
                w = w * d

                # need to make sure the problem is not underdetermined
                # pick only values where w != 0

                ti2 = ti[w!=0]
                loni2 = loni[w!=0]
                lati2 = lati[w!=0]
                w = w[w!=0]

                q = np.arange(len(w)) # only a dummy variable
                #P = len(w)

                xi2, yi2 = latlon2xy_tangent(lati2, loni2, lati[k], loni[k])

                X = np.zeros((len(ti2), p+1))
                z = ti2 - ti[k]
                for j in range(p+1):
                    X[:, j] = z**j

                W = np.diag(w, 0)
                Y = xi2 + 1j*yi2
                R = np.dot(np.dot(X.T, W), X)

                # condition for inversion (ratio of largest singular value to the smallest)
                cR = np.linalg.cond(R)

                if cR < 1/np.finfo(float).eps:
                    g = 1

                    A = np.dot(np.dot(np.linalg.inv(R), X.T), W)
                    beta = np.dot(A, Y)

                    lat_e[k], lon_e[k] = xy2latlon_sphere(beta[0].real, beta[0].imag, lati[k], loni[k])

                    if m<N:
                        res[k] = distance.distance(Point(lat_e[k], lon_e[k]), Point(lati[k], loni[k])).km

                co = (len(q) < pp) + 2 * (g < 1)

                if co == 1: #undetermined system
                    b2 += 1
                    c += 1
                elif co ==2: # cond(R)
                    b2 += 1
                    c +=1 
                elif co == 3: #cond(R) and undetermined system
                    b2 += 1
                    c += 1

                # limit on possible infinite while statement
                if c > 10:
                    break
        if m<N:
            M = np.nanmedian(abs(res))
            d = kernelB(res/(6*M))
            d[np.isnan(d)] = 1

    return lat_e, lon_e