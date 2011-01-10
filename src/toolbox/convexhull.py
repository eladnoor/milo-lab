#!/usr/bin/python

import random
import time
import numpy as n
import pylab as p


def _angle_to_point(point, centre):
    """calculate angle in 2-D between points and x axis"""
    delta = point - centre
    res = n.arctan(delta[1] / delta[0])
    if delta[0] < 0:
        res += n.pi
    return res


def _draw_triangle(p1, p2, p3, **kwargs):
    tmp = n.vstack((p1,p2,p3))
    x,y = [x[0] for x in zip(tmp.transpose())]
    p.fill(x,y, **kwargs)
    #time.sleep(0.2)


def area_of_triangle(p1, p2, p3):
    """calculate area of any triangle given co-ordinates of the corners"""
    return n.linalg.norm(n.cross((p2 - p1), (p3 - p1)))/2.


def QuickHull(points):
    """Randomized divide and conquer convex hull.
    
    Args:
        points: NxD matrix of points in dimension D.
    """
    N, D = points.shape
    dim = random.randint(0, D-1)
    min_dim = p.amin(points.T, dim)
    max_dim = p.amax(points.T, dim)
     

def convex_hull(points, graphic=True, smidgen=0.0075):
    """Calculate subset of points that make a convex hull around points

    Recursively eliminates points that lie inside two neighbouring points until only convex hull is remaining.
    
    :Parameters:
        points : ndarray (2 x m)
            array of points for which to find hull
        graphic : bool
            use pylab to show progress?
        smidgen : float
            offset for graphic number labels - useful values depend on your data range
    
    :Returns:
        hull_points : ndarray (2 x n)
            convex hull surrounding points
    """
    if graphic:
        p.clf()
        p.plot(points[0], points[1], 'ro')
    n_pts = points.shape[1]
    assert(n_pts > 5)
    centre = points.mean(1)
    if graphic: p.plot((centre[0],),(centre[1],),'bo')
    angles = n.apply_along_axis(_angle_to_point, 0, points, centre)
    pts_ord = points[:,angles.argsort()]
    if graphic:
        for i in xrange(n_pts):
            p.text(pts_ord[0,i] + smidgen, pts_ord[1,i] + smidgen, \
                   '%d' % i)
    pts = [x[0] for x in zip(pts_ord.transpose())]
    prev_pts = len(pts) + 1
    k = 0
    while prev_pts > n_pts:
        prev_pts = n_pts
        n_pts = len(pts)
        if graphic: p.gca().patches = []
        i = -2
        while i < (n_pts - 2):
            Aij = area_of_triangle(centre, pts[i],     pts[(i + 1) % n_pts])
            Ajk = area_of_triangle(centre, pts[(i + 1) % n_pts], \
                                   pts[(i + 2) % n_pts])
            Aik = area_of_triangle(centre, pts[i],     pts[(i + 2) % n_pts])
            if graphic:
                _draw_triangle(centre, pts[i], pts[(i + 1) % n_pts], \
                               facecolor='blue', alpha = 0.2)
                _draw_triangle(centre, pts[(i + 1) % n_pts], \
                               pts[(i + 2) % n_pts], \
                               facecolor='green', alpha = 0.2)
                _draw_triangle(centre, pts[i], pts[(i + 2) % n_pts], \
                               facecolor='red', alpha = 0.2)
            if Aij + Ajk < Aik:
                if graphic: p.plot((pts[i + 1][0],),(pts[i + 1][1],),'go')
                del pts[i+1]
            i += 1
            n_pts = len(pts)
        k += 1
    return n.asarray(pts)


def find_convex_hull(X, num_iter, num_points=None):
    """
        if num_points is set to None, find_convex_hull will return all the points in
        the convex hull (that have been found) sorted according to their sharpness.
        Otherwise, it will return the N-sharpest points.
    """
    (N, D) = X.shape
    if (num_points == None):
        num_points = N

    # randomly choose 'num_iter' direction on the unit sphere.
    # find the maximal point in the chosen direction, and add 1 to its counter.
    # only points on the convex hull will be hit, and 'sharp' corners will 
    # have more hits than 'smooth' corners.
    hits = p.zeros((N, 1))
    for j in xrange(num_iter):
        a = p.randn(D)
        a = a / p.norm(a)
        i = p.dot(X, a).argmax()
        hits[i] += 1
    
    # don't take points with 0 hits
    num_points = min(num_points, sum(p.find(hits)))
    
    # the indices of the n-best points
    o = list(p.argsort(hits, 0)[xrange(-1, -(num_points+1), -1)].flat)
    
    return X[o, :]

def sample_unit_triangle(N, D):
    X = p.ones((N, D))
    for i in xrange(N):
        while (sum(X[i, :]) >= 1):
            X[i, :] = p.random(D)
    return X

def test():
    X = sample_unit_triangle(100, 3)
    Z = find_convex_hull(X, 1000, 3)
    p.plot(X[:,0], X[:,1], 'y.')
    p.plot(Z[:,0], Z[:,1], 'r.')
    p.show()
    
if (__name__ == "__main__"):
    test()
