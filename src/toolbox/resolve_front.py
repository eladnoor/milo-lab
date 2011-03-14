#!/usr/bin/python

import numpy as n
import pylab as p
from scipy.cluster.hierarchy import distance 


def SampleUnitTriangle(N, D):
    """Randomly sample from the unit triangle.
    
    Args:
        N: the number of points to sample.
        D: the dimensionality.
        
    Returns:
        A matrix containing the sampled points.
    
    Stolen from Elad's code. I'm a thief, I know.
    """
    X = p.ones((N, D))
    for i in xrange(N):
        while (sum(X[i, :]) >= 1):
            X[i, :] = p.random(D)
    return X


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


def ConvexHull(points, graphic=False, smidgen=0.0075):
    """Calculate subset of points that make a convex hull around points

    Recursively eliminates points that lie inside two neighboring points
    until only convex hull is remaining. Only for 2d.
    
    NOTE(flamholz): shamelessly stolen from the SciPy cookbook.
        http://www.scipy.org/Cookbook/Finding_Convex_Hull
        
    NOTE(flamholz): multi-dimensional hulls can be computed using the QHull
    algorithm. I haven't found a high-dimensional implementation for Python,
    but a C/C++ implementation can be found at http://www.qhull.org/
    
    Args:
        points: ndarray (2 x m)
            array of points for which to find hull
        graphic: bool
            use pylab to show progress?
        smidgen: float
            offset for graphic number labels - useful values depend on your data range
    
    Returns:
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


def StdDist(points):
    """Returns the standard deviation of the pairwise distances.
    
    Args:
        points: an Nx2 matrix of points.
    """
    ds = distance.pdist(points)
    return p.mean(ds)


def ResolveFront(points, n_archetypes=5, n_iters=100):
    """Attempt to resolve the corners of the Pareto front.
    
    Attempts to resolve the Pareto front with the following
    randomized algorithm:
      - Find the convex hull
      - For each iteration
        - Add normally distributed noise to each hull point.
        - Compute the hull of the hull.
        - Record which points remain on the hull.
      - Return the K points which were most stable on the hull
        under normally distributed hull noise.
    
    NOTE(flamholz): this approach doesn't work if you add noise to 
    all the points rather than only those on the hull. The problem 
    is that the archetypes (because of their corneryness) are more 
    susceptible to being displaced from the hull the hull by interior
    points than points along linear/curve boundaries. Also, this tendency
    doesn't appear to be strong enough to use as a negative signal...
    
    Args:
        points: the points to find the front of.
        n_archetypes: the number of candidate archetypes to return.
        n_iters: the number of randomized iterations to return.
        
    Returns:
        A Kx2 matrix of candidate archetypes where K is the number
        of candidate archetypes requested.
    """
    # Used to scale the additive noise.
    std_dist = StdDist(points)
    
    # Original hull.
    hull = ConvexHull(points.T)
    hull_counts = dict((tuple(p), 0) for p in hull)
    for _ in xrange(n_iters):
        curr_hull = n.array(hull)  # Important! A copy!
        
        noise = n.random.normal(scale=std_dist,
                                size=n.shape(hull))
        curr_hull += noise
        # Map back to original points (before additive noise).
        key_map = dict((tuple(noisy_point), tuple(orig_point))
                       for orig_point, noisy_point in zip(hull, curr_hull))

        # Take the hull of the noisy hull.
        for point in ConvexHull(curr_hull.T):
            point_key = key_map[tuple(point)]
            hull_counts[point_key] += 1
    
    counts = [(count, point) for point, count in hull_counts.iteritems()]
    counts.sort(reverse=True)
    top_points = [point for _, point in counts[:n_archetypes]]
    return n.array(top_points)


def main():
    # Sample 100 points in the unit triangle.
    # Try to resolve their front. Plot them.
    points = SampleUnitTriangle(100, 2)
    front = ResolveFront(points)
        
    xpoints, ypoints = points.T
    xfront, yfront = front.T
    p.scatter(xpoints, ypoints, c='g')
    p.scatter(xfront, yfront, c='r')
    p.legend(('Original points', 'Candidate archetypes'))
    p.show()
    

if __name__ == '__main__':
    main()
