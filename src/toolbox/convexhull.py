from pylab import *

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
    hits = zeros((N, 1))
    for j in xrange(num_iter):
        a = randn(D)
        a = a / norm(a)
        i = dot(X, a).argmax()
        hits[i] += 1
    
    # don't take points with 0 hits
    num_points = min(num_points, sum(find(hits)))
    
    # the indices of the n-best points
    o = list(argsort(hits, 0)[xrange(-1, -(num_points+1), -1)].flat)
    
    return X[o, :]

def sample_unit_triangle(N, D):
    X = ones((N, D))
    for i in xrange(N):
        while (sum(X[i, :]) >= 1):
            X[i, :] = random(D)
    return X
    
def test():
    X = sample_unit_triangle(100, 3)
    Z = find_convex_hull(X, 1000, 3)
    plot(X[:,0], X[:,1], 'y.')
    plot(Z[:,0], Z[:,1], 'r.')
    show()
    
if (__name__ == "__main__"):
    test()
