#!/usr/bin/python

from toolbox import convexhull
import pylab as p
import numpy as n
from scipy.cluster.hierarchy import distance 


def _LayoutPointSet(pset):
    x_coords = p.array([pt[0] for pt in pset])
    y_coords = p.array([pt[1] for pt in pset])
    return p.array([x_coords, y_coords])


def OnionSkin(points, layers=None):
    """Find the n outer layers of points.
    
    Args:
        points: a list of vectors.
    
    Returns:
        A 2 tuple (skin points, remaining points).
    """
    layers = layers or int(p.log(len(points[0]) + 1))
    
    coord_points = zip(points[0], points[1])
    all_points = set(coord_points)
    skin_points = set()
    
    for _ in xrange(layers):
        current_points_x = [pt[0] for pt in all_points]
        current_points_y = [pt[1] for pt in all_points]
        current_points = p.array([current_points_x, current_points_y])
        hull = convexhull.convex_hull(current_points, graphic=False)
        for point in hull:
            pt = tuple(point)
            all_points.discard(pt)
            skin_points.add(pt)
    
    return skin_points, all_points


def _CalcDensities(hull_points, all_points):
    ds = distance.pdist(list(all_points))
    std_d = p.std(ds)
    
    square_ds = distance.squareform(ds)
    densities = {}
    for i, point in enumerate(all_points):
        if point not in hull_points:
            continue
        
        my_ds = square_ds[i]
        density = len([1 for i in my_ds if i <= std_d])
        densities[point] = density
    
    tmp_densities = [(d, pt) for pt,d in densities.iteritems()]
    tmp_densities.sort(reverse=True)
    return tmp_densities, std_d

def _CalcMutualNearestNeighbors(hull_points, all_points):
    all_points_list = list(all_points)
    ds = distance.pdist(list(all_points))
    std_d = p.std(ds)
    
    square_ds = distance.squareform(ds)
    nearest_neighbors = {}
    
    for i, point in enumerate(all_points_list):
        if point not in hull_points:
            continue
        
        my_ds = [(d, j) for j, d in enumerate(square_ds[i])
                 if j != i]
        my_ds.sort()
        nearest_neighbors[point] = set([j for d,j in my_ds[:3]])
    
    no_mutual = set()
    for i, point in enumerate(all_points_list):
        if point not in hull_points:
            continue
        
        no_nbrs = True
        for neighbor_index in nearest_neighbors.get(point, []):
            neighbor = all_points_list[neighbor_index]
            neighbor_set = nearest_neighbors.get(neighbor, [])
            if i in neighbor_set:
                no_nbrs = False
        
        if no_nbrs:
            no_mutual.add(point)
                
    return no_mutual
    

def Test():
    t = convexhull.sample_unit_triangle(100, 2)
    noise = n.random.random_sample((4, 2))
    t = n.concatenate((t, noise))
    t = [tuple(pt) for pt in t]
    tri_coords = _LayoutPointSet(t)
    
    skin, leftover = OnionSkin(tri_coords)
    
    densities, rad = _CalcDensities(skin, t)
    no_nbrs = _CalcMutualNearestNeighbors(skin, t)
    
    least_densest_pts = [pt for _, pt in densities[-10:]]
    points_for_randomized_hull = p.array(list(skin - set(least_densest_pts)))
    no_nbrs = _LayoutPointSet(no_nbrs)
    least_densest_pts = _LayoutPointSet(least_densest_pts)
    
    skin = _LayoutPointSet(skin)
    leftover = _LayoutPointSet(leftover)
    
    
    
    pointies = convexhull.find_convex_hull(
        points_for_randomized_hull, 10000, 3)
    
    p.scatter(skin[0], skin[1], c='g')
    p.scatter(leftover[0], leftover[1], c='r')
    p.scatter(least_densest_pts[0], least_densest_pts[1], c='k',
              s=300, alpha=0.2)
    p.scatter(no_nbrs[0], no_nbrs[1], c='b',
              s=300, alpha=0.2)
    p.scatter(pointies[:,0], pointies[:,1], c='y')
    
    p.show()
    
if __name__ == '__main__':
    Test()
    