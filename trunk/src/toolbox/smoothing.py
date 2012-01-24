
from scipy import linalg
import numpy as np



class WeightedAverageSmoother(object):
    
    def __init__(self, xs, ys, sigma=1):
        self.xs = np.array(xs)
        self.ys = np.array(ys)
        self.sigma = sigma
        
    def _GetY(self, x):
        diffs = np.abs(self.xs - x)
        dists = np.exp(diffs ** 2 / (2*self.sigma ** 2))
        sims = 1.0 / dists
        weights = sims / np.sum(sims)
        return np.sum(np.dot(self.ys, weights))
    
    def __call__(self, xs):
        return np.array(map(self._GetY, xs))