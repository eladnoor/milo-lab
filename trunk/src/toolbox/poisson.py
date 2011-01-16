import pylab
import random

def Simulate(n_buckets, n_balls, n_iterations=1000):
    result_list = []
    for _ in xrange(n_iterations):
        buckets = pylab.zeros((n_buckets))
        for _ in xrange(n_balls):
            k = int(random.random() * n_buckets)
            buckets[k] += 1
        n_single = len(pylab.find(buckets == 1))
        result_list.append(n_single)
        
    return pylab.mean(result_list)

def MultiSimulate(n_buckets, n_iterations=1000):
    averages = pylab.zeros((2*n_buckets))
    for n_balls in xrange(len(averages)):
        averages[n_balls] = Simulate(n_buckets, n_balls, n_iterations)
    pylab.plot(range(len(averages)), averages)
    pylab.show() 
    
if __name__ == "__main__":
    MultiSimulate(49, 1000)
        
        