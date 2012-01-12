#import cv2.cv as cv
import scipy.misc as misc
import scipy.ndimage as ndimage
import numpy as np
import matplotlib.pyplot as plt
import pymorph

def GetCircle(center, radius):
    x_points = []
    y_points = []
    for x in xrange(radius+1):
        l = np.sqrt(2*x*radius - x**2)
        r_y = range(int(np.floor(center[1]-l)), int(np.floor(center[1]+l)))
        y_points += r_y
        x_points += [int(np.floor(center[0] - radius + x))]*len(r_y)
        if x < radius:
            y_points += r_y
            x_points += [int(np.floor(center[0] + radius - x))]*len(r_y)
    return (x_points, y_points)

def NumpyDetectEdge(fname, color_filter=(1, 1, 1), 
                    radius=510, center=(582, 791), show=False):
    plate_circle = GetCircle(center, radius)
    
    im = misc.imread(fname)
    im_gray = np.dot(im, color_filter) / sum(color_filter)

    imf = ndimage.gaussian_filter(im_gray, sigma=2)
    imf = 256 - imf # invert color

    imc = np.zeros(imf.shape, dtype=np.int32)
    imc[plate_circle] = imf[plate_circle]
    imc = imc[center[0]-radius:center[0]+radius, center[1]-radius:center[1]+radius]

    T = 100
    imc[np.where(imc < T)] = 0
    
    rmax = pymorph.regmax(imc)
    seeds, n_colonies = ndimage.label(rmax) # gives a unique integer number to each region and fills it with that number
    
    count = 0
    for i in xrange(n_colonies):
        seed = np.where(seeds==i)
        l = len(seed[0])
        # filter seeds which are touching the perimeter (or 5 pixel away from it) 
        perimeter = False
        for j in xrange(l):
            if np.sqrt((radius-seed[0][j])**2 + (radius-seed[1][j])**2) > radius-5:
                perimeter = True
        if perimeter:
            continue

        seedx = [(x + center[0] - radius) for x in seed[0]]
        seedy = [(y + center[1] - radius) for y in seed[1]]
        meanx = np.mean(seedx)
        meany = np.mean(seedy)
        
        colony_circle = GetCircle((meanx, meany), 3)
        colony_color = np.min(imf[colony_circle])
        if colony_color >= 110:
            count += 1
        if colony_color < 100:
            im[meanx, meany, :] = (255, 0, 0) # red
        elif 100 < colony_color < 110:
            im[meanx, meany, :] = (255, 255, 0) # yellow
        elif 110 < colony_color < 120:
            im[meanx, meany, :] = (0, 255, 0) # green 
        elif 120 < colony_color < 130:
            im[meanx, meany, :] = (0, 255, 255) # cyan
        elif 130 < colony_color < 140:
            im[meanx, meany, :] = (0, 0, 255) # blue
        elif 140 < colony_color < 150:
            im[meanx, meany, :] = (255, 0, 255) # magenta
        elif 150 < colony_color:
            im[meanx, meany, :] = (255, 255, 255) # white
    
    if show:
        plt.imshow(im)
        plt.show()
    return count

if __name__ == "__main__":
    fnames = []
    #fnames += ['/home/eladn/Dropbox/Experiments/ace-plates/2011-12-29 ace- 10^5.jpg']
    #fnames += ['/home/eladn/Dropbox/Experiments/ace-plates/2011-12-29 ace+ 10^5.jpg']
    fnames += ['/home/eladn/Dropbox/Experiments/ace-plates/2012-01-03 ace- 10^3.jpg']
    fnames += ['/home/eladn/Dropbox/Experiments/ace-plates/2012-01-03 ace+ 10^3.jpg']
    
    for fname in fnames:
        print fname, NumpyDetectEdge(fname, color_filter=(0, 0, 1)) # take only the blue color
    #img = cv.LoadImage(fname)
    #DetectHough(img)
    #DetectEdge(img, 7)
