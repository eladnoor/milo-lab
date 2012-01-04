#import cv2.cv as cv
import scipy.misc as misc
import scipy.ndimage as ndimage
import numpy as np
import matplotlib.pyplot as plt
import mahotas
import pymorph

def NumpyLoadImage(fname):
    #im = Image.open(fname)
    #_, _, height, width = im.getbbox()
    #bitmap = np.zeros((width, height, 3))
    #for i, (r, g, b) in enumerate(im.getdata()):
    #    bitmap[i/height, i%height, :] = (1.0-r/256.0, 1.0-g/256.0, 1.0-b/256.0)
    
    bitmap = misc.imread(fname) / 256.0
    greyscale = np.mean(bitmap, 2)
    threshold = 0.7
    greyscale[np.where(greyscale > threshold)] = 1.0
    greyscale[np.where(greyscale <= threshold)] = 0.0
    
    plt.imshow(greyscale, cmap=plt.cm.gray)
    plt.show()

def DetectHough(img):
    grey = cv.CreateImage((img.width, img.height), cv.IPL_DEPTH_8U, 1)
    cv.CvtColor(img, grey, cv.CV_RGB2GRAY)
    
    thresholded = cv.CreateImage((img.width, img.height), cv.IPL_DEPTH_8U, 1)
    cv.InRangeS(grey, 150, 256, thresholded)
    storage = cv.CreateMat(1, 1000, cv.CV_32FC3)      
    cv.HoughCircles(thresholded, storage, cv.CV_HOUGH_GRADIENT, 1, 
                    thresholded.height/4, 100, 40, 20, 200)
    print storage[0, 0]
    #for i in xrange(storage.width):
    #    print storage[i, 0], storage[i, 1], storage[i, 2]

def DetectEdge(img, N=7):
    img_b = cv.CreateImage((img.width+N-1, img.height+N-1), 8, 3)
    img_g = cv.CreateImage((img.width+N-1, img.height+N-1), 8, 1)
    out = cv.CreateImage((img.width+N-1, img.height+N-1), 8, 1)
    
    # Add convolution borders
    offset = ((N-1)/2, (N-1)/2)
    cv.CopyMakeBorder(img, img_b, offset, 0, 0)
    cv.CvtColor(img_b, img_g, cv.CV_RGB2GRAY)
    
    # Edge Detection Variables
    aperture_size = N
    lowThresh = 60.0 * N**2
    highThresh = 80.0 * N**2
    
    cv.Canny(img_g, out, lowThresh, highThresh, aperture_size=aperture_size)
    cv.ShowImage("edge", out)

def NumpyDetectEdge(fname, color_filter=(1, 1, 1), radius=510, center=(582, 791)):
    im = misc.imread(fname)
    im_gray = np.dot(im, color_filter) / sum(color_filter)
    imf = ndimage.gaussian_filter(im_gray, sigma=5)
    imf = 256 - imf # invert color
    
    # cut out all except the circle containing the plate
    for x in xrange(im.shape[0]):
        l_sqr = radius**2 - (x-center[0])**2
        if l_sqr < 0:
            im[x, :, :] = 0
        else:
            l = np.sqrt(l_sqr)
            imf[x, 0:(center[1]-l)] = 0
            imf[x, (center[1]+l):] = 0
    
    T = 100
    imf[np.where(imf < T)] = 0
    
    rmax = pymorph.regmax(imf)
    seeds, n_colonies = ndimage.label(rmax) # gives a unique integer number to each region and fills it with that number
    for i in xrange(n_colonies):
        av = np.mean(imf[np.where(seeds==i)])
        if av < T:
            continue
        x = np.mean(np.where(seeds==i)[0])
        y = np.mean(np.where(seeds==i)[1])
        im[x, y, :] = (255, 0, 0)
    #im_gray = pymorph.overlay(im_gray, rmax)
    plt.imshow(im)
    #plt.imshow(im_gray, plt.cm.gray)
    plt.show()

if __name__ == "__main__":
    fnames = []
    fnames += ['/home/eladn/Dropbox/Experiments/ace-plates/2011-12-29 ace- 10^5.jpg']
    fnames += ['/home/eladn/Dropbox/Experiments/ace-plates/2011-12-29 ace+ 10^5.jpg']
    fnames += ['/home/eladn/Dropbox/Experiments/ace-plates/2012-01-03 ace- 10^5.jpg']
    fnames += ['/home/eladn/Dropbox/Experiments/ace-plates/2012-01-03 ace+ 10^5.jpg']
    NumpyDetectEdge(fnames[3], color_filter=(0, 0, 1)) # take only the blue color
    #img = cv.LoadImage(fname)
    #DetectHough(img)
    #DetectEdge(img, 7)
