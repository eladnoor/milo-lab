import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import numpy as np

def contourf():
    """
        Draws a 2D Gaussian centered at the middle of the range.
    """
    x_values = np.r_[-1:1:100j]
    y_values = np.r_[-1:1:100j]
    sigma = 2
    data = np.exp(-np.add.outer(x_values**2, y_values**2) / 2*sigma**2)
    cdict = {'red': ((0.0, 0.0, 0.0),
                    (0.5, 1.0, 1.0),
                    (1.0, 0.2, 0.2)),
            'green': ((0.0, 0.0, 0.0),
                      (0.5, 0.0, 1.0),
                      (1.0, 0.8, 0.8)),
            'blue': ((0.0, 0.0, 0.0),
                     (0.5, 0.0, 0.0),
                     (1.0, 0.2, 0.2))}
    cmap = colors.LinearSegmentedColormap('red_green', cdict, 256)
    
    
    fig = plt.figure()
    plt.contourf(x_values, y_values, data,
                 levels=np.r_[0:1:11j], cmap=cmap)
    plt.colorbar()
    fig.axes[0].set_xlim(-2, 2)
    fig.axes[0].set_ylim(-2, 2)
    #plt.text(0, 0, 'origin', ha='left', va='center')
    plt.show()
    
if __name__ == "__main__":
    contourf()