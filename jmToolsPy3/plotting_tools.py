"""
plotting_tools: Convenience functions for plotting images and data
------------------------------------------------------------------------
    
 Ver        Date       Who    Comments
------    ----------   ---    ------------------------------------------
0.0.90    2015-08-10   JRM    Initial prototype: complexPolarPlot
0.0.92    2015-12-31   JRM    move to PEP8 and add functions
0.0.935   2015-12-31   JRM    Added plotImageWithHistogram and improved
                              the comments.
"""
# -*- coding: utf-8 -*-
def complexPolarPlot(a, color='#CC0000'):
    """complexPolarPlot(a)
    Plot a list or vector of complex points as a polar plot.
    Defaults to red points.

    Parameters
    ----------
    a     : a vector or list of points
    color : a color string (Default red, '#CC0000')
    
    Returns
    -------
    fig   : a matplotlib Figure
    """
    from matplotlib import pyplot as plt
    import numpy as np
    for x in a:
        plt.polar(np.angle(x),np.abs(x), marker='.',color=color)
    fig = plt.gcf()
    return fig

def plotImageMontage(lData, lTitles, nrows, ncols, title=None, baseSize=5, cmap='gray' ,fSize=12, bare=False, cb=False, **kwargs):
    """plotImageMontage(lData, lTitles, nrows, ncols, title=None, cmap='gray',
                         baseSize=5, fSize=12, bare=False, cb=False, **kwargs)
    Create a figure from an ensemble of images

    Parameters
    ----------
    lData    : A list of images (numpy arrays)
               These are the images to be plotted.

    lTitles  : A list of (strings)
               These are the titles corresponding to the images

    nrows    : An integer
               The number of rows for the figure

    ncols    : An integer
               The number of columns for the figure. 
               Note: nrows x ncols should equal the number of images in lData

    title    : A string
               The title for the figure, if desired. Defaults to None

    cmap     : a matplotlib color map.
               Default is 'gray'. One can use a string or a
               colormap like plt.cm.gray or plt.cm.viridis.

    baseSize : A number 
               The base size in inches for an image. The function will
               calculate the figure size. Default is 5.

    fSize    : An integer
               The font size for title of each image. Default is 12.
               The figure title will be 2 points larger.

    bare     : a boolean flag (default False).
               If True, the axis scales are turned off.

    cb       : a boolean flag (default True).
               If True a color bar is plotted on the right side of the each image.

    **kwargs : Any other kword arguments. Useful arguments are vmin, vmax,
               interpolation, and colormap. It accepts string formates like
               'gray' out of the box. If you import matplotlib.pyplot as plt,
               you can input all theformats like plt.cm.spectal, plt.cm.gray,
               plt.cm.viridis. The defalt is the string 'gray'.

    Returns
    --------
    fig      : A matplotlib figure


    Examples
    --------

    from skimage import color, data
    from jmToolsPy import plotImageMontage

    im1 = color.rgb2gray(data.astronaut())
    im2 = data.coins()

    myFig = plotImageMontage([im1, im2], ["One","Two"], 1, 2, title="Figure", bare=True, cb=True)
    """
    from matplotlib import pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    xSize = ncols*baseSize
    ySize = nrows*baseSize
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(xSize, ySize))
    if title != None:
        fig.suptitle(title, fontsize=fSize+2)
    for ax, im, ti in zip(axes.flat, lData, lTitles):
        ax.set_title(ti, fontsize=fSize)
        cim = ax.imshow(im, cmap=cmap, **kwargs)
        if cb == True:
            # Create divider for existing axes instance
            div = make_axes_locatable(ax)
            # Append axes to the right of ax, with 10% width of ax
            cax = div.append_axes("right", size="10%", pad=0.05)
            # Create colorbar in the appended axes
            cbar = plt.colorbar(cim, cax=cax)
        # Remove xicks from ax
        if bare == True:
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
    plt.show();
    return(fig)


def iter_channels(color_image):
    """Yield color channels of an image."""
    # Roll array-axis so that we iterate over the color channels of an image.
    import numpy as np
    for channel in np.rollaxis(color_image, -1):
        yield channel

def plotImageWithHistogram(im, size, alpha=0.3):
    """plotImageWithHistogram(im, size, alpha=0.3)

    Plot an image alongside its histogram

    Parameters
    ----------

    im      : a numpy array
              The input image

    size    : a number
              The size (in in) for the single figure. The plot will be
              that high and twice that wide.

    alpha   : a float
              The transparency. A value of 0.3 is great for an RGB image,
              a value of 0.8 is a bit better for a grayscale image.

    Returns
    -------
    a tuple   (ax__image, ax_hist) 
              The matplotlib axes for the figure.

    Examples
    --------

    from jmToolsPy3 import plotImageWithHistogram
    import numpy as np
    from skimage import data

    img1 = data.camera()
    axImg1, axHis1 = plotImageWithHistogram(img1, 5, alpha=0.8)

    img2 = data.lena()
    axImg2, axHis2 = plotImageWithHistogram(img2, 5, alpha=0.3)
    """
    from skimage import exposure
    from matplotlib import pyplot as plt
    fig, (ax_image, ax_hist) = plt.subplots(ncols=2, figsize=(2*size, size))
    ax_image.imshow(im, cmap=plt.cm.gray)
    if im.ndim == 2:
        hist, bin_centers = exposure.histogram(im)
        ax_hist.fill_between(bin_centers, hist, alpha=alpha, color='gray')        
    elif im.ndim == 3:
        for channel, channel_color in zip(iter_channels(im), 'rgb'):
            hist, bin_centers = exposure.histogram(channel)
            ax_hist.fill_between(bin_centers, hist, alpha=alpha, color=channel_color)
    ax_hist.set_ylabel('# pixels')
    ax_hist.set_xlabel('intensity')
    ax_hist.set_yticklabels("")
    ax_image.set_axis_off()
    # match_axes_height(ax_image, ax_hist)
    dst = ax_hist.get_position()
    src = ax_image.get_position()
    ax_hist.set_position([dst.xmin, src.ymin, dst.width, src.height])
    return ax_image, ax_hist





