"""
plotting_tools: Convenience functions for plotting images and data
------------------------------------------------------------------------
    
 Ver       Date      Who  Comments
-------  ----------  ---  ------------------------------------------
0.0.90   2015-08-10  JRM  Initial prototype: complexPolarPlot
0.0.92   2015-12-31  JRM  move to PEP8 and add functions
0.0.935  2015-12-31  JRM  Added plotImageWithHistogram and improved
                          the comments.
0.0.940  2016-01-06  JRM  Numpy doc string format
0.0.945  2016-01-06  JRM  More functions
0.0.950  2016-02-23  JRM  Added showImages
0.0.955  2016-03-01  JRM  Added watershedBlobAnalysis
0.0.956  2016-03-01  JRM  Some tweaks to output
0.0.957  2016-03-01  JRM  Some more tweaks to output
0.0.959  2016-03-05  JRM  Add ensureDir and fix_gray_image_to_rgb utils
0.0.960  2016-03-07  JRM  Added clip_img_hi
0.0.965  2016-03-24  JRM  Added flag for colorbar to plotImage
"""
# -*- coding: utf-8 -*-

def clip_img_hi(img, cut=10):
    """
    clip_img_hi(img, cut=10)

    Clip a number of spurious pixels from the upper end of the gray
    level histogram of an image.

    Parameters
    ----------
    img: ndarray (2D)
        Input image, assumed to be grayscale
    cut: number (default 10)
        The number of bright pixels to clip

    Returns

    clipped: the clipped image
    """
    import numpy as np
    his, bins = np.histogram(img, bins=np.max(img)+1)
    n = np.max(img)
    sumCut = 0
    while sumCut <= cut:
        n -= 1
        sumCut += his[n]
    clipped = np.clip(img, 0, n)
    return clipped

def ensureDir(d):
    """ensureDir(d)
    Check if the directory, d, exists, and if not create it."""
    import os
    if not os.path.exists(d):
        os.makedirs(d)

def fix_gray_image_to_rgb(img):
    """fix_gray_image_to_rgb(img)
    Convert an image to float, scale from 0.0 to 255.0,
    and convert to RGB"""

    import numpy as np
    from skimage import img_as_ubyte
    from skimage.color import gray2rgb

    img = img.astype(float)
    maxV = np.max(img)
    minV = np.min(img)
    fact = 255.0/(maxV-minV)
    img -= minV
    img *= fact
    img = img.astype(int)
    img = img_as_ubyte(img)
    img = gray2rgb(img)
    return img

def watershedBlobAnalysis(img, thr, bright=True, showPlot=True, sig=3, pkSz=3, minPx=10, sf=1.79):
    """
    Binarize an input image using the supplied threshold, perform a
    watershed analysis on a smoothed Euclidean Distance Map, optionally
    display boundaries on the image, and return a list of lists of the
    properties of the blobs.

    Parameters
    ----------
    img: ndarray (2D)
        Input image, assumed to be grayscale
    thr: number
        The value to use for the gray level threshold
    bright: boolean (True)
        A flag indicating if the background is bright, and the blobs
        have gray levels less than the threshold value or the reverse.
    showPlot: boolean (True)
        A flag indicating whether to display a plot of the boundaries
        displayed in red on the grayscale input image
    sig: number (3)
        The sigma parameter for a gaussian smooth of the Euclidean
        Distance Map
    pkSz: number (3)
        The peak size of the footprint for the call to find the peak 
        local maxima.
    minPx: integrer (10)
        the minimum number of pixels to conside a "blob"
    sf: float (1.79)
        The scale factor (units/px) for the image. The default is for a
        test image of AgX grains where the scale factor is 1.79 nm/px.

    Returns

    out: a list of lists
        [equivalent circular diameter, centroid, aspect ratio, solidity,
        circircularity, squareness]
        Lists of feature vectors useful for particle size and shape
        analysis. Only the ECD is in dimensions implied by the scale
        factor.
        
    """
    from math import sqrt
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import ndimage
    from skimage import img_as_ubyte
    from skimage.morphology import disk, watershed, remove_small_objects
    from skimage.filters.rank import median
    from skimage.feature import peak_local_max
    from skimage.segmentation import find_boundaries, mark_boundaries
    from skimage.measure import regionprops
    from skimage.color import gray2rgb
    
    if bright:
        bin_img = img < thr
    else:
        bin_img = img > thr
    bin_img = remove_small_objects(bin_img, min_size=minPx, connectivity=1,
                                   in_place=False)
    dist = ndimage.distance_transform_edt(bin_img)
    smooth = ndimage.gaussian_filter(dist, sig)
    # peak_local_max
    # possible below: #,labels=c.binarybackground)
    
    local_maxi = peak_local_max(smooth, indices=False, footprint=np.ones((pkSz, pkSz)))
    markers = ndimage.label(local_maxi)[0]  
    labels = watershed(-smooth, markers, mask=bin_img)
    props = regionprops(labels)
    ecd = []
    cent = []
    ar = []
    solid =[]
    circ = []
    square = []
    for prop in props:
        cent.append(prop.centroid)
        ecd.append(round(sf*prop.equivalent_diameter, 3))
        if prop.minor_axis_length == 0:
            ar.append(float('nan'))
        else:
            ar.append(round(prop.major_axis_length/prop.minor_axis_length, 3))
        solid.append(prop.solidity)
        fArea = float(prop.area)
        perim = prop.perimeter
        cir = 4.0 * np.pi * (fArea / (perim)**2)
        circ.append(cir)
        square.append(0.25*perim/sqrt(fArea))
    #find outline of objects for plotting
    
    if showPlot:
        boundaries = find_boundaries(labels)
        max_g = np.max(img)
        if max_g > 255:
            img.astype(np.float32)
            img = img / float(max_g)
        img_rgb = gray2rgb(img) # transform 16bit to 8
        overlay = np.flipud(mark_boundaries(img_rgb, boundaries, color=(1, 0, 0)))
        # plt.imshow(overlay);
    
        fig = plt.figure(figsize=(7,7))
        ax = fig.add_subplot(1, 1, 1)
        ax.imshow(overlay, cmap='spectral');
        ax.xaxis.set_visible(False);
        ax.yaxis.set_visible(False)
        fig.set_tight_layout(True);
    
    return ([ecd, cent, ar, solid, circ, square])

def showImages(images,titles=None, cmap='gray', bare=False):
    """Display a list of images

    Parameters
    ----------
    images : [images] (a list of numpy arrays)
        These are the images to be plotted.

    titles: A list of strings (None)
        These are the titles corresponding to the images. Default is None

    cmap: a matplotlib color map ('gray').
        Default is 'gray'. One can use a string or a
        colormap like plt.cm.gray or plt.cm.viridis.

    bare: Boolean (False)
        A flag to supress axis numbering
    """
    from matplotlib import pyplot as plt
    import numpy as np

    n_ims = len(images)
    if titles is None: titles = ['(%d)' % i for i in range(1,n_ims + 1)]
    fig = plt.figure()
    n = 1
    for image,title in zip(images,titles):
        a = fig.add_subplot(1,n_ims,n) # Make subplot
        # if image.ndim == 2: # Is image grayscale?
        #     plt.gray() # Only place here where you can't replace 'gray' with 'grey'
        plt.imshow(image, cmap=cmap)
        a.set_title(title)
        if bare:
            a.xaxis.set_visible(False);
            a.yaxis.set_visible(False);
        n += 1
    fig.set_size_inches(np.array(fig.get_size_inches()) * n_ims)
    plt.show()


def complexPolarPlot(a, color='#CC0000'):
    r"""Plot a list or vector of complex points as a polar plot.
    Defaults to red points.

    Parameters
    ----------
    a : a vector or list of complex points

    color : string '#CC0000', optionall
        a color string (Default red, '#CC0000')
    
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
    r"""Create a figure from an ensemble of images

    Parameters
    ----------
    lData : A list of images (numpy arrays)
        These are the images to be plotted.

    lTitles : A list of strings
        These are the titles corresponding to the images

    nrows : int
        The number of rows for the figure

    ncols : int
        The number of columns for the figure. 
        Note: nrows x ncols should equal the number of images in lData

    title : string
        The title for the figure, if desired. Defaults to None

    cmap : a matplotlib color map.
        Default is 'gray'. One can use a string or a
        colormap like plt.cm.gray or plt.cm.viridis.

    baseSize : number 5, optional
        The base size in inches for an image. The function will
        calculate the figure size.

    fSize : int 12, optional
        The font size for title of each image. The figure title will
        be 2 points larger.

    bare : boolean False, optional
        If True, the axis scales are turned off.

    cb : boolean True, optional
        If True a color bar is plotted on the right side of the each image.

    **kwargs : Any other kword arguments.
        Useful arguments are vmin, vmax,
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
    from jmToolsPy3 import plotImageMontage

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
    fig.tight_layout();
    plt.show();
    return(fig)


def iter_channels(color_image):
    r"""Yield color channels of an image."""
    # Roll array-axis so that we iterate over the color channels of an image.
    import numpy as np
    for channel in np.rollaxis(color_image, -1):
        yield channel

def plotImageWithHistogram(im, size, alpha=0.3):
    r"""Plot an image alongside its histogram

    Parameters
    ----------

    im : a numpy array
        The input image

    size : number
        The size (in in) for the single figure. The plot will be
        that high and twice that wide.

    alpha : float
        The transparency. A value of 0.3 is great for an RGB image,
        a value of 0.8 is a bit better for a grayscale image.

    Returns
    -------
    (ax__image, ax_hist) 
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


def plotImage(im, cmap='gray' , figsize=(7,8), cb=False, bare=True):
    r"""Plot an image in a tight layout

    Parameters
    ----------

    im : a numpy array
        The input image

    cmap : a string or plt.cm.colormap
        The colormap

    figsize : tuple (width, height), (7,8). default
        The size (in in) for the single figure.

    cb : boolean False, default
        Add a colorbar

    bare: boolean True, default
        Flag for no axis labels

    Returns
    -------
    None.

    Example
    --------

    from jmToolsPy3 import plotImage
    from skimage import data

    img1 = data.camera()
    plotImage(img1, cmap='viridis', figsize=(5,5), bare=False)
    """
    from matplotlib import pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(1, 1, 1)
    cim = ax.imshow(im, cmap=cmap);
    if bare:
        ax.xaxis.set_visible(False);
        ax.yaxis.set_visible(False);
    if cb == True:
        # Create divider for existing axes instance
        div = make_axes_locatable(ax)
        # Append axes to the right of ax, with 10% width of ax
        cax = div.append_axes("right", size="10%", pad=0.05)
        # Create colorbar in the appended axes
        cbar = plt.colorbar(cim, cax=cax)
    else:
        fig.set_tight_layout(True)


