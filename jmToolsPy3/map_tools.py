"""
map_tools: Convenience functions for processing X-Ray EDS maps
------------------------------------------------------------------------
    
 Ver       Date      Who  Comments
-------  ----------  ---  ------------------------------------------
0.0.959  2016-03-04  JRM  Add ensureDir and fix_gray_image_to_rgb utils
0.0.961  2016-03-10  JRM  Add medianFilterMaps
0.0.962  2016-03-14  JRM  Add rescaleFilteredMaps
"""
# -*- coding: utf-8 -*-

def readMapFilesTif(lNames, path='./tif/', bVerbose=False):
    """readMapFilesTif(lNames, path='./tif/', bVerbose=False)
    
    Read and scale a list of EDS maps with the ROI as the last element.
    Convert to float images, scaling to the such that the maximum counts
    in the most intense map has a gray value of 1.0
    
    """
    import numpy as np
    from skimage.external.tifffile import imread
    imgs = []
    maxCts = []
    for name in lNames:
        fPath = path + name + '.tif'
        if bVerbose:
            print(fPath)
        img = imread(fPath)
        img = img.astype(float)
        maxCts.append(np.max(img))
        imgs.append(img)
    l = len(imgs)
    if bVerbose:
        print(l)
        print(maxCts)
    maxG = np.max(maxCts[0:l-1])
    roiG = maxCts[l-1]
    if bVerbose:
        print(maxG, roiG)
    for i in range(l-1):
        imgs[i] /= maxG
    imgs[l-1] /= roiG

    return imgs, maxCts


def denoiseMaps(imgs, wt = 0.05, ep = 0.02, cycles = 200):
    """denoiseMaps(imgs, wt = 0.05, ep = 0.02, cycles = 200)
    
    A wrapper for denoise_tv_chambolle
    
    returns
    
    A list of denoised images
    """
    from skimage.restoration import denoise_tv_chambolle
    denoised = []
    for img in imgs:
        dn = denoise_tv_chambolle(img, weight=wt, eps=ep, n_iter_max=cycles, multichannel=False)
        denoised.append(dn)

    return denoised

def medianFilterMaps(imgs, px=3):
    """medianFilterMaps(imgs, px=3)
    
    A wrapper for ndimage median_filter
    
    Parameters
    ----------
    imgs: a list of ndarrays (2D)
        Input images, assumed to be grayscale
    px: number (default 3)
        The size of the median filter
    
    A list of median filtered images
    """
    from scipy import ndimage as ndi
    filtered = []
    for img in imgs:
        mf = ndi.filters.median_filter(img, size=(px,px))
        filtered.append(mf)

    return filtered

def rescaleFilteredMaps(imgs, skipBlank=False, bVerbose=False):
    """rescaleMaps(imgs, skipBlank=False, bVerbose=False)

    Rescale the map intensities of noise-filtered maps. Treat
    the ROI by itself. Skip a blank if needed."""
    import numpy as np
    maxCts = []
    for img in imgs:
        maxCts.append(np.max(img))
    l = len(imgs)
    if bVerbose:
        print(l)
        print(maxCts)
    if skipBlank:
        maxG = np.max(maxCts[0:l-2])
        maxB = maxCts[l-2]
        roiG = maxCts[l-1]
        for i in range(l-2):
            imgs[i] /= maxG
        imgs[l-2] /= maxB
        imgs[l-1] /= roiG
    else:
        maxG = np.max(maxCts[0:l-1])
        roiG = maxCts[l-1]
        for i in range(l-1):
            imgs[i] /= maxG
        imgs[l-1] /= roiG
    return imgs



