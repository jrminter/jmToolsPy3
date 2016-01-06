def peak_detect(y_axis, x_axis = None, lookahead = 500, delta = 0):
    r"""
    Detect peaks in a list

    Converted from/based on a MATLAB script at http://billauer.co.il/peakdet.html
    
    From https://gist.githubusercontent.com/schlady/1576079/raw/03449c7af8bec6a4744b4cd69a671359ab6a9afe/peakdetect.py
    
    JM update to Python3 exception handling
    
    Algorithm for detecting local maximas and minmias in a signal.
    Discovers peaks by searching for values which are surrounded by lower
    or larger values for maximas and minimas respectively
    
    Parameters
    ----------
    y_axis : array_like
        A list containg the signal over which to find peaks

    x_axis : array_like {None} (optional)
        A x-axis whose values correspond to the 'y_axis' list and is used
        in the return to specify the postion of the peaks. If omitted the
        index of the y_axis is used.

    lookahead : int {500}, optional
        distance to look ahead from a peak candidate to determine if it
        is the actual peak. '(sample / period) / f' where '4 >= f >= 1.25'
        might be a good value.

    delta : int must be positive {0}, optional
        This specifies a minimum difference between a peak and the
        following points, before a peak may be considered a peak. Useful
        to hinder the algorithm from picking up false peaks towards to end of
        the signal. To work well delta should be set to 'delta >= RMSnoise * 5'.
        Delta function causes a 20% decrease in speed, when omitted
        Correctly used it can double the speed of the algorithm
    
    Returns
    -------
    [maxtab, mintab] : Two lists containing the positive and negative
                       peaks respectively. Each cell of the lists
                       contains a tupple of: (position, peak_value) 
                       to get the average peak value do
                       'np.mean(maxtab, 0)[1]' on the results
    """
    import numpy as np
    maxtab = []
    mintab = []
    dump = []   #Used to pop the first hit which always if false
       
    length = len(y_axis)
    if x_axis is None:
        x_axis = range(length)
    
    #perform some checks
    try:
        length == len(x_axis)
    except ValueError:
        print("Input vectors y_axis and x_axis must have same length")

    try:
        lookahead >= 1
    except ValueError:
        print("Lookahead must be above '1' in value")

    try:
        np.isscalar(delta) and (delta >= 0)
    except ValueError:
        print("delta must be a positive number")
    
    #needs to be a numpy array
    y_axis = np.asarray(y_axis)
    
    #maxima and minima candidates are temporarily stored in
    #mx and mn respectively
    mn, mx = np.Inf, -np.Inf
    
    #Only detect peak if there is 'lookahead' amount of points after it
    for index, (x, y) in enumerate(zip(x_axis[:-lookahead], y_axis[:-lookahead])):
        if y > mx:
            mx = y
            mxpos = x
        if y < mn:
            mn = y
            mnpos = x
        
        ####look for max####
        if y < mx-delta and mx != np.Inf:
            #Maxima peak candidate found
            #look ahead in signal to ensure that this is a peak and not jitter
            if y_axis[index:index+lookahead].max() < mx:
                maxtab.append((mxpos, mx))
                dump.append(True)
                #set algorithm to only find minima now
                mx = np.Inf
                mn = np.Inf
        
        ####look for min####
        if y > mn+delta and mn != -np.Inf:
            #Minima peak candidate found 
            #look ahead in signal to ensure that this is a peak and not jitter
            if y_axis[index:index+lookahead].min() > mn:
                mintab.append((mnpos, mn))
                dump.append(False)
                #set algorithm to only find maxima now
                mn = -np.Inf
                mx = -np.Inf
    
    
    #Remove the false hit on the first value of the y_axis
    try:
        if dump[0]:
            maxtab.pop(0)
            #print "pop max"
        else:
            mintab.pop(0)
            #print "pop min"
        del dump
    except IndexError:
        #no peaks were found, should the function return empty lists?
        pass
    
    return maxtab, mintab
