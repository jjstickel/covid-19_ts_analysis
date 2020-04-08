"""
Provides convenience 1D numerical derivative function
"""
# Jonathan Stickel 2010, 2011, 2012, 2013

# TODO:
# - in 'deriv1_fd', calculate second-order accurate forward and backward
# difference for endpoints (currently only first-order accurate)

import numpy as np

def is_even(val):
    return np.mod(val, 2) == 0

def deriv1_fd(y,x=None,h=None,xmid=False,central=False,e=0.0):
    """
    Computes the approximate derivative of the y vs x data by finite
    difference.

    Parameters
    ___________________
    y : 1D array
        y data values of array length N
    x : 1D array
        monotonically increasing x data values of array length N; need
        not be provided if equally spaced and h is provided.
    h : float
        Uniform grid spacing of x.  If x is not equally spaced, x
        should be provided and h should be None.
    xmid : boolean
        If true, then midpoint values of x are returned along with the
        the derivative values
    central : boolean
        Performs a central difference, correct for unequally spaced x
        values. The endpoints are calculate by forward and backward
        difference so that yp is the same length as x.
    e : float
        Backward/forward adjustment to the central difference
        calculation: -1<e<0 shift to backward differnce, 0<e<1 shift
        to forward difference.
        
    Returns
    _______
    yp : 1D array
        The derivative values; array length N-1 (length N if
        central=True)
    xm : 1D array (optional)
        The midpoint values of x; array length N-1

    see also:  numpy.gradient
    """

    if xmid and central:
        raise Exception('Both xmid and central cannot both be True')

    # calculate direct differences
    if h is not None:
        # this calculation could potentially be avoided if doing central difference
        yp = np.diff(y)/h
    else:
        yp = np.diff(y)/np.diff(x.astype(float))
    
    if central:
        yd21 = np.diff(y[::2]) # every-other difference of y[i+1] - y[i-1]
                               # starting at i = 1
        yd22 = np.diff(y[1::2]) # every-other difference of y[i+1] - y[i-1]
                                # starting at i = 2
        if is_even(y.size):
            # interleave the two 1d arrays
            yd2 = np.reshape(np.column_stack((yd21,yd22)),-1)
        else:
            yd22 = np.hstack((yd22,np.nan))
            yd2 = np.reshape(np.column_stack((yd21,yd22)),-1)
            yd2 = yd2[:-1]
        if h is not None:
            ypc = yd2/(2*h)
        else:
            xd21 = np.diff(x[::2])
            xd22 = np.diff(x[1::2])
            if is_even(x.size):
                xd2 = np.reshape(np.column_stack((xd21,xd22)),-1)
            else:
                xd22 = np.hstack((xd22,np.nan))
                xd2 = np.reshape(np.column_stack((xd21,xd22)),-1)
                xd2 = xd2[:-1]
            ypc = yd2/xd2.astype(float) # central difference
        # paramater e can be used to shift towards backward/forward difference
        ypc = (1-e)*yp[:-1] + (1+e)*yp[1:] - ypc
        return np.hstack( ( yp[0], ypc, yp[-1] ) )
    else:
        if xmid:
            if x is None:
                raise Exception('x must be provided if xmid=True')
            if h is not None:
                xm = x[:-1] + h/2.0
            else:
                xm = ( x[:-1] + x[1:] )/2.0
            return yp,xm
        else:
            return yp

