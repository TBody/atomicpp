from __future__ import division, print_function, absolute_import

__all__ = [
    'UnivariateSpline',
    'InterpolatedUnivariateSpline',
    'LSQUnivariateSpline',
    'BivariateSpline',
    'LSQBivariateSpline',
    'SmoothBivariateSpline',
    'LSQSphereBivariateSpline',
    'SmoothSphereBivariateSpline',
    'RectBivariateSpline',
    'RectSphereBivariateSpline']

_surfit_messages = {1:"""
The required storage space exceeds the available storage space: nxest
or nyest too small, or s too small.
The weighted least-squares spline corresponds to the current set of
knots.""",
                    2:"""
A theoretically impossible result was found during the iteration
process for finding a smoothing spline with fp = s: s too small or
badly chosen eps.
Weighted sum of squared residuals does not satisfy abs(fp-s)/s < tol.""",
                    3:"""
the maximal number of iterations maxit (set to 20 by the program)
allowed for finding a smoothing spline with fp=s has been reached:
s too small.
Weighted sum of squared residuals does not satisfy abs(fp-s)/s < tol.""",
                    4:"""
No more knots can be added because the number of b-spline coefficients
(nx-kx-1)*(ny-ky-1) already exceeds the number of data points m:
either s or m too small.
The weighted least-squares spline corresponds to the current set of
knots.""",
                    5:"""
No more knots can be added because the additional knot would (quasi)
coincide with an old one: s too small or too large a weight to an
inaccurate data point.
The weighted least-squares spline corresponds to the current set of
knots.""",
                    10:"""
Error on entry, no approximation returned. The following conditions
must hold:
xb<=x[i]<=xe, yb<=y[i]<=ye, w[i]>0, i=0..m-1
If iopt==-1, then
  xb<tx[kx+1]<tx[kx+2]<...<tx[nx-kx-2]<xe
  yb<ty[ky+1]<ty[ky+2]<...<ty[ny-ky-2]<ye""",
                    -3:"""
The coefficients of the spline returned have been computed as the
minimal norm least-squares solution of a (numerically) rank deficient
system (deficiency=%i). If deficiency is large, the results may be
inaccurate. Deficiency may strongly depend on the value of eps."""
                    }

import warnings

from numpy import zeros, concatenate, alltrue, ravel, all, diff, array, ones
import numpy as np

# from . import fitpack
from . import dfitpack

class _BivariateSplineBase(object):
    """ Base class for Bivariate spline s(x,y) interpolation on the rectangle
    [xb,xe] x [yb, ye] calculated from a given set of data points
    (x,y,z).

    See Also
    --------
    bisplrep, bisplev : an older wrapping of FITPACK
    BivariateSpline :
        implementation of bivariate spline interpolation on a plane grid
    SphereBivariateSpline :
        implementation of bivariate spline interpolation on a spherical grid
    """

    # def get_residual(self):
    #     """ Return weighted sum of squared residuals of the spline
    #     approximation: sum ((w[i]*(z[i]-s(x[i],y[i])))**2,axis=0)
    #     """
    #     return self.fp

    # def get_knots(self):
    #     """ Return a tuple (tx,ty) where tx,ty contain knots positions
    #     of the spline with respect to x-, y-variable, respectively.
    #     The position of interior and additional knots are given as
    #     t[k+1:-k-1] and t[:k+1]=b, t[-k-1:]=e, respectively.
    #     """
    #     return self.tck[:2]

    def get_coeffs(self):
        """ Return spline coefficients."""
        return self.tck[2]

    # def __call__(self, x, y, mth=None, dx=0, dy=0, grid=True):
    def __call__(self, x, y):
        """
        Evaluate the spline or its derivatives at given positions.

        Parameters
        ----------
        x, y : array_like
            Input coordinates.

            If `grid` is False, evaluate the spline at points ``(x[i],
            y[i]), i=0, ..., len(x)-1``.  Standard Numpy broadcasting
            is obeyed.

            If `grid` is True: evaluate spline at the grid points
            defined by the coordinate arrays x, y. The arrays must be
            sorted to increasing order.
        dx : int
            Order of x-derivative

            .. versionadded:: 0.14.0
        dy : int
            Order of y-derivative

            .. versionadded:: 0.14.0
        grid : bool
            Whether to evaluate the results on a grid spanned by the
            input arrays, or at points specified by the input arrays.

            .. versionadded:: 0.14.0

        mth : str
            Deprecated argument. Has no effect.

        """
        x = np.asarray(x)
        y = np.asarray(y)

        # if mth is not None:
        #     warnings.warn("The `mth` argument is deprecated and will be removed",
        #                   FutureWarning)

        tx, ty, c = self.tck[:3]
        kx, ky = self.degrees
        # if grid:
        #     if x.size == 0 or y.size == 0:
        #         return np.zeros((x.size, y.size), dtype=self.tck[2].dtype)

        #     if dx or dy:
        #         z,ier = dfitpack.parder(tx,ty,c,kx,ky,dx,dy,x,y)
        #         if not ier == 0:
        #             raise ValueError("Error code returned by parder: %s" % ier)
        #     else:
        #         z,ier = dfitpack.bispev(tx,ty,c,kx,ky,x,y)
        #         if not ier == 0:
        #             raise ValueError("Error code returned by bispev: %s" % ier)
        # else:
            # standard Numpy broadcasting
        if x.shape != y.shape:
            x, y = np.broadcast_arrays(x, y)

        shape = x.shape
        x = x.ravel()
        y = y.ravel()

        if x.size == 0 or y.size == 0:
            return np.zeros(shape, dtype=self.tck[2].dtype)

        # if dx or dy:
        #     z,ier = dfitpack.pardeu(tx,ty,c,kx,ky,dx,dy,x,y)
        #     if not ier == 0:
        #         raise ValueError("Error code returned by pardeu: %s" % ier)
        # else:
        z,ier = dfitpack.bispeu(tx,ty,c,kx,ky,x,y)
        if not ier == 0:
            raise ValueError("Error code returned by bispeu: %s" % ier)

        z = z.reshape(shape)
        return z

class BivariateSpline(_BivariateSplineBase):
    """
    Base class for bivariate splines.

    This describes a spline ``s(x, y)`` of degrees ``kx`` and ``ky`` on
    the rectangle ``[xb, xe] * [yb, ye]`` calculated from a given set
    of data points ``(x, y, z)``.

    This class is meant to be subclassed, not instantiated directly.
    To construct these splines, call either `SmoothBivariateSpline` or
    `LSQBivariateSpline`.

    See Also
    --------
    UnivariateSpline : a similar class for univariate spline interpolation
    SmoothBivariateSpline :
        to create a BivariateSpline through the given points
    LSQBivariateSpline :
        to create a BivariateSpline using weighted least-squares fitting
    SphereBivariateSpline :
        bivariate spline interpolation in spherical cooridinates
    bisplrep : older wrapping of FITPACK
    bisplev : older wrapping of FITPACK

    """

    @classmethod
    def _from_tck(cls, tck):
        """Construct a spline object from given tck and degree"""
        self = cls.__new__(cls)
        if len(tck) != 5:
            raise ValueError("tck should be a 5 element tuple of tx, ty, c, kx, ky")
        self.tck = tck[:3]
        self.degrees = tck[3:]
        return self

    # def ev(self, xi, yi, dx=0, dy=0):
        """
        Evaluate the spline at points

        Returns the interpolated value at ``(xi[i], yi[i]),
        i=0,...,len(xi)-1``.

        Parameters
        ----------
        xi, yi : array_like
            Input coordinates. Standard Numpy broadcasting is obeyed.
        dx : int, optional
            Order of x-derivative

            .. versionadded:: 0.14.0
        dy : int, optional
            Order of y-derivative

            .. versionadded:: 0.14.0
        """
        # return self.__call__(xi, yi, dx=dx, dy=dy, grid=False)

    # def integral(self, xa, xb, ya, yb):
        """
        Evaluate the integral of the spline over area [xa,xb] x [ya,yb].

        Parameters
        ----------
        xa, xb : float
            The end-points of the x integration interval.
        ya, yb : float
            The end-points of the y integration interval.

        Returns
        -------
        integ : float
            The value of the resulting integral.

        """
        # tx,ty,c = self.tck[:3]
        # kx,ky = self.degrees
        # return dfitpack.dblint(tx,ty,c,kx,ky,xa,xb,ya,yb)

class RectBivariateSpline(BivariateSpline):
    """
    Bivariate spline approximation over a rectangular mesh.

    Can be used for both smoothing and interpolating data.

    Parameters
    ----------
    x,y : array_like
        1-D arrays of coordinates in strictly ascending order.
    z : array_like
        2-D array of data with shape (x.size,y.size).
    bbox : array_like, optional
        Sequence of length 4 specifying the boundary of the rectangular
        approximation domain.  By default,
        ``bbox=[min(x,tx),max(x,tx), min(y,ty),max(y,ty)]``.
    kx, ky : ints, optional
        Degrees of the bivariate spline. Default is 3.
    s : float, optional
        Positive smoothing factor defined for estimation condition:
        ``sum((w[i]*(z[i]-s(x[i], y[i])))**2, axis=0) <= s``
        Default is ``s=0``, which is for interpolation.

    See Also
    --------
    SmoothBivariateSpline : a smoothing bivariate spline for scattered data
    bisplrep : an older wrapping of FITPACK
    bisplev : an older wrapping of FITPACK
    UnivariateSpline : a similar class for univariate spline interpolation

    """

    # def __init__(self, x, y, z, bbox=[None] * 4, kx=3, ky=3, s=0):
    def __init__(self, x, y, z):
        kx=3
        ky=3
        x, y = ravel(x), ravel(y)

        if not all(diff(x) > 0.0):
            raise ValueError('x must be strictly increasing')
        if not all(diff(y) > 0.0):
            raise ValueError('y must be strictly increasing')
        if not ((x.min() == x[0]) and (x.max() == x[-1])):
            raise ValueError('x must be strictly ascending')
        if not ((y.min() == y[0]) and (y.max() == y[-1])):
            raise ValueError('y must be strictly ascending')
        if not x.size == z.shape[0]:
            raise ValueError('x dimension of z must have same number of '
                            'elements as x')
        if not y.size == z.shape[1]:
            raise ValueError('y dimension of z must have same number of '
                             'elements as y')
        z = ravel(z)
        xb, xe, yb, ye = [None]
        nx, tx, ny, ty, c, fp, ier = dfitpack.regrid_smth(x, y, z, xb, xe, yb,
                                                          ye, kx, ky, s)

        if ier not in [0, -1, -2]:
            msg = _surfit_messages.get(ier, 'ier=%s' % (ier))
            raise ValueError(msg)

        self.fp = fp
        self.tck = tx[:nx], ty[:ny], c[:(nx - kx - 1) * (ny - ky - 1)]
        self.degrees = kx, ky


