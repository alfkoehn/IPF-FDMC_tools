#!/usr/bin/env python
# coding=utf-8

"""
Wrapper for the interpolation function found in myInterpolation tailored
for computational grids used in IPF-FDMC.
"""

__author__      = 'Alf Köhn-Seemann'
__email__       = 'koehn@igvp.uni-stuttgart.de'
__copyright__   = 'University of Stuttgart'
__license__     = 'MIT'


import numpy as np
import scipy.interpolate as interp
import myInterpolation

def get_new_spatial_axis( 
                          x0, x1, 
                          l_0=299792458./2.45e9, period=128,
                          add_absorber=False,
                          silent=False,
                        ):
#{{{
    """
    Create spatial axis for full-wave grid used in IPF-FDMC.

    This function simply returns an evenly spaced 1D array corresponding to
    one of the spatial axis of the computational grid used in IPF-FDMC. 

    Parameters
    ----------
    x0, x1 : float
        `x0` and `x1` are the starting and end point of the axis un units of
        meters.
    l_0 : float, optional
        Vacuum wavelength of the microwave simulated in IPF-FDMC in units of
        meters.
    period : int (or float)
        Number of grid points per vacuum wavelength. Note that this 
        corresponds to the actually used grid, i.e. in units of Yee cells it
        is period/2.
    add_absorber : bool, optional
        Optionally add absorbers of the default width to the beginning and to
        the end of the axis. Warning: the default value is hard-coded here
        which is a bad idea and should be changed.
    silent : bool, optional
        Optionally print some output to the console while running.

    Returns
    -------
    1D numpy array
    """

    # get size of new axis
    N_x = int( round( (abs(x1-x0)/l_0 * period), 0 ) )

    # optionally, add absorber at beginning and at end
    if add_absorber:
        # default value
        d_absorb = 3*period

        N_x += int(2*d_absorb)
        x0  += -d_absorb/period * l_0
        x1  +=  d_absorb/period * l_0

    # N_x has to be even numbers
    if ( N_x%2 > 0 ):
        N_x += -1
    
    x_new = np.linspace( x0, x1, N_x )

    return x_new

#}}}


def make_fullwave_interpolation( x_old, y_old, z_old, 
                                 x_new=[], y_new=[],
                                 x_new_lim=[], y_new_lim=[],
                                 f_0=2.45e9, period=128,
                                 add_absorber=False, kill_negatives=False,
                                 suppress_extrapolation=False,
                                 fill_value=[],
                                 silent=False, debug=False
                               ):
#;{{{

    errValue = -1

    if not silent:
        print( 'make_fullwave_interpolation' )

    if add_absorber:
        print( "    Note: 'add_absorber' is set, this can be realized in 3 ways" )
        print( "          1) duplicate last column or row" )
        print( "          2) smoothly set everything to zero" )
        print( "          3) include absorber-region in interpolation" )

        include_absorbers_in_interpolation = False
        if include_absorbers_in_interpolation:
            print( "        --> set to option 3" )
        else:
            print( "        --> set to option 1 (2 needs to be implemented)" )
    else:
        # if no absorbers are added, this flag also has to be set to False
        include_absorbers_in_interpolation = False


    # calculate vacuum wavelength in meters
    c_0 = 299792458.
    l_0 = c_0 / f_0

    # get sizes of old array
    N_x_old = np.size(x_old)
    N_y_old = np.size(y_old)

    # get sizes of new arrays (by creating spatial axis)
    # N_x
    if np.size(x_new) != 2:
        if np.size( x_new_lim ) != 2:
            print( "ERROR: if 'x_new' is not set, 'x_new_lim' must be set containing 2 elements (x0, x1)" )
            print( "       make_fullwave_interpolation will return -1" )
            print( np.size(x_new_lim) )
            return errValue
        else:
            x_new = get_new_spatial_axis( 
                                          x_new_lim[0], x_new_lim[1], 
                                          l_0=l_0, period=period,  
                                          add_absorber=include_absorbers_in_interpolation,
                                          silent=silent,
                                        )
    N_x = np.size(x_new)
    # N_y
    if np.size(y_new) != 2:
        if np.size( y_new_lim ) != 2:
            print( "ERROR: if 'y_new' is not set, 'y_new_lim' must be set containing 2 elements (y0, y1)" )
            print( "       make_fullwave_interpolation will return -1" )
            return errValue
        else:
            y_new = get_new_spatial_axis( 
                                          y_new_lim[0], y_new_lim[1], 
                                          l_0=l_0, period=period, 
                                          add_absorber=include_absorbers_in_interpolation,
                                          silent=silent,
                                        )
    N_y = np.size(y_new)

    if not silent:
        print( "    input: x_old = {0:13.6e} ... {1:13.6e}".format( x_old[0], x_old[-1] ) )
        print( "           y_old = {0:13.6e} ... {1:13.6e}".format( y_old[0], y_old[-1] ) )
        print( "           N_x_old = {0:5d}, N_y_old = {1:5d}".format( N_x_old, N_y_old ) )
        print( "             ==> dx = {0:13.6e}, dy = {1:13.6e}".format( 
               (x_old[-1]-x_old[0])/(N_x_old-1), (y_old[-1]-y_old[0])/(N_y_old-1) ))
        print( "                 following should be the same (otherwise spatial axis not regular)" )
        print( "                 dx = {0:13.6e}, dy = {1:13.6e}".format( 
               (x_old[1]-x_old[0]), (y_old[1]-y_old[0]) ))
        print( "           f_0 = {0:13.6e} Hz (lambda_0 = {1:13.6e})".format( f_0, l_0 ) )
        print( "           min | max of data2interpolate: {0:13.6e} | {1:13.6e}".format( 
               np.amin(z_old),np.amax(z_old) ) )
        print( "    fw-grid: x_new = {0:13.6e} ... {1:13.6e}".format( x_new[0], x_new[-1] ) )
        print( "             y_new = {0:13.6e} ... {1:13.6e}".format( y_new[0], y_new[-1] ) )
        print( "             N_x_new = {0:5d}, N_y_new = {1:5d}".format( N_x, N_y ) )
        print( "             ==> dx = {0:13.6e}, dy = {1:13.6e}".format( 
               (x_new[-1]-x_new[0])/(N_x-1), (y_new[-1]-y_new[0])/(N_y-1) ))
        print( "                 following should be the same (otherwise spatial axis not regular)" )
        print( "                 dx = {0:13.6e}, dy = {1:13.6e}".format( 
               (x_new[1]-x_new[0]), (y_new[1]-y_new[0]) ))
 
        if include_absorbers_in_interpolation:
            print( "             (absorbers are included in above values for fw-grid)" )
        elif add_absorber:
            print( "             (absorbers will be added after interpolation)" )


    # function for interpolation
    interp_method = 4
    if interp_method == 1:
        f_interp = interp.interp2d( x_old, y_old, z_old, 
                                         kind='linear',           # 'linear', 'cubic', 'quintic'
                                       )
        # do the interpolation
        z_new = f_interp( x_new, y_new )
    elif interp_method == 2:
        # radial basis function interpolation, requires all array to be of same size
        xx_old, yy_old = np.meshgrid( x_old, y_old )
        f_interp = interp.Rbf( xx_old, yy_old, z_old, 
                               function='cubic',   # 'multiquadric' (default), 'inverse', 'gaussian', 'linear', 'cubic', 'quintic', 'thin_plate'
                               smooth=0,    # 0:interpolation (default)|greater than 0: increase smoothness (allowing to deveate from original points (?))
                             )
        xx_new, yy_new = np.meshgrid( x_new, y_new )
        # do the interpolation
        z_new = f_interp( xx_new, yy_new )
    elif interp_method == 3:
        xx_old, yy_old = np.meshgrid( x_old, y_old )
        xx_new, yy_new = np.meshgrid( x_new, y_new )
        print( "AAAAAAA" )  
        print( np.min(z_old), np.amin(z_old) )
        print( len(np.where( z_old == np.amin(z_old) )) )
        print( len(np.argwhere( z_old == np.amin(z_old) )) )
        z_old_tmp = np.copy(z_old)
        z_old_tmp[z_old_tmp==np.amin(z_old)] = np.nan
        print( np.isnan(z_old_tmp).sum() )
        print( z_old_tmp[20,:] )
        print( "AAAAAAA" )  
        z_new = interp.griddata( np.array([xx_old.ravel(), yy_old.ravel()]).T, 
                                 z_old_tmp.ravel(),
                                 (xx_new, yy_new),
                                 method='cubic',   # 'linear' (default), 'cubic'
                                 fill_value=np.amin(z_old),
                               )
        print( z_new[300,:] )
        print( np.isnan(z_new).sum() )
    elif interp_method == 4:
        z_new = myInterpolation.pchip_2D( x_old, y_old, z_old, x_new, y_new, debug_plot=debug )

    if (len(fill_value) > 0) or suppress_extrapolation:
        print( "    fill_value provided during function call or suppress_extrapolation=True" )
        print( "    => will fill area where no original data is provided" )
        if x_new[-1] > x_old[-1]:
            print( "    extrapolation for large R-values detected, will be overwritten" )
            x_extrapol_1_id = np.argmin( np.abs(x_old[-1] - x_new) )

            z_new[:,x_extrapol_1_id:]   = np.amax( z_old )


    # optionally, remove negative values in interpolated array
    if (kill_negatives) and (np.amin(z_new) < .0):
        print( "ATTENTION: negative values in interpolated array" )
        print( "           will be set to zero..." )
        z_new[ z_new<.0 ] = .0
        print( "           ...done; minimum value in interpolated array: {0:13.6e}".format(np.amin(z_new)) )

    # optionally, add absobers at every side
    if add_absorber and (not include_absorbers_in_interpolation):
        if not silent: 
            print( "    starting to add absorbers (z_new.shape = {0})...".format(z_new.shape) )

        # default size of absorbers 
        # (changing this values requires changing of 'eco' in full-wave code)
        d_absorb = int( round(3.*period) )

        # by default, set top and bottom absorbers as same values as last row of non-absorber grid
        absorber_top    = z_new[ -1, : ]
        absorber_bottom = z_new[ 0,  : ]

        if not silent:
            print( "        absorber_top.shape = {0}, absorber_bottom.shape = {1}, z_new.shape = {2}".format( 
                   absorber_top.shape, absorber_bottom.shape, z_new.shape )
                 )

        # add absorbers row by row
        for ii in range(d_absorb):
            z_new = np.vstack( (absorber_bottom, z_new, absorber_top) )

        # by default, left and right absorbers have same values as last column of non-absorber grid
        # note: np.hstack cannot concatenate two arrays with different numbers of rows
        #       which is why we need to reshape it first (correspond to a transpose)
        absorber_right = z_new[ :,-1 ].reshape( (z_new.shape)[0], 1 )
        absorber_left  = z_new[ :,0  ].reshape( (z_new.shape)[0], 1 )

        if not silent:
            print( "        absorber_right.shape = {0}, absorber_left.shape = {1}, z_new.shape = {2}".format( 
                   absorber_right.shape, absorber_left.shape, z_new.shape )
                 )

        # add absorbers line by line
        for ii in range(d_absorb):
            z_new = np.hstack( (absorber_left, z_new, absorber_right) )

        if not silent:
            print( "    ...finished with adding absorbers (z_new.shape = {0}".format(z_new.shape) )

    if not silent:
        print( "    finished with interpolation" )
        print( "      shape of interpolated grid: {0}".format( z_new.shape ) )
        print( "      min = {0:13.6e} (old = {1:13.6e})".format( np.amin(z_new), np.amin(z_old) ) )
        print( "      max = {0:13.6e} (old = {1:13.6e})".format( np.amax(z_new), np.amax(z_old) ) )

    return z_new
#;}}}




