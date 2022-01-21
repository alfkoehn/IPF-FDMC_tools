#!/usr/bin/env python
# coding=utf-8

"""
Core/heart of the interpolation usually used for grids in IPF-FDMC.

In addition to the main interpolation function (pchip_2D), a few functions
used for debugging/plotting purposes are included.
"""

__author__      = 'Alf Köhn-Seemann'
__email__       = 'koehn@igvp.uni-stuttgart.de'
__copyright__   = 'University of Stuttgart'
__license__     = 'MIT'


import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp


def debug_plot_cuts( x_old, y_old, z_old, x_new, y_new, z_new,  
                     fig_title='',
                   ):
#;{{{

    debug_print = False

    # get the colors from the default color-cycle
    col_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    #col_cycle = np.append( col_cycle, col_cycle )
    N_plotCut = 10

    # get coordinate in x_new corresponding to x_old
    x_ids_new2old   = None
    for ii in range(len(x_old)):
        if (x_old[ii] >= np.amin(x_new)) and (x_old[ii] <= np.amax(x_new)):
            # get old x-coordinate and corresponding ID in new array
            x_tmp_new_id = (np.abs(x_new - x_old[ii])).argmin()
            # store ids in array
            if x_ids_new2old is None:
                x_old_restr2new = np.array( [x_old[ii]] )
                x_ids_new2old = np.array( [x_tmp_new_id] )
            else:
                x_old_restr2new = np.append( x_old_restr2new, x_old[ii] )
                x_ids_new2old = np.append( x_ids_new2old, x_tmp_new_id )

    # get coordinates in y_new corresponding to y_old
    y_ids_new2old   = None
    for ii in range(len(y_old)):
        if (y_old[ii] >= np.amin(y_new)) and (y_old[ii] <= np.amax(y_new)):
            # get ID in new y-array corresponding to old y-array
            y_tmp_new_id = (np.abs(y_new - y_old[ii])).argmin()
            # store ids in array
            if y_ids_new2old is None:
                y_old_restr2new = np.array( [y_old[ii]] )
                y_ids_new2old = np.array( [y_tmp_new_id] )
            else:
                y_old_restr2new = np.append( y_old_restr2new, y_old[ii] )
                y_ids_new2old = np.append( y_ids_new2old, y_tmp_new_id )
    #for ii in range(len(y_old_restr2new)):
    #    print( 'ii = {0:4d}: y_old = {1:.3f}, y_new@old = {2:.3f}'.format(
    #        ii, y_old_restr2new[ii], y_new[y_ids_new2old[ii]] ))

    fig1      = plt.figure( figsize=(10,10) )
    plt_grid1 = plt.GridSpec( 2,2, hspace=.1, wspace=.1 )
    main_ax1  = fig1.add_subplot( plt_grid1[1,1], aspect='equal' )
    cuth_ax1  = fig1.add_subplot( plt_grid1[0,1], sharex=main_ax1 )
    cutv_ax1  = fig1.add_subplot( plt_grid1[1,0], sharey=main_ax1 )

    XX_new, YY_new = np.meshgrid( x_new, y_new )

    plot_interp_grid_style  = 2
    # for a nice explanation with handling grids (and shading) in pcolormesh, see
    # https://matplotlib.org/3.3.3/gallery/images_contours_and_fields/pcolormesh_grids.html
    if plot_interp_grid_style == 1:
        if np.amax(z_new.shape) > 1000:
            plt_red_fact    = 2
        if np.amax(z_new.shape) > 4000:
            plt_red_fact    = 4
        if np.amax(z_new.shape) > 6000:
            plt_red_fact    = 8
        if np.amax(z_new.shape) > 10000:
            plt_red_fact    = 12
        main_ax1.pcolormesh( XX_new[::plt_red_fact,::plt_red_fact], 
                             YY_new[::plt_red_fact,::plt_red_fact], 
                             z_new[::plt_red_fact,::plt_red_fact], 
                             shading='nearest' )
    elif plot_interp_grid_style == 2:
        XX_old, YY_old = np.meshgrid( x_old, y_old )
        main_ax1.pcolormesh( XX_old, 
                             YY_old, 
                             z_old,
                             shading='nearest' )
        main_ax1.contour( XX_new,
                          YY_new,
                          z_new,
                          levels=16, colors='black', linestyles='solid', 
                        )
    main_ax1.set_xlabel( 'x' )
    main_ax1.set_ylabel( 'y' )

    #y_old_plotCut = np.linspace( y_old[0], y_old[-1], N_plotCut, dtype=y_old.dtype )
    #y_old_plotCut = np.linspace( y_new[0], y_new[-1], N_plotCut, dtype=y_old.dtype )
    y_old_plotCut = np.linspace( y_old_restr2new[0], y_old_restr2new[-1], N_plotCut, dtype=y_old.dtype )
    for ii in range(len(y_old_plotCut)):
        # note: assumes that y_old_restr2new and y_ids_new2old have exactly same length
        # get ID in y_old-array corresponding to plot_cut-value
        y_plotCut_id = (np.abs(y_old - y_old_plotCut[ii])).argmin()
        # get ID in y_new-array corresponding to plot_cut-value
        y_new_plotCut_id    = (np.abs(y_new - y_old_plotCut[ii])).argmin()
        print( 'ii = {0:2d}, y_old_plotCut[ii] = {1:.3f}, y_old[plotCut_id] = {2:.3f}, y_new[y_new_plotCut_id] = {3:.3f}'.format(
            ii, y_old_plotCut[ii], y_old[y_plotCut_id], y_new[y_new_plotCut_id] ))
        cuth_ax1.plot( x_old, z_old[y_plotCut_id,:],
                       marker='x', linestyle='none',
                       color=col_cycle[ii],
                       label='y_old={0:5.3f}'.format(y_old_plotCut[ii])
                     )
        cuth_ax1.plot( x_new, z_new[y_new_plotCut_id,:],
                       linestyle='solid',
                       color=col_cycle[ii],
                       label='y_new={0:5.3f}'.format(y_new[ y_new_plotCut_id ])
                     )
        cuth_ax1.legend( fontsize=5 )
        cuth_ax1.set_xlabel('x')
        cuth_ax1.set_ylabel('z')

    #x_old_plotCut = np.linspace( x_old[0], x_old[-1], N_plotCut, dtype=x_old.dtype )
    #x_old_plotCut = np.linspace( x_new[0], x_new[-1], N_plotCut, dtype=y_old.dtype )
    x_old_plotCut = np.linspace( x_old_restr2new[0], x_old_restr2new[-1], N_plotCut, dtype=y_old.dtype )
    for ii in range(len(x_old_plotCut)):
        # note: assumes that x_old_restr2new and x_ids_new2old have exactly same length
        # get ID in x_old-array corresponding to plot_cut-value
        x_plotCut_id = (np.abs(x_old - x_old_plotCut[ii])).argmin()
        # get ID in x_new-array corresponding to plot_cut-value
        x_new_plotCut_id    = (np.abs(x_new - x_old_plotCut[ii])).argmin()
        print( 'ii = {0:2d}, x_old_plotCut[ii] = {1:.3f}, x_old[plotCut_id] = {2:.3f}, x_new[x_new_plotCut_id] = {3:.3f}'.format(
            ii, x_old_plotCut[ii], x_old[x_plotCut_id], x_new[x_new_plotCut_id] ))
        cutv_ax1.plot( z_old[:,x_plotCut_id], y_old, 
                       marker='x', linestyle='none',
                       color=col_cycle[ii],
                       label='x_old={0:5.3f}'.format(x_old_plotCut[ii])
                     )
        cutv_ax1.plot( z_new[:,x_new_plotCut_id], y_new, 
                       linestyle='solid',
                       color=col_cycle[ii],
                       label='x_new={0:5.3f}'.format(x_new[ x_new_plotCut_id])
                     )
        #cutv_ax1.set_xlim( 1.05*np.amax(z_old), (.0-.05*np.amax(z_old)) )
        cutv_ax1.legend( fontsize=5 )
        cutv_ax1.set_xlabel('z')
        cutv_ax1.set_ylabel('y')

    if debug_print:
        for ii in range( 0, len(y_new), 10 ):
            # get nearest old coordinates corresponding to coordinates new array
            y_tmp_id = (np.abs(y_new[ii] - y_old)).argmin()
            x_tmp_id = (np.abs(x_new[x_new_plotCut_id] - x_old)).argmin()

            print( 'ii = {0:4d}, y_new = {1:.4f}, y_old (nearest) = {2:.4f}, z_new = {3:.4f}, z_old (nearest) = {4:.4f}'.format(
                ii, y_new[ii], y_old[y_tmp_id], 
                z_new[ii,x_new_plotCut_id], z_old[y_tmp_id,x_tmp_id] ))

    if len(fig_title) > 0:
        fig1.suptitle( fig_title )

    plt.show()
#;}}}


def debug_plot_griddiffs( xArr1D, yArr1D, z_real, z_diffs=[],
                          fig_title='',
                        ):
#;{{{

    XX, YY = np.meshgrid( xArr1D, yArr1D )

    fig1      = plt.figure( figsize=(10,8) )
    plt_grid1 = plt.GridSpec( ncols=(len(z_diffs)+1), nrows=2, hspace=.1, wspace=.1 )

    ax1_zreal   = fig1.add_subplot( plt_grid1[0,0], aspect='equal' )

    ax1_zcompare1 = fig1.add_subplot( plt_grid1[0,1], aspect='equal' )
    ax1_zdiff1    = fig1.add_subplot( plt_grid1[1,1], aspect='equal' )

    if len(z_diffs) >= 2:
        ax1_zcompare2 = fig1.add_subplot( plt_grid1[0,2], aspect='equal' )
        ax1_zdiff2    = fig1.add_subplot( plt_grid1[1,2], aspect='equal' )
    if len(z_diffs) >= 3:
        ax1_zcompare3 = fig1.add_subplot( plt_grid1[0,3], aspect='equal' )
        ax1_zdiff3    = fig1.add_subplot( plt_grid1[1,3], aspect='equal' )
    if len(z_diffs) >= 4:
        ax1_zcompare4 = fig1.add_subplot( plt_grid1[0,4], aspect='equal' )
        ax1_zdiff4    = fig1.add_subplot( plt_grid1[1,4], aspect='equal' )

    # plot reference array
    ax1_zreal.pcolormesh( XX, YY, z_real, shading='nearest'  )
    ax1_zreal.set_title( 'reference data' )

    # plot 1st array to compare with reference array
    ax1_zcompare1.pcolormesh( XX, YY, z_diffs[0], shading='nearest'  )
    ax1_zcompare1.set_title( '1' )
    # calculate mean difference to reference array
    diff1_mean = np.mean( abs(z_real-z_diffs[0]) )
    # plot difference to reference array
    ax1_zdiff1.pcolormesh( XX, YY, (z_real-z_diffs[0]), shading='nearest'  )
    ax1_zdiff1.set_title( 'mean(|diff1|)={0:10.3e}'.format(diff1_mean) )

    if len(z_diffs) >= 2:
        ax1_zcompare2.pcolormesh( XX, YY, z_diffs[1], shading='nearest'  )
        ax1_zcompare2.set_title( '2' )
        diff2_mean = np.mean( abs(z_real-z_diffs[1]) )
        ax1_zdiff2.pcolormesh( XX, YY, (z_real-z_diffs[1]), shading='nearest'  )
        ax1_zdiff2.set_title( 'mean(|diff2|)={0:10.3e}'.format(diff2_mean) )
    if len(z_diffs) >= 3:
        ax1_zcompare3.pcolormesh( XX, YY, z_diffs[2], shading='nearest'  )
        ax1_zcompare3.set_title( '3' )
        diff3_mean = np.mean( abs(z_real-z_diffs[2]) )
        ax1_zdiff3.pcolormesh( XX, YY, (z_real-z_diffs[2]), shading='nearest'  )
        ax1_zdiff3.set_title( 'mean(|diff3|)={0:10.3e}'.format(diff3_mean) )
    if len(z_diffs) >= 4:
        ax1_zcompare4.pcolormesh( XX, YY, z_diffs[3], shading='nearest'  )
        ax1_zcompare4.set_title( '4' )
        diff4_mean = np.mean( abs(z_real-z_diffs[3]) )
        ax1_zdiff4.pcolormesh( XX, YY, (z_real-z_diffs[3]), shading='nearest'  )
        ax1_zdiff4.set_title( 'mean(|diff4|)={0:10.3e}'.format(diff4_mean) )

    if len(fig_title) > 0:
        fig1.suptitle( fig_title )

    plt.show()

#;}}}


def debug_plot_gridOldVsNew( x_old, y_old, z_old, x_new, y_new, z_new, 
                             fig_title='', fname_plot='', colorbar=True,
                             cmap='',
                           ):
#;{{{

    XX_old, YY_old = np.meshgrid(x_old, y_old)
    XX_new, YY_new = np.meshgrid(x_new, y_new)

    fig1      = plt.figure( figsize=(10,6) )    # (widht, heigt) in inches
    plt_grid1 = plt.GridSpec( ncols=2, nrows=1, hspace=.2, wspace=.1 )

    # set colormap
    if len(cmap) == 0:
        if np.amin(z_old) >= .0:
            my_cmap = plt.cm.YlOrRd
        else:
            my_cmap = plt.cm.bwr
    else:
        my_cmap = cmap

    # plot old data 
    ax1_old  = fig1.add_subplot( plt_grid1[0,0], aspect='equal' )
    cont_old = ax1_old.pcolormesh( XX_old, YY_old, z_old, cmap=my_cmap, shading='nearest'  )
    ax1_old.set_title( 'original data' )
    # by default, add colorbar
    if colorbar:
        cb_format = '%.2f'
        cb_old = fig1.colorbar( cont_old, format=cb_format,
                                orientation='horizontal',
                                fraction=0.046, pad=0.06,
                              )
        cb_old.ax.tick_params( direction='in' )

    # plot new data
    ax1_new  = fig1.add_subplot( plt_grid1[0,1], aspect='equal' )
    cont_new = ax1_new.pcolormesh( XX_new, YY_new, z_new, cmap=my_cmap, shading='nearest'  )
    ax1_new.set_title( 'interpolated data' )
    # by default, add colorbar
    if colorbar:
        cb_new = fig1.colorbar( cont_new, format=cb_format,
                                orientation='horizontal',
                                fraction=0.046, pad=0.06,
                              )
        cb_new.ax.tick_params( direction='in' )

    if len(fig_title) > 0:
        fig1.suptitle( fig_title )

    # force ticks to point inwards
    ax1_old.tick_params( axis='both', which='both', direction='in', top=True, right=True )
    ax1_new.tick_params( axis='both', which='both', direction='in', top=True, right=True )

    # show or save plot
    if len(fname_plot) > 0:
        plt.savefig( fname_plot, bbox_inches='tight', dpi=600 )
        print( '    plot written into file {0}'.format( fname_plot ) )
    else:
        plt.show()

#;}}}


def pchip_2D( x_old, y_old, z_old, x_new, y_new, 
              debug_plot=True,
              silent=False,
            ):
#;{{{
    """"
    runs scipy.interpolate.PchipInterpolator (which is only 1D) over a 2D array
    PCHIP: Piecewise Cubic Hermite Interpolating Polynomial
    corresponds to pchip_2D_v2 from sandbox/interpolation/interpolation.py

    use two interpolation orders:
    1) interpolate along old_rows in old_data using old_data 
         ==> new_data
       interpolate along old_cols in new_data using new_data 
         ==> new_data (overwriting intersections)
       ==> new_data_row1st
    2) interpolate along old_cols in old_data using old_data
         ==> new_data
       interpolate along old_rows in new_data using new_data
         ==> new_data (overwriting intersections)
       ==> new_data_cols1st
    return mean of those two
    """

    if not silent:
        print( "pchip_2D" )

    debug_print = debug_plot

    z_new = np.zeros( (y_new.shape[0], x_new.shape[0]), dtype=float )

    if not silent:
        print( "    x_old.dtype={0}, y_old.dtype={1}, z_old.dtype={2}".format(
                x_old.dtype, y_old.dtype, z_old.dtype ) )
        print( "    x_new.dtype={0}, y_new.dtype={1}, z_new.dtype={2}".format(
                x_new.dtype, y_new.dtype, z_new.dtype ) )
        print( "    x_old.shape={0}, y_old.shape={1}, z_old={2}".format(
                x_old.shape, y_old.shape, z_old.shape ) )
        print( "    x_new.shape={0}, y_new.shape={1}, z_new={2}".format(
                x_new.shape, y_new.shape, z_new.shape ) )

    #############################################################################
    # interpolate along rows in original data, then cols in full data
    #############################################################################
    # interpolate row
    y_ids_new2old   = None
    for ii in range(len(y_old)):
        # make sure that interpolation is only performed within new grid
        if (y_old[ii] >= np.amin(y_new)) and (y_old[ii] <= np.amax(y_new)):
            # get ID in new y-array corresponding to old y-array
            y_tmp_new_id = (np.abs(y_new - y_old[ii])).argmin()
            # store ids in array
            if y_ids_new2old is None: 
                y_old_restr2new = np.array( [y_old[ii]] )
                y_ids_new2old   = np.array( [y_tmp_new_id] )
            else:
                y_old_restr2new = np.append( y_old_restr2new, y_old[ii] )
                y_ids_new2old = np.append( y_ids_new2old, y_tmp_new_id )

            if debug_print:
                print( "    row1st: ii = {0:5d}, y_old = {1:10.3e}, y_new = {2:10.3e} (ii_new = {3:5d})".format(
                       ii, y_old[ii], y_new[y_tmp_new_id], y_tmp_new_id ) )
            # PCHIP 1D monotonic cubic interpolation
            f_interp_row = interp.PchipInterpolator( x_old, z_old[ii,:] )
            z_new[ y_tmp_new_id, : ] = f_interp_row( x_new )
    # interpolate every column
    for ii in range(len(x_new)):
        # PCHIP 1D monotonic cubic interpolation
        f_interp_col = interp.PchipInterpolator( y_old_restr2new, z_new[y_ids_new2old,ii] )
        z_new[:,ii]  = f_interp_col( y_new )
    #############################################################################

    z_row1st = np.copy(z_new)
    z_new   *= .0


    #############################################################################
    # interpolate along cols in original data, then row in full data
    #############################################################################
    # interpolate every column
    x_ids_new2old   = None
    counter         = 0
    for ii in range(len(x_old)):
        # make sure that interpolation is only performed within new grid
        if (x_old[ii] >= np.amin(x_new)) and (x_old[ii] <= np.amax(x_new)):
            # get old x-coordinate and corresponding ID in new array
            x_tmp_new_id = (np.abs(x_new - x_old[ii])).argmin()
            # store ids in array
            if x_ids_new2old is None:
                x_old_restr2new = np.array( [x_old[ii]] )
                x_ids_new2old = np.array( [x_tmp_new_id] )
            else:
                x_old_restr2new = np.append( x_old_restr2new, x_old[ii] )
                x_ids_new2old = np.append( x_ids_new2old, x_tmp_new_id )

            if debug_print:
                print( "    col1st: ii = {0:5d}, x_old = {1:10.3e}, x_new = {2:10.3e} (ii_new = {3:5d})".format(
                       ii, x_old[ii], x_new[x_tmp_new_id], x_tmp_new_id ) )
            # PCHIP 1D monotonic cubic interpolation
            f_interp_col = interp.PchipInterpolator( y_old, z_old[:,ii] )
            z_new[ :, x_tmp_new_id ] = f_interp_col( y_new )
            counter += 1
    # interpolate every row
    for ii in range(len(y_new)):
        # PCHIP 1D monotonic cubic interpolation
        f_interp_row = interp.PchipInterpolator( x_old_restr2new, z_new[ii,x_ids_new2old] )
        z_new[ii,:]  = f_interp_row( x_new )
    #############################################################################

    z_col1st = np.copy(z_new)

    z_mean = .5*(z_col1st + z_row1st)

    if not silent:
        print( '    min/max in original z-array: {0:13.4e}/{1:13.4e}'.format(
            np.amin(z_old), np.amax(z_old) ))
        print( '    min/max in new z-array     : {0:13.4e}/{1:13.4e}'.format(
            np.amin(z_new), np.amax(z_new) ))

    if debug_plot:
    
        debug_plot_cuts( x_old, y_old, z_old, x_new, y_new, z_row1st, 
                         fig_title='pchip_2D, row1st' )

        debug_plot_cuts( x_old, y_old, z_old, x_new, y_new, z_col1st, 
                         fig_title='pchip_2D, col1st' )
        
        debug_plot_cuts( x_old, y_old, z_old, x_new, y_new, z_mean, 
                         fig_title='pchip_2D, (row1st+col1st)/2' )

        debug_plot_griddiffs( x_new, y_new, 
                              z_col1st, z_diffs=[z_row1st,z_mean], 
                              fig_title='pchip_2D, col1st vs. row1st vs. mean'
                            )

        debug_plot_gridOldVsNew( x_old, y_old, z_old, x_new, y_new, z_mean, 
                                 fig_title='z_old vs. z_mean',
                               )

    return z_mean

#;}}}


