# IPF-FDMC_tools
Tools used for pre- or post-processing data for/from the full-wave code IPF-FDMC.

Note that the full commit history before 2022-01-21 can be found on 
https://github.tik.uni-stuttgart.de/
There is an internal repository located and the files included here are copies
from that repository, they will be regularly mirrored. Not all files from 
that internal repository can be found in this open repository due to license
and/or copyright reasons.

## Example how to interpolate an input grid to a grid used in IPF-FDMC
```
make_fullwave_interpolation( x_old=R_old, y_old=Z_old, z_old=ne_old,
                             x_new_lim=[ x0_fullwave, x1_fullwave ],
                             y_new_lim=[ y0_fullwave, y1_fullwave ],
                             f_0=28e9, period=32,
                             add_absorber=False, 
                             kill_negatives=True,
                             silent=silent, 
                             debug=False
                           )

```
