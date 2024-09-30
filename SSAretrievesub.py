from multi_grid import multidimgrid
from splinterp1D import spline, splint
import numpy as np


def bilint_new(BAM, nbands, ddcols, icolx, dxf, dyf, incr, numu, mxband, mxwaves, mxcols, mxline, mxsamp, mxcolssm, out):
    mxmat = 2
    mat = np.zeros((mxmat, mxmat))
    vector = np.array([dxf, dyf])
    axes = np.array([[dxf - 0.1, dyf - 0.1], [dxf + 0.1, dyf + 0.1]])
    
    icol1 = icolx

    for i in range(ddcols):
        for k in range(nbands):
            #print("numu", numu)
            mat[0, 0] = BAM[k][icol1 + 1]
            mat[0, 1] = BAM[k][icol1 + 2]
            mat[1, 0] = BAM[k][icol1 +int(numu) + 1]
            mat[1, 1] = BAM[k][icol1 +int(numu) + 2]

            multi_interp = multidimgrid(2, mat, axes, vector)
            out[i + 1][k] = multi_interp

        icol1 += incr

    return out

def rayterpolate(trdrsp, disortout, mlst, nlst, endboundsubband, pressure, xpres, 
                 nsurfpres, alb, nalb, mxcolssm, mxline, mxband, mxsamp, mxwaves, max_val, ssa_out):
    maxsm = 50
    sx = 0.0
    IOFi = np.zeros(max_val)
    z = np.zeros((maxsm, maxsm))
    t = np.zeros(maxsm)
    y2 = np.zeros(maxsm)

    nalb = int(nalb)
    nsurfpres = int(nsurfpres)
    endboundsubband = int(endboundsubband)

    # Sum xpres values to calculate sx
    sx = np.sum(xpres[:nsurfpres])

    # Loop over wavelengths for interpolation
    for i in range(endboundsubband):
        # Build z for this wavelength
        for j in range(nsurfpres):
            for k in range(nalb):
                z[j, k] = disortout[int(mlst[j]) + k, i]

        # Loop over albedos to compute the IOFi array
        for k in range(nalb):
            sy = 0.0
            b = 0.0
            st2 = 0.0
            for j in range(nsurfpres):
                sy += z[j, k]  # sum up z values
                t[j] = xpres[j] - (sx / nsurfpres)
                b += t[j] * z[j, k]  # sum up t(j) * z(j, k)
                st2 += t[j] ** 2     # sum up t * t

            # Compute b and a
            b /= st2
            a = (sy - sx * b) / nsurfpres
            IOFi[k] = pressure * b + a

        # Set up and run interpolation
        yp1 = (alb[1] - alb[0]) / (IOFi[1] - IOFi[0])  # y-prime(1)
        ypn = (alb[nalb - 1] - alb[nalb - 2]) / (IOFi[nalb - 1] - IOFi[nalb - 2])  # y-prime(n)

        # print(f"The value of IOFi is: {IOFi}")
        # print(f"The value of alb is: {alb}")
        # print(f"The value of nalb is: {nalb}")
        # print(f"The value of yp1 is: {yp1}")
        # print(f"The value of ypn is: {ypn}")
        # print(f"The value of y2 is: {y2}")

        # print(f"The value of trdrsp is: {trdrsp[i]}")
        # print(f"The value of ssa_out is: {ssa_out[i]}")

        # Use scipy's CubicSpline for interpolation
        y2 = spline(IOFi, alb, nalb, yp1, ypn)
        ssa_out[i] = splint(IOFi, alb, y2, nalb, trdrsp[i])
        #spline = CubicSpline(IOFi[:nalb], alb[:nalb], bc_type=((1, yp1), (1, ypn)))
        #ssa_out[i] = spline(trdrsp[i])
    return ssa_out

