
import sys
import os
import numpy as np
import time
import numpy as np


def multidimgrid(ndim_interp, unit_grid, axes, vector):

    nodes = 2 ** ndim_interp
    axes_offset = axes[0, :]
    axes_scale = axes[1] - axes[0]
    unit_axes = np.zeros((2,ndim_interp))
    unit_axes[0, :] = 0.0 #initialize it
    unit_axes[1, :] = 1.0

    unit_vector = (vector-axes_offset)/axes_scale
    #print("unit vector:", vector-axes_offset, axes_scale)
    reform_grid = np.zeros(nodes)
    ii = 0

    for i in range(2):
        for j in range(2):
            reform_grid[ii]  =unit_grid[i][j]
            ii += 1


    multi_interp = 0.0
    for i in range(0, nodes):
        index_grid_value = reform_grid[i]
        interpolate_product=1.0
        index_product=1.0

        for j in range(ndim_interp):
            i1 = int((((i-1)/2**(j-1)) % 2))
            i2 = int((((i-1)*ndim_interp + (j-1)) % ndim_interp ))
            index_product = ( 1.0 - abs(unit_vector[j] - unit_axes[i1,i2]) )
            interpolate_product = interpolate_product * index_product

        index_sum = index_grid_value * interpolate_product
        multi_interp = multi_interp + index_sum
        #print("multi:", multi_interp)

    return multi_interp


def spline(x, y, n, yp1, ypn, y2):
    NMAX = 50
    u = [0]* NMAX
    
    if(yp1 > 0.99e30):
        y2[0] = 0.0
        u[0] = 0.0
    else:
        y2[0] = -0.5
        u[0] = (3.0/(x[1] - x[0])) * (y[1] - y[0]) / ((x[1] - x[0]) - yp1 )

    for i in range(1, n-1):
        sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1])
        p = sig * y2[i - 1] + 2
        y2[i] = (sig - 1.0)/p
        u[i] = (6.0*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1]) / (x[i]-x[i-1]))/(x[i+1]-x[i-1])-sig*u[i-1])/p

    if(ypn > 0.99e30):
        qn = 0
        un = 0

    else:
        qn = 0.5
        un = (3.0/(x[n-1]-x[n-2])) * (ypn-(y[n]-y[n-1])/(x[n]-x[n-1]))

    y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0)

    for k in range(n-2, -1, -1):
        y2[k] = y2[k] * y2[k + 1] + u[k]

    return y2


def splint(xa, ya, y2a, n, x, y):

    klo = 0
    khi = n - 1

    while True:

        if (khi - klo) <= 1:
            break
        k = (khi + klo)/2

        if(xa[k] > x):
            khi = k

        else:
            klo = k

        
    h = xa[khi] - xa[klo]
    if h == 0:
        print("Bad xa input in splint")
        print("input xa value is 0")
        return
    
    a = (xa[khi] - x)/h
    b = (x - xa[klo])/h
    y = a * ya[klo] + b * ya[khi] + ((a**3 - a) * y2a[klo] + (b**3 - b)* y2a[khi]) * (h**2) /6
    
    return y


def bilint_new(BAM, nbands, ddcols, icolx, dxf, dyf, incr, numu, mxband, mxwaves, mxcols, mxline, mxsamp, mxcolssm, out):
    mxmat = 2
    mat = np.zeros((mxmat, mxmat))
    vector = np.zeros(2)
    axes = np.zeros((2, 2))
    
    icol1 = icolx
    #print("icol1:", icol1)

    for i in range(ddcols):
        for k in range(nbands):
            mat[0][0] = BAM[k][icol1 + 2 - 1]
            mat[0][1] = BAM[k][icol1 + 3 - 1]
            mat[1][0] = BAM[k][icol1 + numu + 2 - 1]
            mat[1][1] = BAM[k][icol1 + numu + 3 - 1]

            axes[0][0]=dxf - 0.1
            axes[1][0]=dxf+.1
            axes[0][1]=dyf-.1
            axes[1][1]=dyf+.1
            vector[0]=dxf
            vector[1]=dyf

            multi_interp = multidimgrid(2,mat,axes,vector)
            #print(multi_interp)
            out[i+1][k]=multi_interp

        icol1 += incr

    return out

def rayterpolate(trdrsp,disortout,mlst,nlst,endboundsubband,pressure,xpres, nsurfpres,alb,nalb,mxcolssm,mxline,mxband,mxsamp,mxwaves,max,ssa_out):

    maxsm = 50
    z = np.zeros((maxsm, maxsm))
    t = np.zeros(maxsm)
    IOFi = np.zeros(max)
    y2 = np.zeros(maxsm)

    sx = 0.0
    for i in range(nsurfpres):
        sx += xpres[i]

    for i in range(endboundsubband):
        for j in range(nsurfpres):
            for k in range(nalb):
                z[j][k] = disortout[mlst[j]+k][i]

        
        for k in range(nalb):
            sy = 0.0
            b = 0.0
            st2 = 0.0

            for j in range(nsurfpres):
                sy+= z[j][k]
                t[j] = xpres[j] = (sx/nsurfpres)
                b += (t[j]*z[j][k])
                st2 += (t[j]**2)

            b = b/st2
            a = (sy-(sx*b))/nsurfpres
            IOFi[k] = pressure * b + a

        yp1 = (alb[1] - alb[0])/(IOFi[1] - IOFi[0])
        ypn = (alb[nalb] - alb[nalb - 1])/(IOFi[nalb] - IOFi[nalb - 1])# Doubt about nalb value
        

        y2 = spline(IOFi,alb,nalb,yp1,ypn,y2)
        ssa_out[i] = splint(IOFi,alb,y2,nalb,trdrsp[i],ssa_out[i])
    
    return ssa_out


def BAMmain():
    start_bandc=1
    end_bandc=437
    start_param=1
    end_param=10
    incr1=0.0
    incr2=0.0
    incr4=0.0
    incr9=0.0
    incrA=0

    DtoR=0.017453292
    RtoD=57.29578049

    max = 200
    mxsamp = 1000
    mxline = 1000
    mxband = 1000
    mxcols = 30000
    mxwaves = 3000
    mxparam = 20
    crmxbands = 500
    crmxparam = 15
    crmxsamp = 700
    mxsmcols = 5
    mxcolssm = 200

    xmu = np.zeros(max)
    _lambda = np.zeros(mxwaves)
    xphi = np.zeros(max)
    alb = np.zeros(max)
    pres = np.zeros(max)
    mlst = np.zeros(max)
    nlst = np.zeros(max)

    # Initialize arrays
    ina = np.zeros((mxline,mxsamp))
    ema = np.zeros((mxline,mxsamp))
    pha = np.zeros((mxline,mxsamp))
    skiptmp = np.zeros(mxline)
    sini = np.zeros((mxline,mxsamp))
    cosi = np.zeros((mxline,mxsamp))
    cose = np.zeros((mxline,mxsamp))
    sine = np.zeros((mxline,mxsamp))
    phaN = np.zeros((mxline,mxsamp))
    arcphi = np.zeros((mxline,mxsamp))
    phi = np.zeros((mxline,mxsamp))
    x = np.zeros((mxline,mxsamp))
    x0 = np.zeros((mxline,mxsamp))
    dxf = np.zeros((mxline,mxsamp))
    y = np.zeros((mxline,mxsamp))
    y0 = np.zeros((mxline,mxsamp))
    dyf = np.zeros((mxline,mxsamp))
    icol = np.zeros((mxline,mxsamp))
    BAM = np.zeros((mxwaves, mxcols))
    CRwalx = np.zeros((mxwaves,mxcols))
    CRwasub = np.zeros((crmxbands,crmxsamp))
    CRsbsub = np.zeros((crmxbands,crmxparam,crmxsamp))
    trdrsub = np.zeros((mxline,mxband,mxsamp))
    GaussS = np.zeros((mxband,mxwaves))
    out = np.zeros((mxcolssm, mxband))
    trdrsp = np.zeros(mxband)
    ssa = np.zeros((crmxsamp,mxline,mxband))
    xssa = np.zeros(mxband)
    zwave = np.zeros(mxband)


    
    arg_len = len(sys.argv)
    #print(sys.argv)
    if arg_len != 4:
        raise ValueError("Not enough ")

    p_file = sys.argv[1]
    ssaofile = sys.argv[2]
    hdrofile = sys.argv[3]


    print("##############################################")
    print("### SSAretriev starting, version MKFortran ###")
    print("##############################################")

    #This part was tested it works
    with open(p_file, "r") as file:
        trdrfile = file.readline()[1:-2]
        ddrfile = file.readline()[1:-2]
        disortfile = file.readline()[1:-2]
        crsblxfile = file.readline()[1:-2]
        crwalxfile = file.readline()[1:-2]
        numu, umu0, umu1, dumu = tuple(map(float, file.readline().split()))
        nphi, phi0, phi1, dphi = tuple(map(float, file.readline().split()))
        numu0, umu00, umu01, dumu0 = tuple(map(float, file.readline().split()))
        nalb, alb0, alb1, dalb = tuple(map(float, file.readline().split()))
        ntaud, taud0, taud1, dtaud = tuple(map(float, file.readline().split()))
        ntaui, taui0, taui1, dtaui = tuple(map(float, file.readline().split()))
        nozone, ozone0, ozone1, dozone = tuple(map(float, file.readline().split()))
        nsurftemp, surftemp0, surftemp1, dsurftemp = tuple(map(float, file.readline().split()))
        nsurfpres, surfpres0, surfpres1, dsurfpres = tuple(map(float, file.readline().split()))
        start_sample, end_sample, nsamples = tuple(map(int, file.readline().split()))
        start_line, end_line, nlines = tuple(map(int, file.readline().split()))
        woffset = tuple(map(float, file.readline().split()))
        pressure = tuple(map(float, file.readline().split()))
        wmin, wmax = tuple(map(float, file.readline().split()))

    numu = int(numu)
    nsurfpres = int(nsurfpres)
    #Compute vectors of phi and cose values in DISORT model
    for i in range(int(numu)):
        xmu[i] = umu0 + incr1*dumu
        incr1 += 1.0

    for i in range(int(nphi)):
        xphi[i] = phi0 + incr2*dphi
        incr2 += 1.0

    #Compute vectors of pressure and ssa values in DISORT model
    for i in range(int(nalb)):
        alb[i] = alb0 + incr4*dalb
        incr4 += 1.0

    for i in range(int(nsurfpres)):
        pres[i] = surfpres0 + incr9*nsurfpres
        incr9 += 1.0

    #Compute constants for subset BAM
    incr = int(numu * nphi)
    dcols = int(numu0*nalb*ntaud*ntaui*nozone*nsurftemp*nsurfpres)
    for i in range(int(nsurfpres)):
        mlst[i] = nalb*incrA+1
        nlst[i] = mlst[i] + (nalb - 1)
        incrA += 1.0

    endboundsamp=end_sample-start_sample+1
    endboundline=end_line-start_line+1
    #endboundparam=end_param-start_param+1

    #Read ddr and extract sub-cube based on input parameters
    ############### Paste Section ######################################
    with open(ddrfile, 'r') as file:
        next(file)  # skip line
        next(file)  # skip line
        w = next(file).split()
        #print(w[3])
        nsamp = int(w[3])  # extract #of samples
        nline = int(w[6])  # extract #of lines
        
        if nsamp > mxsamp:
            #print(f"nsamp > mxsamp, nsamp, mxsamp = {nsamp}, {mxsamp}")
            print('must increase mxsamp')
            raise ValueError('Increase mxsamp')
        
        if nline > mxline:
            print(f'nline > mxline, nline, mxline = {nline}, {mxline}')
            print('must increase mxline')
            raise ValueError('Increase mxline')
        
        if nsamp != nsamples:
            print('Sample dimension of ddr file do not match dimension of parameter file')
            raise ValueError('Dimension mismatch')
        
        if nline != nlines:
            print('Line dimension of ddr file do not match dimension of parameter file')
            raise ValueError('Dimension mismatch')
        
        
        w = next(file).split()
        #print(w)
        fmt_ddr = next(file).split()[0]  # extract format code
        next(file)  # skip line
        
        for _ in range(int(start_line) - 1):
            next(file)  # skips to info at start_line
        
        # ina = np.zeros((end_line - start_line + 1, end_sample - start_sample + 1))
        # ema = np.zeros((end_line - start_line + 1, end_sample - start_sample + 1))
        # pha = np.zeros((end_line - start_line + 1, end_sample - start_sample + 1))
        # skiptmp = np.zeros(end_sample)
        
        for i in range(int(start_line), int(end_line) + 1):
            line = file.readline().split()
            #print(line)
            skiptmp = [float(line[j]) for j in range(int(end_sample))]
            for j in range(int(start_sample), int(end_sample) + 1):
                ina[i - int(start_line), j - int(start_sample)] = skiptmp[j - 1]
        
        for i in range(int(end_line) + 1, int(nline) + 1):
            next(file)  # skips to end of ina
        next(file)  # skip line
        
        for _ in range(int(start_line) - 1):
            next(file)  # skips to info at start_line
        
        for i in range(int(start_line), int(end_line) + 1):
            line = file.readline().split()
            #print(line)
            skiptmp = [float(line[j]) for j in range(int(end_sample))]
            for j in range(int(start_sample), int(end_sample) + 1):
                ema[i - int(start_line), j - int(start_sample)] = skiptmp[j - 1]
        
        for i in range(int(end_line) + 1, int(nline) + 1):
            next(file)  # skips to end of ema
        next(file)  # skip line
        
        for _ in range(int(start_line) - 1):
            next(file)  # skips to info at start_line
        
        for i in range(int(start_line), int(end_line) + 1):
            line = file.readline().split()
            #print(line)
            skiptmp = [float(line[j]) for j in range(int(end_sample))]
            for j in range(int(start_sample), int(end_sample) + 1):
                pha[i - int(start_line), j - int(start_sample)] = skiptmp[j - 1]
        
    sini = np.sin(DtoR * ina)
    cosi = np.cos(DtoR * ina)
    cose = np.cos(DtoR * ema)
    sine = np.sin(DtoR * ema)
    phaN = np.cos(DtoR * pha)
    
    arcphi = (phaN - (cosi * cose)) / (sini * sine)
    arcphi = np.clip(arcphi, -1.0, 1.0)
    phi = 180.0 - RtoD * np.arccos(arcphi)
    
    x = (phi - phi0) / dphi
    x0 = np.floor(x).astype(int)
    x0[x0 == nphi - 1] = nphi - 2
    dxf = x - x0
    
    y = (cose - umu0) / dumu
    y0 = np.floor(y).astype(int)
    y0[y0 == numu - 1] = numu - 2
    dyf = y - y0
    
    icol = x0 * numu + y0
        

    with open(disortfile, 'r') as file:
        cols, waves = map(int, file.readline().split())
        if cols > mxcols:
            print(f'cols > mxcols; cols, mxcols = {cols}, {mxcols}')
            print('must increase mxcols')
            raise ValueError('Increase mxcols')
        
        if waves > mxwaves:
            print(f'waves > mxwaves; waves, mxwaves = {waves}, {mxwaves}')
            print('must increase mxwaves')
            raise ValueError('Increase mxwaves')
        
        BAM = np.zeros((waves, cols))
        for i in range(waves):
            BAM[i, :] = [float(x) for x in file.readline().split()]
    
    skiptmp = np.zeros(int(end_sample))
    
    with open(crwalxfile, 'r') as file:
        for i in range(start_bandc, end_bandc + 1):
            line = file.readline().split()
            skiptmp[:] = [float(line[j]) for j in range(end_sample)]
            for j in range(start_sample, end_sample + 1):
                CRwalx[i - start_bandc, j - start_sample] = skiptmp[j - 1]
        
    ccol = endboundsamp // 2
    
    if wmin > 0.0 and wmax > 0.0:
        wave_start = wmin - 0.04
        wave_end = wmax + 0.04
    else:
        wave_start = BAM[0, 0]
        wave_end = BAM[-1, 0]
    for i in range(end_bandc, 0, -1):
        if (wave_start + 0.04) * 1000.0 < CRwalx[i - 1, ccol - 1]:
            #print((wave_start + 0.04), CRwalx[i - 1, ccol - 1])
            subbandhi = i

    for i in range(1, end_bandc+1):
        if (wave_end - 0.04) * 1000.0 > CRwalx[i - 1, ccol - 1]:
            #print((wave_end - 0.04), CRwalx[i - 1, ccol - 1])
            subbandlow = i
    #print("here", (wave_end - 0.04), CRwalx[2 - 1, ccol - 1])
    #subbandhi = next(i for i in range(end_bandc, 0, -1) if (wave_start + 0.04) * 1000.0 < CRwalx[i, ccol])
    #subbandlow = next(i for i in range(1, end_bandc + 1) if (wave_end - 0.04) * 1000.0 > CRwalx[i, ccol])
    endboundsubband=subbandhi-subbandlow

    #print("subband values: ", subbandhi, subbandlow, ccol)
    subbandlow = 2
    CRwasub = CRwalx[subbandlow - 1:subbandhi, :endboundsamp]
        

    #CRsbsub = np.zeros((subbandhi - subbandlow + 1, end_param - start_param + 1, end_sample - start_sample + 1))
    skiptmp = np.zeros(end_sample)
    
    with open(crsblxfile, 'r') as file:
        for _ in range(subbandlow - 1):
            for j in range(start_param, end_param + 1):
                next(file)
                if j == end_param:
                    next(file)
        
        for i in range(subbandlow, subbandhi + 1):
            for j in range(start_param, end_param + 1):
                line = file.readline().split()
                #print(line)
                skiptmp[:] = [float(line[k]) for k in range(end_sample)]
                for k in range(start_sample, end_sample + 1):
                    CRsbsub[i - subbandlow, j - start_param, k - start_sample] = skiptmp[k - 1]
        

    w = np.zeros(10)
    with open(trdrfile, 'r') as file:
        next(file)  # skip line
        next(file)  # skip line
        w = next(file).split()
        nsamp = int(w[3])  # extract #of samples
        nline = int(w[6])  # extract #of lines
        nband = int(w[9])  # extract #of bands
        
        if nband > mxband:
            print(f'nband > mxband, nband, mxband = {nband}, {mxband}')
            print('must increase mxband')
            raise ValueError('Increase mxband')
        
        if nsamp > mxsamp:
            print(f'nsamp > mxsamp, nsamp, mxsamp = {nsamp}, {mxsamp}')
            print('must increase mxsamp')
            raise ValueError('Increase mxsamp')
        
        if nline > mxline:
            print(f'nline > mxline, nline, mxline = {nline}, {mxline}')
            print('must increase mxline')
            raise ValueError('Increase mxline')
        
        if nsamp != nsamples:
            print('Sample dimension of trdr file do not match dimension of parameter file')
            raise ValueError('Dimension mismatch')
        
        if nline != nlines:
            print('Line dimension of trdr file do not match dimension of parameter file')
            raise ValueError('Dimension mismatch')
        
        w = next(file).split()
        fmt_trdr = next(file).split()[0]  # extract format code
        next(file)  # skip line
        
        for _ in range(start_line - 1):
            for j in range(nband):
                next(file)  # skips to start_line
                if j == nband-1:
                    next(file)
        
        #trdrsub = np.zeros((end_line - start_line + 1, subbandhi - subbandlow + 1, end_sample - start_sample + 1))
        #skiptmp = np.zeros(end_sample)
        
        for i in range(start_line, end_line + 1):
            for _ in range(subbandlow - 1):
                next(file)  # skips to start_band
            
            for j in range(subbandlow, subbandhi + 1):
                line = file.readline().split()
                skiptmp = [float(line[k]) for k in range(end_sample)]
                for k in range(start_sample, end_sample + 1):
                    trdrsub[i - start_line, j - subbandlow, k - start_sample] = skiptmp[k - 1]
            
            for _ in range(subbandhi + 1, nband + 1):
                next(file)
            
            if i != end_line:
                next(file)  # skips to next line
        
    ####################################################################
    dsubset = np.zeros((mxcolssm,mxwaves))
    _lambda = np.zeros(mxwaves)
    #Create output and working arrays
    for i in range(waves):
        for j in range(cols):
            dsubset[0][i] = BAM[i][0]

        _lambda[i] = BAM[i][0] * (1000.0) - woffset

    stepsize = (BAM[1][0] - BAM[0][0]) * 1000.0
    xcols = dcols+1
    print("Total Loops: ", endboundsamp)
    
    for i in range(1):#endboundsamp):
        print("Current i", i)
        start_time = time.time()
        for k in range(waves):
            for j in range(endboundsubband):
                if (abs(_lambda[k] - CRwasub[j][i])) > 35.0:
                    GaussS[j][k] = 0.0
                
                else:
                    GaussS[j][k] = (CRsbsub[j][2][i] * np.exp( -1 * CRsbsub[j][3][i]* ((_lambda[k]-CRsbsub[j][4][i])**2)))+ \
                                    (CRsbsub[j][2][i] * np.exp( -1 * CRsbsub[j][3][i]* ((_lambda[k]-CRsbsub[j][4][i])**2)))+ \
                                    (CRsbsub[j][2][i] * np.exp( -1 * CRsbsub[j][3][i]* ((_lambda[k]-CRsbsub[j][4][i])**2)))

        #Find closest phi and cose for this pixel
        for l in range(endboundline):
            icolx = int(icol[l][i])
            dsubset = bilint_new(BAM,waves,dcols,icolx,dxf[l,i],dyf[l,i],incr,numu,mxband,mxwaves,\
                     mxcols,mxline,mxsamp,mxcolssm,dsubset)
            
            #Resample interpolated data to CRISM
            for j in range(1, xcols):
                for k in range(endboundsubband):
                    out[j][k] = 0.0
                    for m in range(waves):
                        out[j][k] += (GaussS[k][m] * dsubset[j][m] * stepsize)

            #Extract CRISM spectrum for this location
            for j in range(endboundsubband):
                trdrsp[j] = trdrsub[l][j][i]

            xssa = rayterpolate(trdrsp,out,mlst,nlst,endboundsubband,pressure,pres,nsurfpres, \
                        alb,nalb,mxcolssm,mxline,mxband,mxsamp,mxwaves,max,xssa)
            
            ssa[i][l][:endboundsubband] = xssa[j]
            
            # for j in range(endboundsubband):
            #     ssa[i][l][j] = xssa[j]
            #     print(xssa[j])
            #     print("[INFO]: SAVING NPY FILE.....")
                #np.save("SSA_Numpy.npy", ssa)

                


        
        print("Placement something") #Recheck
        print("Time taken: ", time.time() - start_time)

    np.save("SSA_Numpy.npy", ssa)
    return
    #Compute parameters for beginning of SSA output file
    lcount = 5
    for i in range(endboundsubband):
        for j in range(endboundline):
            lcount += 1
        lcount += 1

    totalline = lcount - endboundsubband - 5


    #Save SSA results to a text file that is compatible with ENVI
    with open(ssaofile, 'w') as f:
        f.write(';\n')
        f.write(';This file contains the SSA output calculations\n')
        f.write(f' ;File dimensions={endboundsamp:5} samples x{endboundline:5} lines x{endboundsubband:5} bands\n')
        f.write(f' ;Total file lines={lcount:7}    Total readable data lines={totalline:7}\n')
        f.write(';\n')
        
        for i in range(endboundsubband):
            for j in range(endboundline):
                # Write each row of data
                data_row = ssa[:, j, i]
                formatted_row = ' '.join(f'{value:37.6f}' for value in data_row)
                print("Data row", data_row)
                print("formatted row:", formatted_row)
                return
                f.write(f'{formatted_row}\n')
            f.write('\n')  # Blank line after each subband
    
    #Create ASCII text file that will serve as a header file for labeling the bands.
    with open(hdrofile, 'w') as f:
        if endboundsubband > mxband:
            print("must increase mxband for function zwave")
        for i in range( endboundsubband):
            zwave[i] = (CRwasub[i][ccol+1] + woffset) / 1000.0
            f.write(f'Band {i:4d} ({zwave[i]:6.4f})\n')

    

if __name__ == "__main__":
    BAMmain()


