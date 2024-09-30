import multiprocessing as mp
import numpy as np
import time
import sys
from SSAretrievesub import bilint_new, rayterpolate

def process_sample(i, endboundline, endboundsubband, waves, _lambda, CRwasub, CRsbsub, BAM, dcols, icol, dxf, dyf, incr, numu, mxband, mxwaves, mxcols, mxline, mxsamp, mxcolssm, trdrsub, mlst, nlst, pressure, pres, nsurfpres, alb, nalb, max):
    GaussS = np.zeros((mxband, mxwaves))
    dsubset = np.zeros((mxcolssm, mxwaves))
    out = np.zeros((mxcolssm, mxband))
    trdrsp = np.zeros(mxband)
    xssa = np.zeros(mxband)
    ssa_sample = np.zeros((endboundline, endboundsubband))

    for k in range(waves):
        for j in range(endboundsubband):
            if abs(_lambda[k] - CRwasub[j][i]) > 35.0:
                GaussS[j][k] = 0.0
            else:
                GaussS[j][k] = CRsbsub[j][2][i] * np.exp(-1 * CRsbsub[j][3][i] * ((_lambda[k] - CRsbsub[j][4][i]) ** 2)) * 3

    for l in range(endboundline):
        icolx = int(icol[l][i])
        dsubset = bilint_new(BAM, waves, dcols, icolx, dxf[l,i], dyf[l,i], incr, numu, mxband, mxwaves, mxcols, mxline, mxsamp, mxcolssm, dsubset)

        for j in range(1, dcols + 1):
            for k in range(endboundsubband):
                out[j][k] = np.sum(GaussS[k][:] * dsubset[j][:] * (_lambda[1] - _lambda[0]))

        for j in range(endboundsubband):
            trdrsp[j] = trdrsub[l][j][i]

        xssa = rayterpolate(trdrsp, out, mlst, nlst, endboundsubband, pressure, pres, nsurfpres, alb, nalb, mxcolssm, mxline, mxband, mxsamp, mxwaves, max, xssa)
        ssa_sample[l][:endboundsubband] = xssa[:endboundsubband]

    return ssa_sample

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
        woffset = tuple(map(float, file.readline().split()))[0]
        pressure = tuple(map(float, file.readline().split()))[0]
        wmin, wmax = tuple(map(float, file.readline().split()))

    
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
    endboundparam=end_param-start_param+1


    #Read ddr and extract sub-cube based on input parameters
    ############### Paste Section ######################################
    with open(ddrfile, 'r') as file:
        next(file)  # skip line
        next(file)  # skip line
        w = next(file).split()
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
        #print("ina values:")
        count = 0
        for i in range(int(start_line), int(end_line) + 1):
            line = file.readline().split()
            #print(line)
            skiptmp = [float(line[j]) for j in range(int(end_sample))]
            for j in range(int(start_sample), int(end_sample) + 1):
                ina[i - int(start_line), j - int(start_sample)] = skiptmp[j - 1]
                #print(ina[i - int(start_line), j - int(start_sample)])
                #count += 1
        

        
        for i in range(int(end_line) + 1, int(nline) + 1):
            next(file)  # skips to end of ina
        next(file)  # skip line
        
        for _ in range(int(start_line) - 1):
            next(file)  # skips to info at start_line
        
        #print("pha values:")
        for i in range(int(start_line), int(end_line) + 1):
            line = file.readline().split()
            #print(line)
            skiptmp = [float(line[j]) for j in range(int(end_sample))]
            for j in range(int(start_sample), int(end_sample) + 1):
                ema[i - int(start_line), j - int(start_sample)] = skiptmp[j - 1]
                #print(ema[i - int(start_line), j - int(start_sample)])
                #count += 1
        
        
        
        for i in range(int(end_line) + 1, int(nline) + 1):
            next(file)  # skips to end of ema
        next(file)  # skip line
        
        for _ in range(int(start_line) - 1):
            next(file)  # skips to info at start_line
        
        #print("pha values:")
        for i in range(int(start_line), int(end_line) + 1):
            line = file.readline().split()
            #print(line)
            skiptmp = [float(line[j]) for j in range(int(end_sample))]
            for j in range(int(start_sample), int(end_sample) + 1):
                pha[i - int(start_line), j - int(start_sample)] = skiptmp[j - 1]
                #print(pha[i - int(start_line), j - int(start_sample)])
                #count += 1
        #print("count: ", count)
        #return
        
        
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
        #print("cols, waves", cols, waves) #checked till here
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
            #print("BAM:", BAM[i, :])
    
    skiptmp = np.zeros(int(end_sample))
    
    with open(crwalxfile, 'r') as file:
        for i in range(start_bandc, end_bandc + 1):
            line = file.readline().split()
            skiptmp[:] = [float(line[j]) for j in range(end_sample)]
            for j in range(start_sample, end_sample + 1):
                CRwalx[i - start_bandc, j - start_sample] = skiptmp[j - 1]
                
    #Checked till here    
    ccol = endboundsamp // 2
    
    if wmin > 0.0 and wmax > 0.0:
        wave_start = wmin - 0.04
        wave_end = wmax + 0.04
    else:
        wave_start = BAM[0, 0]
        wave_end = BAM[-1, 0]

    for i in range(end_bandc, 0, -1):
        if (wave_start + 0.04) * 1000.0 < CRwalx[i-1, ccol-1]:  # Adjusting for 0-based indexing
            subbandhi = i
            #print(wave_start + 0.04, CRwalx[i-1, ccol-1])
            break

    for i in range(1, end_bandc+1):
        if (wave_end - 0.04) * 1000.0 > CRwalx[i - 1, ccol - 1]:
            #print((wave_end - 0.04), CRwalx[i - 1, ccol - 1])
            subbandlow = i
            break

    endboundsubband=subbandhi-subbandlow+1
    #CRwasub = CRwalx[subbandlow - 1:subbandhi, :endboundsamp]

    for i in range(subbandlow, subbandhi + 1):
        for j in range(1, endboundsamp + 1):
            CRwasub[i - subbandlow, j - 1] = CRwalx[i - 1, j - 1]  # Adjust for 0-based indexing
            
    
    #CRsbsub = np.zeros((subbandhi - subbandlow + 1, end_param - start_param + 1, end_sample - start_sample + 1))
    skiptmp = np.zeros(end_sample)
    
    with open(crsblxfile, 'r') as file:
    # Skip the first part of the file
        for i in range(1, subbandlow):
            for j in range(start_param, end_param + 1):
                _ = file.readline()  # Skip a line
                if j == end_param:
                    _ = file.readline()  # Skip another line if j == end_param
        
        # Read and resample the data
        for i in range(subbandlow, subbandhi + 1):
            for j in range(start_param, end_param + 1):
                # Read the entire line of band samples into skiptmp
                skiptmp = np.array(file.readline().split(), dtype=float)
                
                # Process the specified sample range
                for k in range(start_sample, end_sample + 1):
                    CRsbsub[i - subbandlow, j - start_param, k - start_sample] = skiptmp[k - 1]
                    

                if j == end_param:
                    _ = file.readline()  # Skip another line if j == end_param

    

    subbandhi += 1
    subbandlow += 1

    # Check condition
    if (subbandhi - subbandlow + 1) != endboundsubband:
        print('Error in incrementing bands for trdr, check endboundsubband value')
        raise SystemExit  # Equivalent to Fortran's `stop`
    
    #Checked till here
    w = np.zeros(10)
    # with open(trdrfile, 'r') as file:
    #     next(file)  # skip line
    #     next(file)  # skip line
    #     w = next(file).split()
    #     nsamp = int(w[3])  # extract #of samples
    #     nline = int(w[6])  # extract #of lines
    #     nband = int(w[9])  # extract #of bands
        
    #     if nband > mxband:
    #         print(f'nband > mxband, nband, mxband = {nband}, {mxband}')
    #         print('must increase mxband')
    #         raise ValueError('Increase mxband')
        
    #     if nsamp > mxsamp:
    #         print(f'nsamp > mxsamp, nsamp, mxsamp = {nsamp}, {mxsamp}')
    #         print('must increase mxsamp')
    #         raise ValueError('Increase mxsamp')
        
    #     if nline > mxline:
    #         print(f'nline > mxline, nline, mxline = {nline}, {mxline}')
    #         print('must increase mxline')
    #         raise ValueError('Increase mxline')
        
    #     if nsamp != nsamples:
    #         print('Sample dimension of trdr file do not match dimension of parameter file')
    #         raise ValueError('Dimension mismatch')
        
    #     if nline != nlines:
    #         print('Line dimension of trdr file do not match dimension of parameter file')
    #         raise ValueError('Dimension mismatch')
        
    #     w = next(file).split()
    #     fmt_trdr = next(file).split()[0]  # extract format code
    #     next(file)  # skip line
        
    #     for _ in range(start_line - 1):
    #         for j in range(nband):
    #             next(file)  # skips to start_line
    #             if j == nband-1:
    #                 next(file)
        
    #     #trdrsub = np.zeros((end_line - start_line + 1, subbandhi - subbandlow + 1, end_sample - start_sample + 1))
    #     #skiptmp = np.zeros(end_sample)
        
    #     print("trdrsub:")
    #     print("start sample end sample", start_sample, end_sample)

    #     for i in range(start_line, end_line + 1):
    #         # Skip over unnecessary bands
    #         for j in range(1, subbandlow):
    #             _ = file.readline()  # Skipping the line

    #         # Start reading matrix
    #         for j in range(subbandlow, subbandhi + 1):
    #             # Read entire samples of band into skiptmp
    #             skiptmp = np.array(file.readline().split(), dtype=float)

    #             # Grab the specified sample range
                
    #             for k in range(start_sample, end_sample + 1):
    #                 trdrsub[i - start_line, j - subbandlow, k - start_sample] = skiptmp[k - 1]
    #                 count += 1
    #                 print(trdrsub[i - start_line, j - subbandlow, k - start_sample])

    #         # Skip bands after subbandhi
    #         for j in range(subbandhi + 1, nband + 1):
    #             _ = file.readline()  # Skipping the line

    #         # Skip additional line unless it's the last one
    #         if i != end_line:
    #             _ = file.readline()
        
    #     print("count:", count)


    # Read and resample TRDR, and extract sub-cube based on input parameters
    with open(trdrfile, 'r') as file:
        _ = file.readline()  # skip line
        _ = file.readline()  # skip line

        # Read the w array
        w = file.readline().split()
        nsamp = int(w[3])  # extract number of samples
        nline = int(w[6])  # extract number of lines
        nband = int(w[9])  # extract number of bands

        # Prevent matrix overflow
        if nband > mxband:
            raise ValueError(f'nband > mxband, nband, mxband = {nband}, {mxband}. Must increase mxband.')

        if nsamp > mxsamp:
            raise ValueError(f'nsamp > mxsamp, nsamp, mxsamp = {nsamp}, {mxsamp}. Must increase mxsamp.')

        if nline > mxline:
            raise ValueError(f'nline > mxline, nline, mxline = {nline}, {mxline}. Must increase mxline.')

        # Assure dimensions match
        if nsamp != nsamples:
            raise ValueError('Sample dimension of TRDR file does not match dimension of parameter file.')

        if nline != nlines:
            raise ValueError('Line dimension of TRDR file does not match dimension of parameter file.')

        # Skip another line
        #_ = file.readline()  
        fmt_trdr = file.readline().strip()  # extract format code
        _ = file.readline()  # skip line

        # Skip over unnecessary lines
        for i in range(1, start_line):
            for j in range(1, nband + 1):
                _ = file.readline()
                if j == nband:
                    _ = file.readline()

        

        # Read and process the TRDR data
        count = 0
        for i in range(start_line, end_line + 1):
            # Skip unnecessary bands
            for j in range(1, subbandlow):
                _ = file.readline()

            # Start reading the matrix
            for j in range(subbandlow, subbandhi + 1):
                skiptmp = np.array(file.readline().split(), dtype=float)  # reads entire samples of band

                # Grab the specified sample range
                for k in range(start_sample, end_sample + 1):
                    trdrsub[i - start_line, j - subbandlow, k - start_sample] = skiptmp[k - 1]
                    #print(trdrsub[i - start_line, j - subbandlow, k - start_sample])

            # Skip bands after subbandhi
            for j in range(subbandhi + 1, nband + 1):
                _ = file.readline()

            # Skip additional line unless it's the last one
            if i != end_line:
                _ = file.readline()

    # Close file automatically with 'with' block

    
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
    def process_pixel(i, waves, endboundsubband, endboundline, xcols, CRwasub, _lambda, CRsbsub, trdrsub, BAM, \
                  dcols, icol, dxf, dyf, incr, numu, mxband, mxwaves, mxcols, mxline, mxsamp, mxcolssm, stepsize, \
                  pressure, pres, nsurfpres, alb, nalb, mlst, nlst, max, rayterpolate):
        print(f"Calculating Loop {i}")
        GaussS = np.zeros((endboundsubband, waves))
        out = np.zeros((xcols, endboundsubband))
        trdrsp = np.zeros(endboundsubband)
        ans = []
        
        # Step 1: Loop over wavelengths
        for k in range(waves):
            for j in range(endboundsubband):
                if (abs(_lambda[k] - CRwasub[j][i])) > 35.0:
                    GaussS[j][k] = 0.0
                else:
                    GaussS[j][k] = (
                        CRsbsub[j, 1, i] * np.exp(-1 * CRsbsub[j, 2, i] * (_lambda[k] - CRsbsub[j, 3, i])**2)
                        + CRsbsub[j, 4, i] * np.exp(-1 * CRsbsub[j, 5, i] * (_lambda[k] - CRsbsub[j, 6, i])**2)
                        + CRsbsub[j, 7, i] * np.exp(-1 * CRsbsub[j, 8, i] * (_lambda[k] - CRsbsub[j, 9, i])**2)
                    )

        # Step 2: Find closest phi and cose for this pixel
        for l in range(endboundline):
            icolx = int(icol[l][i])
            dsubset = bilint_new(BAM, waves, dcols, icolx, dxf[l, i], dyf[l, i], incr, numu, mxband, mxwaves, \
                                mxcols, mxline, mxsamp, mxcolssm)
            
            # Step 3: Resample interpolated data to CRISM
            for j in range(1, xcols):
                for k in range(endboundsubband):
                    out[j][k] = 0.0
                    for m in range(waves):
                        out[j][k] += (GaussS[k][m] * dsubset[j][m] * stepsize)

            # Step 4: Extract CRISM spectrum for this location
            for j in range(endboundsubband):
                trdrsp[j] = trdrsub[l][j][i]

            # Step 5: Perform rayterpolate and save the result
            xssa = rayterpolate(trdrsp, out, mlst, nlst, endboundsubband, pressure, pres, nsurfpres, \
                                alb, nalb, mxcolssm, mxline, mxband, mxsamp, mxwaves, max)
            
            ans.append(xssa[:endboundsubband])

        return ans

    def compute_gauss(k, j, i):
        if abs(_lambda[k] - CRwasub[j][i]) > 35.0:
            return 0.0
        else:
            gauss_term = CRsbsub[j][2][i] * np.exp(-1 * CRsbsub[j][3][i] * ((_lambda[k] - CRsbsub[j][4][i]) ** 2))
            return gauss_term * 3  # Merged repeated expression
        
    
    print("Total Loops: ", endboundsamp)
    ans = []

    # Create a pool of workers
    with mp.Pool(mp.cpu_count()//2) as pool:
        results = pool.starmap(process_pixel, [(i, waves, endboundsubband, endboundline, xcols, CRwasub, _lambda, \
                                                CRsbsub, trdrsub, BAM, dcols, icol, dxf, dyf, incr, numu, mxband, \
                                                mxwaves, mxcols, mxline, mxsamp, mxcolssm, stepsize, pressure, pres, \
                                                nsurfpres, alb, nalb, mlst, nlst, max, rayterpolate) 
                                               for i in range(3)])

        # Aggregate results
        for result in results:
            ans.extend(result)
    
    print("Parallel processing complete.")
    return ans

    # ... (keep the rest of the code as is)

if __name__ == "__main__":
    BAMmain()