# '''
# This is the Python Implementation of the 
# splinter1D program.
# '''

# def spline(x, y, n, yp1, ypn, y2):
#     NMAX = 50
#     u = [0]* NMAX
    
#     if(yp1 > 0.99e30):
#         y2[0] = 0.0
#         u[0] = 0.0
#     else:
#         y2[0] = -0.5
#         u[0] = (3.0/(x[1] - x[0])) * (y[1] - y[0]) / ((x[1] - x[0]) - yp1 )

#     for i in range(1, n-1):
#         sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1])
#         p = sig * y2[i - 1] + 2
#         y2[i] = (sig - 1.0)/p
#         u[i] = (6.0*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1]) / (x[i]-x[i-1]))/(x[i+1]-x[i-1])-sig*u[i-1])/p

#     if(ypn > 0.99e30):
#         qn = 0
#         un = 0

#     else:
#         qn = 0.5
#         un = (3.0/(x[n-1]-x[n-2])) * (ypn-(y[n]-y[n-1])/(x[n]-x[n-1]))

#     y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0)

#     for k in range(n-2, -1, -1):
#         y2[k] = y2[k] * y2[k + 1] + u[k]

#     return y2


# def splint(xa, ya, y2a, n, x, y):

#     klo = 0
#     khi = n - 1

#     while True:

#         if (khi - klo) <= 1:
#             break
#         k = (khi + klo)/2

#         if(xa[k] > x):
#             khi = k

#         else:
#             klo = k

        
#     h = xa[khi] - xa[klo]
#     if h == 0:
#         print("Bad xa input in splint")
#         print("input xa value is 0")
#         return
    
#     a = (xa[khi] - x)/h
#     b = (x - xa[klo])/h
#     y = a * ya[klo] + b * ya[khi] + ((a**3 - a) * y2a[klo] + (b**3 - b)* y2a[khi]) * (h**2) /6
    
#     return y

    

import numpy as np

def spline(x, y, n, yp1, ypn):
    y2 = np.zeros(n)
    u = np.zeros(n-1)

    if yp1 > 0.99e30:
        y2[0] = 0.0
    else:
        y2[0] = -0.5
        u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1)

    for i in range(1, n-1):
        sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1])
        p = sig * y2[i-1] + 2.0
        y2[i] = (sig - 1.0) / p
        u[i] = (6.0 * ((y[i+1] - y[i]) / (x[i+1] - x[i]) - (y[i] - y[i-1]) / (x[i] - x[i-1])) / (x[i+1] - x[i-1]) - sig * u[i-1]) / p

    if ypn > 0.99e30:
        qn = 0.0
        un = 0.0
    else:
        qn = 0.5
        un = (3.0 / (x[n-1] - x[n-2])) * (ypn - (y[n-1] - y[n-2]) / (x[n-1] - x[n-2]))

    y2[n-1] = (un - qn * u[n-2]) / (qn * y2[n-2] + 1.0)

    # Back-substitution loop
    for k in range(n-2, -1, -1):
        y2[k] = y2[k] * y2[k+1] + u[k]

    return y2

def splint(xa, ya, y2a, n, x):
    klo = 0
    khi = n - 1

    # Binary search to find the correct interval
    while khi - klo > 1:
        k = (khi + klo) // 2
        if xa[k] > x:
            khi = k
        else:
            klo = k

    h = xa[khi] - xa[klo]
    if h == 0:
        raise ValueError("Bad xa input in splint: xa values must not be equal")

    # Compute spline interpolated value
    a = (xa[khi] - x) / h
    b = (x - xa[klo]) / h
    y = (a * ya[klo] + b * ya[khi] +
         ((a**3 - a) * y2a[klo] + (b**3 - b) * y2a[khi]) * (h**2) / 6.0)
    
    return y

