import numpy as np
k = 4.495

def pto_medio():
    return

def simpson(func, a, b, N):
    dx = (b-a)/(2*N)
    f0 = func(a)
    par = 0
    impar = 0
    f2N = func(b)

    for i in range(0,N):
        par += func(a+2*i*dx)
        impar += func(a+(2*i+1)*dx)

    sum = f0 + 2*par + 4*impar + f2N
    sum = sum*dx/3
    return sum

def Gamma(x):
    return

def chi2(x):
    res = 1/(np.exp2(k/2)*Gamma(k/2))*x**(k/2 - 1)*np.exp(-x/2) #definicion de chi^2(x)
    return res