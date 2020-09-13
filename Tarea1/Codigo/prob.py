import numpy as np
import scipy.special as sci

k = 4.495

def pto_medio(func, a, b):
    dx = b-a
    med = (a+b)/2
    rec = func(med)
    return rec*dx

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
    return sci.gamma(x) #mientras se me ocurre el cambio de variable apropiado

def chi2(x):
    res = 1/(np.exp2(k/2)*Gamma(k/2))*x**(k/2 - 1)*np.exp(-x/2) #definicion de chi^2(x)
    return res

#calculo de int_0^a chi^2(x) dx, x<0 no es parte del dominio de chi^2
def prob_chi2(a, N=1000):
    prob = simpson(chi2, 0, a, N)
    return prob

#el problema prob_chi2(a) = 0.95 se puede reescribir como el problema del cero de una funcion
# definiendo g(a) = prob_chi2(a) - 0.95 = 0
def biseccion(func, a, b, tol=1e-7, lim=1000):
    p = (a+b)/2
    counter = 0
    while np.abs(func(b)-func(a)) > tol:
        p = (a+b)/2
        if counter > lim: 
            print('limite de iteraciones excedido')
            return p
        prod = func(p)*func(a)
        
        if prod > 0:
            a = p
        elif prod < 0:
            b = p
        counter += 1
    return p

def prob_menos_95(a):
    res = prob_chi2(a) - 0.95
    return res

a = biseccion(prob_menos_95, 5, 15) #el valor debe ser cercano a 10.3