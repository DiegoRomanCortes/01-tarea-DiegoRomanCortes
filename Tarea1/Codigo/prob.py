import numpy as np

def pto_medio(func, z, a, b):
    dx = b-a
    med = (a+b)/2
    rec = func(z, med)
    return rec*dx

def trapecio(func, z, a, b, N = 100):
    dx = (b-a)/(N)
    f0 = func(z, a)
    fN = func(z, b)

    sum = (f0 + fN)/2
    for i in range(1,N):
        sum += func(z, a+i*dx)
    sum *= dx
    return sum


def simpson(func, z, a, b, N = 100, rel_tol = 1e-6):
    sum = trapecio(func, z, a, b, N)    
    sum2 = trapecio(func, z, a, b, 2*N)
    simp = (4*sum2 - sum)/3
    dif = np.abs((simp - sum)/sum)
    if dif < rel_tol:
        return sum
    else:
        return simpson(func, z, a, b, 2*N, rel_tol)

def Gamma_du(z, u): #integando de la funcion gamma con u=e^(-x)
    gdu = np.log(1/u)**(z-1)
    return gdu

def Gamma(z, dx = 1e-7):
    #nuevos extremos gracias al C.V.
    a = 0
    b = 1
    g = simpson(Gamma_du, z, a + dx, b)
    g += pto_medio(Gamma_du, z, a, a + dx) 
    return g

K = 4.495 #rut 20299495-4
gamma_k_medios = Gamma(K/2)

def chi2(k, x):
    res = 1/(np.exp2(k/2)*gamma_k_medios)*x**(k/2 - 1)*np.exp(-x/2) #definicion de chi^2(x)
    return res

#calculo de int_0^a chi^2(x) dx, x<0 no es parte del dominio de chi^2
def prob_chi2(k, a):
    prob = simpson(chi2, k, 0, a)
    return prob

#el problema prob_chi2(a) = 0.95 se puede reescribir como el problema del cero de una funcion
# definiendo g(a) = prob_chi2(a) - 0.95 = 0
def biseccion(func, k, a, b, tol=1e-6, lim=1000):
    counter = 0
    dif = np.abs(func(k, b)-func(k, a))
    while  dif > tol:
        dif = np.abs(func(k, b)-func(k, a))
        p = (a+b)/2
        if counter > lim: 
            print('limite de iteraciones excedido')
            return p
        fp = func(k, p)
        fa = func(k, a)
        prod = fp*fa
        
        if prod > 0:
            a = p
        elif prod < 0:
            b = p
        counter += 1
    return p

def prob_menos_95(k, a):
    res = prob_chi2(k, a) - 0.95 #se quiere encontrar res = 0
    return res

#print(Gamma(2))
#print(Gamma(3)) #Gamma(n) = (n-1)!
                #como k/2 es cercano a 2 se han ajustado los paramentros para que gamma(2) aprox 1

a = biseccion(prob_menos_95, K, 5, 15) #a debiera ser cercano a 10.28
print("El valor de a es:", a)