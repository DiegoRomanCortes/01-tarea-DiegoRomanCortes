import numpy as np

def simpson(func, z, a, b, N = 100, rel_tol = 1e-6):
    def trapecio(func, z, a, b, N = 100):
        dx = (b-a)/(N)
        f0 = func(z, a)
        fN = func(z, b)

        sum = (f0 + fN)/2
        for i in range(1,N):
            sum += func(z, a+i*dx)
        sum *= dx
        return sum

    dif = np.infty #primera vez siempre entra al bucle
    
    while dif > rel_tol:
        N *= 2
        sum = trapecio(func, z, a, b, N)    
        sum2 = trapecio(func, z, a, b, 2*N)
        simp = (4*sum2 - sum)/3
        dif = np.abs((simp - sum)/sum)
    #se sale del ciclo porque dif < rel_tol
    return simp

def Gamma_du(z, u): #integrando de la funcion gamma con u=e^(-x)
    gdu = np.log(1/u)**(z-1)
    return gdu

def Gamma(z, dx = 1e-6):
    def pto_medio(func, z, a, b):
        dx = b-a
        med = (a+b)/2
        rec = func(z, med)
        return rec*dx
    #nuevos extremos gracias al C.V.
    a = 0
    b = 1
    g = simpson(Gamma_du, z, a + dx, b, rel_tol=dx)
    g += pto_medio(Gamma_du, z, a, a + dx)
    return g

K = 4.495 #rut 20299495-4
GAMMA_K_MEDIOS = Gamma(K/2, dx=1e-8)

def chi2(k, x):
    res = 1/(np.exp2(k/2)*GAMMA_K_MEDIOS)*x**(k/2 - 1)*np.exp(-x/2) #definicion de chi^2(x)
    return res

#calculo de int_0^a chi^2(x) dx, x<0 no es parte del dominio de chi^2
def prob_chi2(k, a, rel_tol=1e-6):
    prob = simpson(chi2, k, 0, a, rel_tol=rel_tol)
    return prob

#el problema prob_chi2(a) = 0.95 se puede reescribir como el problema del cero de una funcion
# definiendo g(a) = prob_chi2(a) - 0.95 = 0
def biseccion(func, k, a, b, tol_abs=1e-6, lim=1e6): #lim permite evitar que el programa colapse 
    counter = 0
    dif = np.abs(func(k, b)-func(k, a))
    while  dif > tol_abs:
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

def prob_menos_95(k, a, rel_tol=1e-8):
    res = prob_chi2(k, a, rel_tol=rel_tol) - 0.95 #se quiere encontrar res = 0
    return res

def newton(func, derivada, k, x_0, tol_abs = 1e-6):
    #f(x) = \int_0^x chi^2(t)dt - 0.95
    #f'(x) = chi^2(x)
    #x_{i+1} = x_i - f(x_i)/f'(x_i)
    dif = np.infty
    while dif > tol_abs:
        x_1 = x_0 - func(k, x_0)/derivada(k, x_0)
        dif = np.abs(func(k, x_1)-func(k, x_0))
        x_0 = x_1
    return x_1

tol = 1e-5
a1 = biseccion(prob_menos_95, K, 0, 20, tol_abs=tol)
a2 = newton(prob_menos_95, chi2, K, 1, tol_abs=tol)
print('a por Biseccion:', a1)
print('a por Newton:', a2)
print('Diferencia:', np.abs(a1-a2))
a3 = newton(prob_menos_95, chi2, K, a1, tol_abs=1e-12)
print('a por Newton dando como punto inicial el de biseccion: ', a3)