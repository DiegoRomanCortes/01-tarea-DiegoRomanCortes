import numpy as np
import matplotlib.pyplot as plt

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
    #nuevos extremos gracias al C.V.
    a = 0
    b = 1
    g = simpson(Gamma_du, z, a + dx, b, rel_tol=dx)
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
def biseccion(func, k, a, b, tol=1e-6, lim=1e6): #lim permite evitar que el programa colapse 
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

#codigo para graficar el integrando original de Gamma
""" x = np.linspace(0, 10, num=500)
y = np.zeros(len(x))
counter = 0
for xi in x:
    y[counter] = xi**(K/2-1)*np.exp(-xi)
    counter += 1
plt.figure(1)
plt.clf()
plt.plot(x, y)
plt.fill_between(x, y)
plt.xlabel('x')
plt.ylabel('$y = x^{k/2 - 1}e^{-x}$', fontsize='x-large')
plt.show() """

#codigo para graficar el integrando modificado de Gamma
""" x = np.linspace(0, 1)
y = np.zeros(len(x))
counter = 0
for xi in x:
    y[counter] = Gamma_du(K/2, xi)
    counter += 1
plt.figure(2)
plt.clf()
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel(r'$y = (-\ln (x))^{K/2 - 1}$', fontsize='x-large')
plt.fill_between(x, y)
plt.show() """

#codigo para graficar Gamma
""" x = np.linspace(1, 3)
y = np.zeros(len(x))
counter = 0
for xi in x:
    y[counter] = Gamma(xi, dx = 1e-3)
    counter += 1
plt.clf()
plt.plot(x, y)
plt.show()
 """

#codigo para graficar distintas tolerancias relativas y determinar el valor a utilizar 
""" tolerancias_de_prueba = np.logspace(-6, -1)
error = np.zeros(len(tolerancias_de_prueba))

plt.clf()

counter = 0
for tol in tolerancias_de_prueba:
    dif1 = np.abs(Gamma(1, dx = tol) - 1) #Gamma(1) = 1
    dif2 = np.abs(Gamma(2, dx = tol) - 1) #Gamma(2) = 1
    dif3 = np.abs(Gamma(3, dx = tol) - 2) #Gamma(3) = 2
    dif4 = np.abs(Gamma(4, dx = tol) - 6) #Gamma(4) = 6
    error[counter] = np.mean(np.array([dif1, dif2, dif3, dif4]))
    counter += 1
plt.plot(tolerancias_de_prueba, error)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('tolerancia relativa para Simpson')
plt.ylabel('promedio de errores absolutos cometido')
plt.show() """

""" tolerancias_de_prueba_2 = np.logspace(-12, -6)
error_2 = np.zeros(len(tolerancias_de_prueba_2))

plt.clf()

counter = 0
for tol_2 in tolerancias_de_prueba_2:
    a = biseccion(prob_menos_95, K, 5, 15, tol=tol_2)
    dif = np.abs(prob_menos_95(K, a))
    error_2[counter] = dif
    counter += 1

plt.plot(tolerancias_de_prueba_2, error_2)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('tolerancia absoluta para bisección')
plt.ylabel('error absoluto cometido')
plt.show() """

#a = biseccion(prob_menos_95, K, 5, 15, tol=1e-12) #a debiera ser cercano a 10.28

#finalmente se graficará la forma de la función chi2, así como el valor a obtenido
""" x = np.linspace(0, 1.1*a)
y = chi2(K ,x)

x1 = np.linspace(0, a)
y1 = chi2(K, x1)

plt.clf()
plt.title(r'$\chi ^2 (x)$ para $k=4.495$')
plt.xlabel('x')
plt.ylabel(r'y = $\chi ^2 (x)$')
plt.axvline(a, label = 'a = {}'.format(a), color = 'r')
plt.fill_between(x1, 0, y1, label = r'$P(x < a) = \int_0^a \chi^2 (x) dx$')
plt.plot(x, y, label=r'$\chi^2(x)$')
plt.legend()
plt.show() """