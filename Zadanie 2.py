# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 19:18:53 2023

@author: wiecz
"""

import numpy
from scipy.interpolate import CubicSpline
import matplotlib as plt

x = numpy.linspace(0.0, 10.0, 100)
data = numpy.array([[1.0,3.0], [2.0,1.0], [3.5,4.0], [5.0,0.0], [6.0,0.5], [9.0,-2.0], [9.5,-3.0]])
iksy = [1.0, 2.0, 3.5, 5.0, 6.0, 9.0, 9.5]
igreki = [3.0, 1.0, 4.0, 0.0, 0.5, -2.0, -3.0]
n = data.shape[0]-1

def Funkcja_szescienna(x,data):
    #obliczanie h
    h = numpy.zeros(n)
    for i in range(n):
        h[i] = (data[i+1,0] - data[i,0])
    
    #obliczanie b
    b = numpy.zeros(n)   
    for i in range(n):
        b[i] = 6/h[i]*(data[i+1,1]-data[i,1])
    
    #obliczanie u
    u = numpy.zeros(n-1)
    u[0] = 2*(h[0]+h[1])
    for i in range(2,n):
            u[i-1] += (2*(h[i-1]+h[i])) - ((h[i-1]**2)/u[i-2])  
    
    #obliczanie v
    v = numpy.zeros(n-1)
    v[0] = b[1]-b[0]
    for i in range(2,n):
        v[i-1] += b[i]-b[i-1] - (h[i-1]*(v[i-2]/u[i-2]))
      
    #obliczanie z
    z = numpy.zeros(n+1)
    i = n
    while i > -1:
        if i == n:
            pass
        if i > 0 and i < n:
            z[i] = (v[i-1] - (h[i]*z[i+1]))/u[i-1]
        if i == 0:
            pass
        
        i -=1
    
    #obliczanie A
    A = numpy.zeros(n) 
    for i in range(n):
        A[i] += (1/(6*h[i])) * (z[i+1] - z[i])
    
    #oblilczanie B
    B = numpy.zeros(n)
    for i in range(n-1):
        if z[i] == 0:
            B[i] = z[i]
        else:
            B[i] = z[i]/2
        
    #obliczanie C
    C = numpy.zeros(n)  
    for i in range(n):
        C[i] += (((-1*h[i])/6)*(z[i+1]+ (2*z[i]))) + ((1/(h[i]) * (data[i+1,1]-data[i,1])))
    
    
    
    s_linear = numpy.zeros(x.shape[0])
    
    for i in range(n):
        if i == 0:
            s_linear += (data[i,1] + (x-data[i,0])*(C[i] + ((x - data[0,0])*(B[i] +((x-data[i,0])*A[i]) )))) * (x<=data[i+1,0])
        elif i == n-1:
            s_linear += (data[i,1] + (x-data[i,0])*(C[i] + ((x - data[i,0])*(B[i] +((x-data[i,0])*A[i]) )))) * (x > data[i,0])
        else:
            s_linear += (data[i,1] + (x-data[i,0])*(C[i] + ((x - data[i,0])*(B[i] +((x-data[i,0])*A[i]) )))) * (x<=data[i+1,0]) * (x > data[i,0])
    
    return s_linear



def Funkcja_liniowa(x,data):
    a = 0
    b = 1
    s_linear = numpy.zeros(x.shape[0])
    N = data.shape[0] - 1
    
    while b < 7:
        c = data[b,0]
        d = data[a,0]
        
        if a ==0:
            s_linear += ((((data[b,1]-data[a,1])/(data[b,0]-data[a,0])) * (x-data[a,0])) + data[a,1])*(x<=c)
        elif a==(N-1):
            s_linear += ((((data[b,1]-data[a,1])/(data[b,0]-data[a,0])) * (x-data[a,0])) + data[a,1])*(x>d)
        else:
            s_linear += ((((data[b,1]-data[a,1])/(data[b,0]-data[a,0])) * (x-data[a,0])) + data[a,1])*(x<=c)*(x>d)
          
        b +=1
        a +=1
    return s_linear

def Lagrange(x,data):
    def l(x, data):
        poly = numpy.ones((data.shape[0], x.shape[0]))
        for i in range(data.shape[0]):
            for j in range(data.shape[0]):
                if i != j:
                    poly[i, :] *= (x - data[j, 0]) / (data[i, 0] - data[j, 0])
        return poly
    
    p = numpy.zeros(x.shape[0])
    basis = l(x, data)
    for n in range(data.shape[0]):
        p += basis[n, :] * data[n, 1]
    return p

#Rysowanie wykresu

fig = plt.figure()
axes = fig.add_subplot(1, 1, 1)

axes.plot(x, Funkcja_szescienna(x,data), 'b', label='Funkcja sklejana szescienna') 
axes.plot(x, Funkcja_liniowa(x, data), 'm', label='Funkcja sklejana liniowa')
axes.plot(x, Lagrange(x,data), 'k', label="Wielomian Lagrange'a")
axes.plot(data[:,0], data[:,1], 'ro')
axes.set_title('Interpolacja na różne sposoby')
axes.set_ylim([-5.0, 8.0])
axes.legend(loc=1)
plt.show()
