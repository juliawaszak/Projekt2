# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 19:22:47 2023

@author: wiecz
"""

import numpy
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt


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
    
    #obliczanie a
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
            s_linear += (data[0,1] + (x-data[0,0])*(C[0] + ((x - data[0,0])*(B[0] +((x-data[0,0])*A[0]) )))) * (x<=data[i+1,0])
        elif i == n-1:
            s_linear += (data[i,1] + (x-data[i,0])*(C[i] + ((x - data[i,0])*(B[i] +((x-data[i,0])*A[i]) )))) * (x > data[i,0])
        else:
            s_linear += (data[i,1] + (x-data[i,0])*(C[i] + ((x - data[i,0])*(B[i] +((x-data[i,0])*A[i]) )))) * (x<=data[i+1,0]) * (x > data[i,0])
    
    return s_linear

cs = CubicSpline(iksy,igreki, bc_type='natural')
y = cs(x)
fig2 = plt.figure()
axes2 = fig2.add_subplot(1,1,1)
axes2.plot(x, Funkcja_szescienna(x,data), 'b', label='Funkcja sklejana szescienna')
axes2.plot(x, y, 'k', label='Interpolacja CubicSpline')
axes2.set_title('Porównanie dwóch metod interpolacji')
axes2.plot(data[:,0], data[:,1], 'ro')
axes2.legend(loc=1)
axes2.set_ylim([-5.0,8.0])
plt.show()

'''
Wykresy te się nie pokrywają. Wynika to z faktu, że wyjsciowym parametrem funkcji CubicSpline jest
bc_type = ‘not-a-knot’. W ustawieniu tym funkcje z pierwszego i drugiego segmentu są tym samym wielomianem,
co powoduje różnicę - w programie napisanym przeze mnie takiego założenia nie ma. Gdyby bc_type zmienić na 
'natural', wykresy by się pokrywały.                                                                                                                                     
'''
