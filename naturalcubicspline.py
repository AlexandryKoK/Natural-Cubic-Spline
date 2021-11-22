import numpy as np
import matplotlib.pyplot as plt
import scipy
import pandas as pd



#x_points = [ 0.9, 1.3, 1.9, 2.1, 2.6, 3.0, 3.9, 4.4, 4.7, 5.0, 6.0, 7.0, 8.0, 9.2, 10.5, 11.3, 11.6, 12.0, 12.6, 13.0, 13.3]
#y_points = [1.3, 1.5, 1.85, 2.1, 2.6, 2.7, 2.4, 2.15, 2.05, 2.1, 2.25, 2.3, 2.25, 1.95, 1.4, 0.9, 0.7, 0.6, 0.5, 0.4, 0.25]

x_points = [0.24, 0.72, 1.2, 1.6, 2.04, 2.74, 3.16, 3.48, 3.82, 4.14, 4.42, 4.78, 5.16, 5.48, 6.00, 6.48, 7.00, 7.4, 7.74, 8.00, 8.26, 8.76, 9.12, 9.46, 9.78, 10.18, 10.58, 11.00, 11.42, 11.76, 12.44]
y_points = [3.1, 2.78, 2.5, 2.3, 2.14, 2.00, 2.04, 2.38, 2.74, 2.98, 3.26, 3.52, 3.76, 3.92, 4.00, 4.06, 4.00, 4.00, 3.94, 4.00, 3.78, 3.66, 3.6, 3.5, 3.42, 3.2, 3.14, 3.00, 2.84, 2.66, 2.64]

#pontos = pd.DataFrame({'x': x_points, 'f(x)': y_points})
#print(pontos.to_latex(index=False))

#x_points = [0,0.5,1, 1.5, 2]
#y_points = [3, 1.8616, -0.5571, -4.1987, -9.0536]

#x_points = [0.9, 1.3, 1.9, 2.1]
#y_points = [1.3, 1.5, 1.85, 2.1]


n = len(x_points)
#print(n)
def h(x):
    return x_points[x] - x_points[x-1]
#for i in range(1,n):
    #print(h(i))


def alpha(x):
    #print(h(x+1))
    return 6*(y_points[x+1]- y_points[x])/h(x+1)- 6*(y_points[x]-y_points[x-1])/h(x)

#for i in range(n):
    #print(i)
    #print(alpha(i))

a = np.zeros( (n-2, n-2) )
#print(a)
b = np.zeros( (n-2, 1) )
#print(b)
    
for i in range(n-2):
    #print(i)
    for j in range(n-2):
        #print(j)
        if i == j:
            a[i][j] = 2*(h(j+1)+ h(j+2))
        if j == i+1 :
            a[i][j] = h(j+1)
        if j == i-1:
            a[i][j] = h(j+2)

for j in range(0,n-2):
    b[j][0] = alpha(j+1)

#print(a)
#print(b)

g = np.linalg.solve(a, b)
#print(g)

#print(g[1][0])
l = np.zeros( (n-2, 1) )
for i in range(n-2):
    if i == 0:
        #print(i)
        #print(h(i+1))
        l[i][0] = (g[i][0])/(6*h(i+1))
    if i != 0 & i!=n:
        #print(i)
        #print(h(i+1))
        l[i][0] = (g[i][0] - g[i-1][0])/(6*h(i+1))
#print(l)

k = np.zeros( (n-2, 1) )
for i in range(n-2):
    if i == 0:
        k[i][0] = g[i][0]/2
    if i != 0 & i!=n:
        k[i][0] = g[i][0]/2
#print(k)

m = np.zeros( (n-2, 1) )
for i in range(n-2):
    #print(y_points[i+1])
    #print(y_points[i])
    if i == 0:
        m[i][0] = (y_points[i+1]-y_points[i])/(h(i+1)) + (2*h(i+1)*g[i][0])/(6)
    if i != 0 :
        m[i][0] = (y_points[i+1]-y_points[i])/(h(i+1)) + (2*h(i+1)*g[i][0] + g[i-1][0]*h(i+1))/(6)
#print(m)


def splinecubic(x, x_points, l, k, m, y_points,j):
        y = l[j][0]*(x-x_points[j+1])**3+k[j][0]*(x-x_points[j+1])**2+m[j][0]*(x-x_points[j+1]) + y_points[j+1]
        return y

#coeff = pd.DataFrame({'a': np.round(l[0:,0], decimals=2), 'b': np.round(k[0:,0], decimals=2), 'c': np.round(m[0:,0], decimals=2)})
#print(coeff.to_latex(index=False))

for j in range(n-2):
    x = np.linspace(x_points[j],x_points[j+1], 1000)
    #print(x)
    print('%.2f*(x-%.2f)^3 + %.2f(x-%.2f)^2 + %.2f*(x-%.2f) + %.2f {%.2f<=x<=%.2f}' % (l[j][0], x_points[j+1], k[j][0],x_points[j+1], m[j][0], x_points[j+1], y_points[j+1], x_points[j], x_points[j+1]))
    plt.plot(x, splinecubic(x, x_points, l, k, m, y_points,j),lw=3)
plt.savefig('resultado.png')

