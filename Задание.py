import math
import numpy as np
import matplotlib.pyplot as plt

x_0,x_n = 0,1 					# 0 <= x <= 1					
t_0,t_m = 0,1 					# 0 <= t <= 1
n = 6						# number of x_i
m = 1001						# number of t_i
h = (x_n - x_0)/(n - 1)			# step of x
tau = (t_m - t_0)/(m - 1)		# step of t
x = np.linspace(x_0, x_n, n)	# massiv of x
t = np.linspace(t_0, t_m, m)	# massiv of t
print(h)
print(tau)
a = 0.018
sigma = 0.5 - (h*h)/(12*a*tau)

def u(x,t):
	return 2*x**5 - 3*t**3 + 3*x*t**2 - 2*np.exp(x)

def f(x,t):
	return -9*t**2 + 6*t*x - 40*a*x**3 + 2*a*np.exp(x)

def mu(x):
	return 2*x**5 - 2*np.exp(x)

def mu_1(t):
	return -3*t**3 - 2

def mu_2(t):
	return 3*t**2 - 3*t**3 + 2 - 2 * np.e

def d2f(x,t):
	return -240*a*x + 2*a*np.exp(x)

f_k_j = np.zeros((n,m))
for k in range(1,n-1):
	for j in range(m):
		f_k_j[k][j] = f(x[k],t[j]+tau/2)+(h*h/12)*d2f(x[k],t[j]+tau/2)

u_k_j = np.zeros((n,m))
for k in range(n):
	for j in range(m):
		u_k_j[k][j] = u(x[k],t[j])

u_approx = np.zeros((n,m))
for k in range(n):
	u_approx[k][0] = mu(x[k])
for j in range(m):
	u_approx[0][j] = mu_1(t[j])
	u_approx[n-1][j] = mu_2(t[j])

def Progonka_method(right, C):
    Alphai = np.zeros(n-2)
    Betai = np.zeros(n-2)
    Alphai[1] = (a*tau*sigma)/(h*h + 2*a*tau*sigma)
    Betai[1] = (right[0])/(h*h + 2*a*tau*sigma)
    for i in range(1, n-3):
        Ai = -a*tau*sigma
        Bi = h*h + 2*a*tau*sigma
        Ci = -a*tau*sigma
        Alphai[i + 1] = -Ci/(Ai*Alphai[i]+Bi)
        Betai[i + 1] = (right[i] - Ai*Betai[i])/(Ai*Alphai[i]+Bi)
    An = -a*tau*sigma
    Bn = h*h + 2*a*tau*sigma
    C[n - 3] = (right[n - 3] - An*Betai[n - 3])/(Bn + An*Alphai[n - 3])
    for i in reversed(range(n - 3)):
        C[i] = Alphai[i + 1]*C[i + 1] + Betai[i + 1]
    return C

for j in range(m-1):
	right = np.zeros(n-2)
	right[0] = a*tau * (1-sigma)*  (u_approx[2][j]+u_approx[0][j]) + (h*h - 2*a*tau*(1 - sigma))*u_approx[1][j]+tau*h*h*f_k_j[1][j]+a*tau*sigma*mu_1(t[j+1])
	right[n-3] = a*tau * (1-sigma) * (u_approx[n-3][j]+u_approx[n-1][j]) + (h*h - 2*a*tau *(1 - sigma))*u_approx[n-2][j]+tau*h*h*f_k_j[n-2][j]+a*tau*sigma*mu_2(t[j+1])
	for k in range(1,n-3):
		right[k] = a*tau * (1-sigma) * (u_approx[k][j]+u_approx[k+2][j]) + (h*h - 2*a*tau*(1 - sigma))*u_approx[k+1][j]+tau*h*h*f_k_j[k+1][j]
	cur = np.zeros(n-2)
	cur = Progonka_method(right, u_approx[1:n-1,j+1])

print(np.max(abs(u_k_j - u_approx)))
