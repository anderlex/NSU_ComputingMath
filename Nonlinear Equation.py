import math

a = 1
b = 100000
eps = 1.0e-05
N = 10000000

def f(x):
    return -0.01 + 1/x

def df(x):
    return -1/(x*x)

def Chord_method(a, b, eps):
    i = 0
    c = a - f(a)*(b - a)/(f(b) - f(a))
    if f(a)-2*f((a+b)/2)+f(b)> 0: 
        x = b
    else: x = a
    if f(a)*f(b)<0:
        while abs(f(c)) > eps * abs(df(x)) and i < N:
            i = i + 1
            c = a - f(a)*(b - a)/(f(b) - f(a))
            if f(a)*f(c)<0:
                b = c
            if f(c)*f(b)<0:
                a = c
    print("Approximate solution =", c)
    print("Number of iterations =", i)
    print("Achieved accuracy =", abs(f(c)))

def Secant_method(a, b, eps):
    i = 0
    if f(a)-2*f((a+b)/2)+f(b) > 0: 
        x0 = a
        x = b
    else: 
        x0 = b
        x = a
    x1 = x0 - f(x0)/df(x0)
    print('Begin = ',a,'End = ', b ,',\n First appr = ',x0,',Second appr = ',x1)
    while abs(f(x1)) > eps*abs(df(x)) and i < N:
        i = i + 1
        x0, x1 = x1, x1 - f(x1)*(x1-x0)/(f(x1)-f(x0))
    print("Approximate solution =", x1)
    print("Number of iterations =", i)
    print("Achieved accuracy =", abs(f(x1)))

#main
print("Chord method:")
Chord_method(a, b, eps)
print("Secant method:")
Secant_method(a, b, eps)
