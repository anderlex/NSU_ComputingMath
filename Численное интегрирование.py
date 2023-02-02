import math
import numpy as np

a = 0
b = 1
N = 2
znac = 2/math.pi

def func(x):
	return math.sin(math.pi * x)

def left_rectangle(a, b):
	return (b - a) * func(a)

def centre_rectangle(a, b):
	return (b - a) * func((a + b)/2)

def trapezoid(a, b):
	return (b - a) * ((func(a) + func(b))/2)

def Simpson(a, b):
	return ((b - a)/6) * (func(a) + 4 * func((a + b)/2) + func(b))

def Gauss_two(a, b):
	x1 = (1./2)*(a + b) - (math.sqrt(3)/6) * (b - a)
	x2 = (1./2)*(a + b) + (math.sqrt(3)/6) * (b - a)
	return ((b - a)/2) * (func(x1) + func(x2))

def Gauss_three(a, b):
	Xipl = [-math.sqrt(15)/5, 0, math.sqrt(15)/5]
	Wi = [5/9, 8/9, 5/9]
	Ksi = [0] * 3
	for i in range (len(Ksi)):
		Ksi[i] = (b + a)/2 + ((b - a)/2)*Xipl[i]
	sum = 0
	for i in range (3):
		sum += Wi[i]*func(Ksi[i])
	return sum*((b - a)/2)

def Calculate(function, nodes):
	answer = 0
	current = [0] * (len(nodes) - 1)
	for i in range(len(nodes) - 1):
		current[i] = function(nodes[i], nodes[i + 1])
		answer += current[i]
	return answer

def Intervals(N, a, b):
	Xi = [0] * (N + 1)
	for i in range(N + 1):
		Xi[i] = a + ((b - a)/(N))*i
	return Xi

def Upgrade(In, I2n, p):
	alpha = 1/(1 - 2**p)
	U = alpha*In + (1 - alpha)*I2n
	return U


Xi1 = Intervals(N, a, b)
left_answer1 = Calculate(left_rectangle, Xi1)
centre_answer1 = Calculate(centre_rectangle, Xi1)
trapezoid_answer1 = Calculate(trapezoid, Xi1)
Simpson_answer1 = Calculate(Simpson, Xi1)
Gauss_two_answer1 = Calculate(Gauss_two, Xi1)
Gauss_three_answer1 = Calculate(Gauss_three, Xi1)

print(Gauss_three_answer1)

la1 = znac - left_answer1
ca1 = znac - centre_answer1
ta1 = znac - trapezoid_answer1
sa1 = znac - Simpson_answer1
g2a1 = znac - Gauss_two_answer1
g3a1 = znac - Gauss_three_answer1

Xi2 = Intervals(N * 2, a, b)
left_answer2 = Calculate(left_rectangle, Xi2)
centre_answer2 = Calculate(centre_rectangle, Xi2)
trapezoid_answer2 = Calculate(trapezoid, Xi2)
Simpson_answer2 = Calculate(Simpson, Xi2)
Gauss_two_answer2 = Calculate(Gauss_two, Xi2)
Gauss_three_answer2 = Calculate(Gauss_three, Xi2)

la2 = znac - left_answer2
ca2 = znac - centre_answer2
ta2 = znac - trapezoid_answer2
sa2 = znac - Simpson_answer2
g2a2 = znac - Gauss_two_answer2
g3a2 = znac - Gauss_three_answer2

UL = znac - Upgrade(left_answer1, left_answer2, 1)
UC = znac - Upgrade(centre_answer1, centre_answer2, 2)
UT = znac - Upgrade(trapezoid_answer1, trapezoid_answer2, 2)
US = znac - Upgrade(Simpson_answer1, Simpson_answer2, 4)
UG2 = znac - Upgrade(Gauss_two_answer1, Gauss_two_answer2, 4)
UG3 = znac - Upgrade(Gauss_three_answer1, Gauss_three_answer2, 6)


print('                         ' '{:16s} {:16s} {:16s} {:16s}'.format('N', '2N', 'N/2N', 'upg'))
print('Left Rectangle:   ' '{: .14f} {: .14f} {: .14f} {: .14f}'.format(la1, la2, la1/la2, UL))
print('Center Rectangle: ' '{: .14f} {: .14f} {: .14f} {: .14f}'.format(ca1, ca2, ca1/ca2, UC))
print('Trapezoid:        ' '{: .14f} {: .14f} {: .14f} {: .14f}'.format(ta1, ta2, ta1/ta2, UT))
print('Simpson:          ' '{: .14f} {: .14f} {:.14f} {: .14f}'.format(sa1, sa2, sa1/sa2, US))
print('Gauss two:        ' '{: .14f} {: .14f} {:.14f} {: .14f}'.format(g2a1, g2a2, g2a1/g2a2, UG2))
print('Gauss three:      ' '{: .14f} {: .14f} {:.14f} {: .14f}'.format(g3a1, g3a2, g3a1/g3a2, UG3))
