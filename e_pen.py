import numpy as np 
from numpy.linalg import inv
from matplotlib import pyplot as plt
from math import cos, sin, tan, pi


def G(y,t): 
	x_d, θ_d, x, θ = y[0], y[1], y[2], y[3]

	x_dd = (l0+x) * θ_d**2 - k/m*x + g*cos(θ)
	θ_dd = -2.0/(l0+x) * x_d * θ_d - g/(l0+x) * sin(θ)	

	return np.array([x_dd, θ_dd, x_d, θ_d])


def RK4_step(y, t, dt):
	k1 = G(y,t)
	k2 = G(y+0.5*k1*dt, t+0.5*dt)
	k3 = G(y+0.5*k2*dt, t+0.5*dt)
	k4 = G(y+k3*dt, t+dt)

	return dt * (k1 + 2*k2 + 2*k3 + k4) /6

# variables
m = 2.0
l0 = 1.0
g = 9.81
k = 100.0

delta_t = 0.01
time = np.arange(0.0, 5.0, delta_t)

# initial state
y = np.array([0, 0.0, 0.0 , 1.0])   # [velocity, displacement]

Y1 = []
Y2 = []

# time-stepping solution
for t in time:
	y = y + RK4_step(y, t, delta_t) 

	Y1.append(y[2])
	Y2.append(y[3])


# plot the result
plt.plot(time,Y1)
plt.plot(time,Y2)
plt.grid(True)
plt.legend(['x', 'θ'], loc='lower right')
plt.show()