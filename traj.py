


import numpy as np
import math as mt
import datetime as dt
import matplotlib
#matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.animation as animation

import parameters as par

print("è un cazzo di script\n")


#xd = fun.DirectKinematics(vertcat(mt.radians(70.), mt.radians(-70.)))

def DirectKinematics(q):	
    x = par.a1*mt.cos(q[0]) + par.a2*mt.cos(q[0]+q[1])
    y = par.a1*mt.sin(q[0]) + par.a2*mt.sin(q[0]+q[1])
	#x = par.a1*mt.cos(q0) + par.a2*mt.cos(q0+q1))
	#y = par.a1*mt.sin(q0) + par.a2*mt.sin(q0+q1)
    p = np.transpose(np.array([x,y]))
    return p

def InverseKinematics(x):
	cosq2 = (x[0]**2 + x[1]**2 - par.a1**2 - par.a2**2) / (2 * par.a1 * par.a2)
	dif = (1.0 - (cosq2**2))
	try:
		sinq2 = -mt.sqrt(dif)
	except ValueError:
		sinq2 = 0
	q2 = mt.atan2(sinq2, cosq2)
	q1 = mt.atan2(x[1], x[0]) - mt.atan2((par.a2*sinq2), (par.a1 + par.a2*cosq2))
	q1 = round(q1,3)
	q2 = round(q2,3)
	result = np.array([q1,q2])
	return result

def polyInterp(xi, dxi, xf, dxf):
	tf = par.N+1	
	p0 = xi
	p1 = dxi
	p2 = 3.0/(tf**2) * (xf - xi) - 1.0/tf * (2.0*dxi + dxf)
	p3 = 2.0/(tf**3) * (xi - xf) + 1.0/(tf**2) * (dxi + dxf)
	return np.array([p0, p1, p2, p3])


q0_vector = np.radians(np.array([0,-40]))
dq0_vector = np.radians(np.array([0,0]))
qf_vector = np.radians(np.array([70,-70]))
dqf_vector = np.radians(np.array([0,0]))
print("q0: {} qf:{}".format(q0_vector,qf_vector))
p0 = DirectKinematics(q0_vector)
pf = DirectKinematics(qf_vector)
print("p0: {} pf:{}".format(p0,pf))
q0_vector_IK = InverseKinematics(p0)
qf_vector_IK = InverseKinematics(pf)
print("q0_IK: {} qf_IK:{}".format(q0_vector_IK,qf_vector_IK))

print("Calcolo Traiettorie")

pol0 = polyInterp(q0_vector[0], dq0_vector[0], qf_vector[0],  dqf_vector[0])
pol1 = polyInterp(q0_vector[1], dq0_vector[1], qf_vector[1],  dqf_vector[1])

q0_traj = q0_vector_IK[0]
q1_traj = q0_vector_IK[1]
dq0_traj = 0
dq1_traj = 0
x_traj = p0[0]
y_traj = p0[1]
i = 0

while i < par.N:
	q0_now = pol0[3]*i**3+pol0[2]*i**2+pol0[1]*i + pol0[0]
	q1_now = pol1[3]*i**3+pol1[2]*i**2+pol1[1]*i + pol1[0]
	q0_traj = np.append(q0_traj,q0_now)
	q1_traj = np.append(q1_traj,q1_now)

	p_now = DirectKinematics([q0_now,q1_now])
	x_traj = np.append(x_traj,p_now[0])
	y_traj = np.append(y_traj,p_now[1])

	dq0_traj = np.append(dq0_traj,3*pol0[3]*i**2+2*pol0[2]*i+pol0[1])
	dq1_traj = np.append(dq1_traj,3*pol1[3]*i**2+2*pol1[2]*i+pol1[1])
	i = i+1

plt.ion()
fig = plt.figure(figsize=(10, 10))
#ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(-3, 3), ylim=(-3, 3))
#ax.grid()
plt.title("Trajectory")
plt.plot( q0_traj, 'r--', label='q0_sol')
plt.plot( q1_traj, 'b--', label='q1_sol')
plt.legend()
plt.savefig('q1q2.png')
fig = plt.figure(figsize=(10, 10))
plt.title("d_Trajectory")
plt.plot( dq0_traj, 'r--', label='dq0_sol')
plt.plot( dq1_traj, 'b--', label='dq1_sol')
plt.legend()
fig = plt.figure(figsize=(10, 10))
plt.title("Trajectory nello spazio operativo")
plt.plot( x_traj, y_traj)
plt.savefig('traj_spazio_op_da_q.png')
plt.show()



print("Calcolo Traiettorie su spazio operativo p0 e pf")

pol0_so = polyInterp(p0[0], 0, pf[0],  0)
pol1_so = polyInterp(p0[1], 0, pf[1],  0)

x_traj = p0[0]
y_traj = p0[1]
dx_traj = 0
dy_traj = 0

i = 0

while i < par.N:
    x_traj = np.append(x_traj,pol0_so[3]*i**3+pol0_so[2]*i**2+pol0_so[1]*i + pol0[0])
    y_traj = np.append(y_traj,pol1_so[3]*i**3+pol1_so[2]*i**2+pol1_so[1]*i + pol1_so[0])
    dx_traj = np.append(dx_traj,3*pol0[3]*i**2+2*pol0[2]*i+pol0[1])
    dy_traj = np.append(dy_traj,3*pol1[3]*i**2+2*pol1[2]*i+pol1[1])
    i = i+1

fig = plt.figure(figsize=(10, 10))
#ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(-3, 3), ylim=(-3, 3))
#ax.grid()
plt.title("Trajectory on Operative Space")
#plt.plot( x_traj,y_traj , 'r--')
plt.plot(x_traj,y_traj , 'r--')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.legend()
plt.show()

print("Calcolo Traiettorie su spazio dei giunti invertendo la cinematica")

pol0_so = polyInterp(p0[0], 0, pf[0],  0)
pol1_so = polyInterp(p0[1], 0, pf[1],  0)

x_traj = p0[0]
y_traj = p0[1]
dx_traj = 0
dy_traj = 0

i = 1

while i < par.N:
    x_traj = np.append(x_traj,pol0_so[3]*i**3+pol0_so[2]*i**2+pol0_so[1]*i + pol0[0])
    y_traj = np.append(y_traj,pol1_so[3]*i**3+pol1_so[2]*i**2+pol1_so[1]*i + pol1_so[0])
    dx_traj = np.append(dx_traj,3*pol0[3]*i**2+2*pol0[2]*i+pol0[1])
    dy_traj = np.append(dy_traj,3*pol1[3]*i**2+2*pol1[2]*i+pol1[1])
    i = i+1

i = 1

q_vector_IK = InverseKinematics([x_traj[0],y_traj[0]])
q0_traj = q_vector_IK[0]
q1_traj = q_vector_IK[1]

while i < par.N:
	p_i = np.array([x_traj[i],y_traj[i]])
	q_vector_IK = InverseKinematics(p_i)
	q0_traj = np.append(q0_traj,q_vector_IK[0])
	q1_traj = np.append(q1_traj,q_vector_IK[1])
	i = i + 1

fig = plt.figure(figsize=(10, 10))
#ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(-3, 3), ylim=(-3, 3))
#ax.grid()
plt.title("Trajectory on inverted kinematics con le traiettorie calcolate negli spazi dei giunti")
plt.plot( q0_traj, 'r--', label='q0_sol')
plt.plot( q1_traj, 'b--', label='q1_sol')
plt.legend()
plt.savefig('traj_spazio_op.png')
plt.show(block=True)