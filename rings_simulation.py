import numpy as np
import sympy as smp
from scipy.integrate import solve_ivp
from utils import *
from plotting import *


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation


###############Settings########################

####RING_1####
m_val_a = 1  # Mass ring a
R_val_a = 5.0  # Radius ring a
k_val_a = 1.0  # Spring constant ring a
theta_a_0 = 1.5  # Initial angle (in radians) ring a
theta_a_dot_0 = 0.0  # Initial angular velocity ring a

####RING_2####
m_val_b = 0.5  # Mass ring b
R_val_b = 2.0  # Radius ring b IMPORTANT: (R_a > R_b)
k_val_b = 1  # Spring constant ring b
theta_b_0 = 0.75  # Initial angle (in radians) ring b
theta_b_dot_0 = 0.0  # Initial angular velocity ring b

####Time####
t_0 = 0  # Starting time
t_f = 10 # Final time
t_points=100 #How many time steps

####Analytical solution####
anal_sol=True
###############################################

#symbols
t, m_a, k_a, R_a, m_b, R_b, k_b = smp.symbols('t m_a k_a R_a m_b R_b k_b', real=True)
the_a = smp.symbols('theta_a', cls=smp.Function)(t)
the_a_d=smp.diff(the_a, t)
the_a_dd=smp.diff(the_a_d, t)
the_b = smp.symbols('theta_b', cls=smp.Function)(t)
the_b_d=smp.diff(the_b, t)
the_b_dd=smp.diff(the_b_d, t)
labels_a=(m_a, k_a, R_a, the_a, the_a_d)
labels_b=(m_b, R_b, k_b, the_b, the_b_d)


#constants and initial conditions
constants=(m_val_a, k_val_a, R_val_a, m_val_b, k_val_b, R_val_b)
initial_conditions=(theta_a_0, theta_a_dot_0,theta_b_0, theta_b_dot_0)

#Angular velocity
omega_a = smp.Matrix([0, 0, the_a_d]) #rendila generale
omega_b = smp.Matrix([0,0, the_b_d])

#Inertia tensor
I_a = smp.Matrix([[(1/2) * m_a * R_a**2,0,0],[0,(1/2) * m_a * R_a**2,0],[0,0,m_a * R_a**2]])
I_b = smp.Matrix([[(1/2) * m_b * R_b**2,0,0],[0,(1/2) * m_b * R_b**2,0],[0,0,m_b * R_b**2]])

#Lagrangian
T_a=(1/2) * omega_a.T.dot(I_a*omega_a).simplify()
T_b=(1/2) * omega_b.T.dot(I_b*omega_b).simplify()
V_a=(1/2) * k_a * the_a**2
V_b=(1/2) * k_b * (the_b-the_a)**2
L= T_a + T_b - V_a - V_b

#Solving equation of motion
diff_eq_theta_a = eq_motion(L, the_a, the_a_d, the_a_dd, t)

diff_eq_theta_b = eq_motion(L, the_b, the_b_d, the_b_dd, t)

print(the_a_d)

print("Equation of motion for theta_a:")
smp.pprint(diff_eq_theta_a)
print("\nEquation of motion for theta_b:")
smp.pprint(diff_eq_theta_b)

"""
#Work in progress
if anal_sol:
    eq_a = smp.Eq(the_a.diff(t, 2), diff_eq_theta_a)
    eq_b = smp.Eq(the_b.diff(t, 2), diff_eq_theta_b)
    theta_a0, omega_a0, theta_b0, omega_b0 = smp.symbols('theta_a0 omega_a0 theta_b0 omega_b0')

    # Solve the equations
    anal_sol_a = smp.dsolve(eq_a, the_a, ics={the_a.subs(t, 0): theta_a0, the_a.diff(t).subs(t, 0): omega_a0})
    anal_sol_b = smp.dsolve(eq_b, the_b,ics={the_b.subs(t, 0): theta_b0, the_b.diff(t).subs(t, 0): omega_b0})

    anal_sol_b = anal_sol_b.rhs

    integral_expr = smp.integrate(anal_sol_b, t)  # Integrate theta_b with respect to t
    anal_sol_a = anal_sol_a.rhs.subs(the_b, integral_expr)


    print("Analytical solution for theta_a:")
    smp.pprint(anal_sol_a)
    print("\nAnalytical solution for theta_b:")
    smp.pprint(anal_sol_b)

"""
#Lambdify the equations
dthe_dtdt_a=smp.lambdify(labels_a + labels_b + (t,), diff_eq_theta_a)
dthe_dtdt_b=smp.lambdify(labels_a + labels_b + (t,), diff_eq_theta_b)
sol = solve_system(constants=constants, t_0=t_0, t_f=t_f, t_points=t_points, initial_conditions=initial_conditions, eq_theta_a=dthe_dtdt_a, eq_theta_b=dthe_dtdt_b)

#Solutions
theta_a_sol = sol.y[0]  # theta_a(t)
theta_a_dot_sol = sol.y[1]  # theta_a_dot(t)
theta_b_sol = sol.y[2]  # theta_b(t)
theta_b_dot_sol = sol.y[3]  # theta_b_dot(t)

simulation(R_1=R_val_a, R_2=R_val_b, sol=sol, t_points=t_points)


