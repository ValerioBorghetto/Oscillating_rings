import numpy as np
import sympy as smp
from scipy.integrate import solve_ivp
from utils import *
from plotting_vpython import *
###################System settings####################

#RING_1
m_val_a = 2.5  # Mass ring a (kg)
R_val_a = 0.30  # Radius ring a (m)
k_val_a =  0.4 # Spring constant ring a (N/rad)
theta_a_0 = 0.5  # Initial angle (in radians) ring a
theta_a_dot_0 = 0.0  # Initial angular velocity ring a
 
#RING_2
m_val_b = 0.25  # Mass ring b (kg)
R_val_b = 0.10  # Radius ring b IMPORTANT: (R_a > R_b) (m)
k_val_b = 1  # Spring constant ring b (N/rad)
theta_b_0 = 0.75  # Initial angle (in radians) ring b
theta_b_dot_0 = 0.0  # Initial angular velocity ring b

#TIME
t_0 = 0  # Starting time (s)
t_f = 30 # Final time (s)
t_points=1800 #How many time steps (Usually 60*(t_f-t_0))

#ANALYTICAL SOLUTION
anal_sol=False #If you want an alaytical solution

######################################################

assert(R_val_a>R_val_b), "R_a must be bigger than R_b"

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
omega_a = smp.Matrix([0, 0, the_a_d]) 
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

#Find the differential equation with Euler-Lagrange
diff_eq_theta_a = eq_motion(L, the_a, the_a_d, the_a_dd, t)
diff_eq_theta_b = eq_motion(L, the_b, the_b_d, the_b_dd, t)

"""
print("Equation of motion for theta_a :")
smp.pprint(diff_eq_theta_a)
print("\nEquation of motion for theta_b :")
smp.pprint(diff_eq_theta_b)
"""

#Work in progress, analytical solution
if anal_sol:
    sol_anal(t, the_a, the_b, diff_eq_theta_a, diff_eq_theta_b)

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

plot_theta(theta_b_sol, t_0, t_f, t_points )
plot_theta(theta_a_sol, t_0, t_f, t_points)
