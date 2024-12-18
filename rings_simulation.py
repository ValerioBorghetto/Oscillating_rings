import sympy as smp
from utils import *
from plotting_funcs import *
from magnetic_field import *
###################System settings####################

#constants
mu_0=4*np.pi*10**(-7) #H/m

#RING_1
m_val_a = 0.5  # Mass ring a (kg)
R_val_a = 0.30  # Radius ring a (m)
k_val_a =  0.05 # Spring constant ring a (N/rad)
M_val_a = 0.9  # Magnetization ring a (A/m) (For Nickel in a low magnetic field is 3 A/m)
theta_a_0 = 0  # Initial angle (in radians) ring a
theta_a_dot_0 = 10  # Initial angular velocity ring a (rad/s)
 
#RING_2
m_val_b = 0.25  # Mass ring b (kg)
R_val_b = 0.20  # Radius ring b IMPORTANT: (R_a > R_b) (m)
k_val_b = 0.08  # Spring constant ring b (N/rad)
M_val_b = 3  # Magnetization ring b (A/m) (For Nickel in a low magnetic field is 3 A/m)
theta_b_0 = 1  # Initial angle (in radians) ring b
theta_b_dot_0 = 0     # Initial angular velocity ring b (rad/s)

#B 
B_val_est=0 # External magnetic field (T) 

#TIME
t_0 = 0  # Starting time (s)
t_f = 20 # Final time (s)
t_points=1200 #How many time steps (Usually 60*(t_f-t_0))

#ANALYTICAL SOLUTION (Work in progress)
anal_sol=False # If you want an analytical solution

######################################################

assert(R_val_a>R_val_b), "R_a must be bigger than R_b"

#symbols
t, m_a, k_a, R_a, M_a, m_b, R_b, k_b, M_b, B_est= smp.symbols('t m_a k_a R_a M_a m_b R_b k_b M_b B_est', real=True)
the_a = smp.symbols('theta_a', cls=smp.Function)(t)
the_a_d=smp.diff(the_a, t)
the_a_dd=smp.diff(the_a_d, t)
the_b = smp.symbols('theta_b', cls=smp.Function)(t)
the_b_d=smp.diff(the_b, t)
the_b_dd=smp.diff(the_b_d, t)
labels_a=(m_a, k_a, R_a, M_a, the_a, the_a_d)
labels_b=(m_b, R_b, k_b, M_b, the_b, the_b_d)


#constants and initial conditions
constants=(m_val_a, k_val_a, R_val_a, M_val_a, m_val_b, k_val_b, R_val_b, M_val_b, B_est)
initial_conditions=(theta_a_0, theta_a_dot_0, theta_b_0, theta_b_dot_0)

#Angular velocity
omega_a = smp.Matrix([0, 0, the_a_d]) 
omega_b = smp.Matrix([0, 0, the_b_d])

#Inertia tensor
I_a = smp.Matrix([[(1/2) * m_a * R_a**2,0,0],[0,(1/2) * m_a * R_a**2,0],[0,0,m_a * R_a**2]])
I_b = smp.Matrix([[(1/2) * m_b * R_b**2,0,0],[0,(1/2) * m_b * R_b**2,0],[0,0,m_b * R_b**2]])

#Lagrangian
#kinetic therms
T_a=(1/2) * omega_a.T.dot(I_a*omega_a).simplify()
T_b=(1/2) * omega_b.T.dot(I_b*omega_b).simplify()
#Elastic potential
Vel_a=(1/2) * k_a * the_a**2
Vel_b=(1/2) * k_b * (the_a-the_b)**2
#Magnetic moment alignment
#h_a=2*smp.pi*R_a*(R_a/20)*M_a*B_est # using \mu=M*A
#h_b=2*smp.pi*R_b*(R_b/20)*M_b*B_est
#Hm_a=-h_a*smp.cos(the_a)
#Hm_b=-h_b*smp.cos(the_b)
#DA TOGLIERE
a=10#m_val_a*mu_0/(4*np.pi) #m*mu_0/4pi
b=8
L=0.09  
L_1=0.05
k=R_val_a-R_val_b
m=m_val_b
H_mag=H_magnetic(a, L, L_1, k, m, t)
#H_mag_2=H_magnetic(b, L_1, L, k, m_val_a, t)
#Hm_ab1=-M_a*M_b*2*smp.pi*R_a*(R_a/20)*smp.cos(the_a-the_b)
#Hm_ab2 = -M_a*M_b*2*smp.pi*R_b*(R_b/20)*smp.cos(the_a-the_b)
#L= T_a + T_b - V_a - V_b 
H=T_a + T_b + H_mag + Vel_a + Vel_b# + H_mag_2 #+ Vel_a + Vel_b + H_sin# + Hm_a + Hm_b + Hm_ab1 + Hm_ab2
print(H)

#Find the differential equation with Euler-Lagrange
#el_eq_theta_a = EL_eq_motion(L, the_a, the_a_d, the_a_dd, t)
#el_eq_theta_b = EL_eq_motion(L, the_b, the_b_d, the_b_dd, t)

#Find the differential equation with Hamilton
diff_eq_theta_a = H_eq_motion(H, I_a, the_a, the_a_d, the_a_dd, t)
diff_eq_theta_b = H_eq_motion(H, I_b, the_b, the_b_d, the_b_dd, t)

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
dthe_dtdt_a=smp.lambdify( labels_a + labels_b + (t, B_est), diff_eq_theta_a) ##guarda bene il t e il B_est dove va!!!!!
dthe_dtdt_b=smp.lambdify(labels_a + labels_b + (t, B_est), diff_eq_theta_b)
sol = solve_system(constants=constants, t_0=t_0, t_f=t_f, t_points=t_points, initial_conditions=initial_conditions, eq_theta_a=dthe_dtdt_a, eq_theta_b=dthe_dtdt_b)

#Solutions
theta_a_sol = sol.y[0]  # theta_a(t)
theta_a_dot_sol = sol.y[1]  # theta_a_dot(t)
theta_b_sol = sol.y[2]  # theta_b(t)
theta_b_dot_sol = sol.y[3]  # theta_b_dot(t)

plot_results(sol=sol, t_min=t_0, t_max=t_f, t_points=t_points)

simulation(R_1=R_val_a, R_2=R_val_b, sol=sol, t_points=t_points)
