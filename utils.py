import sympy as smp
import numpy as np
from scipy.integrate import solve_ivp





#idea, le costanti mettile in un array

def eq_motion(L, var, var_d, var_dd, t):
    LE= smp.diff(L, var) - smp.diff(smp.diff(L, var_d), t)
    LE=LE.simplify()
    return(smp.solve(LE, var_dd, simplify = False, rational=False)[0])



# System of first-order ODEs
def system_of_eqs(t, y, m_a, k_a, R_a, m_b, R_b, k_b, eq_theta_a, eq_theta_b):
    theta_a, theta_a_dot, theta_b, theta_b_dot = y
    theta_a_dd = eq_theta_a(m_a, k_a, R_a, theta_a, theta_a_dot, m_b, R_b, k_b, theta_b, theta_b_dot, t)
    theta_b_dd = eq_theta_b(m_a, k_a, R_a, theta_a, theta_a_dot, m_b, R_b, k_b, theta_b, theta_b_dot, t)
    return [theta_a_dot, theta_a_dd, theta_b_dot, theta_b_dd]



def solve_system(initial_conditions, constants, t_0, t_f, t_points, eq_theta_a,eq_theta_b):
    t_span = (t_0, t_f)
    t_eval = np.linspace(t_0, t_f, t_points)
    return(solve_ivp(system_of_eqs, t_span, initial_conditions, args=constants + (eq_theta_a, eq_theta_b), t_eval=t_eval, method='RK45'))

