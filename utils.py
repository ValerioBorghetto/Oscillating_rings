import sympy as smp
import numpy as np
from scipy.integrate import solve_ivp

def eq_motion(L, var, var_d, var_dd, t):
    LE= smp.diff(L, var) - smp.diff(smp.diff(L, var_d), t)
    LE=LE.simplify()
    return(smp.solve(LE, var_dd, simplify = False, rational=False)[0])


def system_of_eqs(t, y, m_a, k_a, R_a, m_b, R_b, k_b, eq_theta_a, eq_theta_b):
    theta_a, theta_a_dot, theta_b, theta_b_dot = y
    theta_a_dd = eq_theta_a(m_a, k_a, R_a, theta_a, theta_a_dot, m_b, R_b, k_b, theta_b, theta_b_dot, t)
    theta_b_dd = eq_theta_b(m_a, k_a, R_a, theta_a, theta_a_dot, m_b, R_b, k_b, theta_b, theta_b_dot, t)
    return [theta_a_dot, theta_a_dd, theta_b_dot, theta_b_dd]


def solve_system(initial_conditions, constants, t_0, t_f, t_points, eq_theta_a,eq_theta_b):
    t_span = (t_0, t_f)
    t_eval = np.linspace(t_0, t_f, t_points)
    return(solve_ivp(system_of_eqs, t_span, initial_conditions, args=constants + (eq_theta_a, eq_theta_b), t_eval=t_eval, method='RK45'))

#Work in progress
def sol_anal(t, the_a, the_b, diff_eq_a, diff_eq_b):
    #for the initial conditions
    theta_a0, omega_a0, theta_b0, omega_b0 = smp.symbols('theta_a0 omega_a0 theta_b0 omega_b0')   
    diff_eqs = [
    smp.Eq(the_a.diff(t, 2), diff_eq_a),
    smp.Eq(the_b.diff(t, 2), diff_eq_b)
    ]
    #with initial conditions not able to solve
    solutions = smp.dsolve(diff_eqs)#, ics={the_a.subs(t, 0): theta_a0, the_a.diff(t).subs(t, 0): omega_a0, the_b.subs(t,0):theta_b0, the_b.diff(t).subs(t,0): omega_b0})
    print("\\Theta_a solution:")
    smp.pprint(solutions[0].simplify())
    print("")
    print("\\Theta_b solution:")
    smp.pprint(solutions[1].simplify()) 

    