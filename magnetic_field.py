import numpy as np
from scipy.integrate import quad
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sympy import symbols, integrate, sqrt, cos, sin, simplify, trigsimp, expand, Function



def H_magnetic(a, L, L_1, k, m, t):

    def B_field(x, y, z):
        r=sqrt(x**2+y**2+z**2)
        r_2=r**2
        r_4=r**4
        num_1=(L/2-x)*(-r_4+(-(L/2)**2+x*L+x**2+3*(z**2))*r_2+((L/2)**2)*(x**2+2*z**2)-L*(x**3+2*x*(z**2))-(x*z)**2)
        #print(num_1)
        den_1=(r_2+L/2*(L/2-2*x))**(3/2)
        #print(den_1)
        num_2=(L/2+x)*(r_4-(-L**2/4-x*L+x**2+3*z**2)*r_2-L**2/4*(x**2+2*z**2)-L*(x**3+2*x*z**2)+(x*z)**2)
        #print(num_2)
        den_2=(r_2+L/2*(L/2+2*x))**(3/2)
        #print(den_2)
        den_3=1/(y**2+z**2)**2
        B= a*(num_1/den_1-num_2/den_2)*den_3
        B_simplified=simplify(B)
        B_simplified=trigsimp(B_simplified)
        #print(B_simplified)
        B_simplified = expand(B_simplified)
        #print(B_simplified)
        return ((-m)*B_simplified)



    # Funzione per calcolare l'integrale
    def B_numerical(phi_vals):
        integrals = []  # Lista per salvare i risultati
        for phi in phi_vals:
            # Definizione della funzione integranda per ogni phi
            def integrand(R):
                return B_field(R * np.cos(phi), R * np.sin(phi), -k)
            
            # Calcolo dell'integrale
            result, _ = quad(integrand, -L_1 / 2, L_1 / 2)
            integrals.append(result)  # Aggiungi il risultato alla lista
        
        return integrals  # Ritorna i risultati

    # Vettore di phi
    phi_vals = np.linspace(0, 2 * np.pi, 100)  # 100 punti tra 0 e 2pi

    # Calcolo degli integrali
    integral_results = B_numerical(phi_vals)

    #modello numpy
    def trig_model(phi, A, B, C, n):
        return A * np.cos(n * phi) + B * np.sin(n * phi) + C
    #modello sympy
    def trig_model_sym(phi, A, B, C, n):
        return A * cos(n * phi) + B * sin(n * phi) + C

    # Stima iniziale dei parametri (A, B, C, n)
    initial_guess = [10, 10, -20, 2]  # Valori iniziali stimati manualmente

    # Fitting
    params, params_covariance = curve_fit(trig_model, phi_vals, integral_results, p0=initial_guess)

    # Estrai i parametri ottimizzati
    A, B, C, n = params
    print(f"Parametri ottimizzati: A={A:.3f}, B={B:.3f}, C={C:.3f}, n={n:.3f}")
    the_a = symbols('theta_a', cls=Function)(t)
    the_b = symbols('theta_b', cls=Function)(t)
    phi = the_a-the_b
    print(trig_model_sym(phi, *params))
    plt.figure(figsize=(8, 6))
    plt.plot(phi_vals, integral_results, label=r'$H_B(\Delta\theta)$', color='blue')
    plt.plot(phi_vals, trig_model(phi_vals, *params), label='Fit Trigonometrico', linestyle='--', color='red')
    plt.xlabel(r'$\Delta\theta (rad)$', fontsize=14)
    plt.ylabel(r'Integrale di $H_B (J)$', fontsize=14)
    plt.title('Integrale di $H_B$ al variare di $\Delta\\theta$ ', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True)
    plt.show()
    return(trig_model_sym(phi, *params))


#H_magnetic(a, L, L_1, k, m)

