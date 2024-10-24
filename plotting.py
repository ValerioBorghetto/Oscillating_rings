
#####################################################################################################################################
#################################### VTK per renderizzare meglio l'immagine, matplotlib sbaglia sulle disanze #######################
#nel video usa plotly

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from functools import partial







######Setting the camera angle###############################

#azim=theta in coord. sferiche///ruota nel piano xy
azim=30
#elev=phi in coord. sferiche///tuota nel piano xz
elev=40
#############################################################



# Funzione per generare l'anello (torus) in forma parametrica
def generate_ring(R, r, u_steps=100, v_steps=100):
    u = np.linspace(0, 2 * np.pi, u_steps)
    v = np.linspace(0, 2 * np.pi, v_steps)
    u, v = np.meshgrid(u, v)
    
    # Equazioni parametriche per l'anello
    x = (R + r * np.cos(v)) * np.cos(u)
    y = r * np.sin(v)
    z = (R + r * np.cos(v)) * np.sin(u)
    
    return x, y, z

# Funzione per applicare rotazioni arbitrarie
def rotate_arbitrary_axis(x, y, z, axis, angle):
    c, s = np.cos(angle), np.sin(angle)
    u, v, w = axis / np.linalg.norm(axis)
    
    rotation_matrix = np.array([
        [c + u*u*(1-c), u*v*(1-c) - w*s, u*w*(1-c) + v*s],
        [v*u*(1-c) + w*s, c + v*v*(1-c), v*w*(1-c) - u*s],
        [w*u*(1-c) - v*s, w*v*(1-c) + u*s, c + w*w*(1-c)]
    ])
    
    xyz = np.vstack([x.flatten(), y.flatten(), z.flatten()])
    rotated_xyz = np.dot(rotation_matrix, xyz)
    
    x_rotated = rotated_xyz[0, :].reshape(x.shape)
    y_rotated = rotated_xyz[1, :].reshape(y.shape)
    z_rotated = rotated_xyz[2, :].reshape(z.shape)
    
    return x_rotated, y_rotated, z_rotated


def update(frame, ax, a_theta_1, a_theta_2, x_1, x_2): #theta 1 e 2 sono array di tutti i valori degli angoli
    # Definire l'asse di rotazione arbitrario
    axis = np.array([0, 0, 1])  # Asse arbitrario per il movimento simile a un giroscopio
    theta_1 = a_theta_1[frame]
    theta_2 = a_theta_2[frame]
    
    # Applicare la rotazione
    x_rotated_1, y_rotated_1, z_rotated_1 = rotate_arbitrary_axis(x_1[0], x_1[1], x_1[2], axis, theta_1)
    x_rotated_2, y_rotated_2, z_rotated_2 = rotate_arbitrary_axis(x_2[0], x_2[1], x_2[2], axis, theta_2)
    
    ax.clear()
    
    # Plottare il nuovo anello
    ax.plot_surface(x_rotated_1, y_rotated_1, z_rotated_1, color='gold', rstride=5, cstride=5, edgecolor='k')
    ax.plot_surface(x_rotated_2, y_rotated_2, z_rotated_2, color='red', rstride=5, cstride=5, edgecolor='k')
    #ax.view_init(elev=60, azim=30)
    return ax,


def update(frame, ax, sol, x_1, x_2, range):
    # Define the arbitrary rotation axis
    #print(f"Frame: {frame}, Type: {type(frame)}")
    axis = np.array([0, 0, 1])  # Arbitrary axis for gyroscope-like motion
    
    # Get the angle theta at the current frame
    theta_1 = sol.y[0][frame]  # theta at the current frame
    theta_2 = sol.y[2][frame]
    
    # Apply the rotation
    x_rotated_1, y_rotated_1, z_rotated_1 = rotate_arbitrary_axis(x_1[0], x_1[1], x_1[2], axis, theta_1)
    x_rotated_2, y_rotated_2, z_rotated_2 = rotate_arbitrary_axis(x_2[0], x_2[1], x_2[2], axis, theta_2)
    
    # Clear the previous plot
    ax.clear()
    
    # Plot the new ring
    ax.plot_surface(x_rotated_1, y_rotated_1, z_rotated_1, color='gold', rstride=5, cstride=5, edgecolor='k')
    ax.plot_surface(x_rotated_2, y_rotated_2, z_rotated_2, color='gold', rstride=5, cstride=5, edgecolor='k')
    # Reset limits and labels for the plot
    ax.set_xlim(range)
    ax.set_ylim(range)
    ax.set_zlim(range)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.view_init(elev=elev, azim=azim)
    
    return ax,



#da aggiungere plotly o vtk 


#R_1>R_2
def simulation(R_1, R_2, sol, t_points): #non c'è la possibilità di settare il centro, è l'origine di default
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x_1, y_1, z_1 = generate_ring(R_1, R_1/10)
    x_a=(x_1, y_1, z_1)
    x_2, y_2, z_2 = generate_ring(R_2, R_2/10)
    x_b=(x_2, y_2, z_2)
    ring_1 = ax.plot_surface(x_1, y_1, z_1, color='gold', rstride=5, cstride=5, edgecolor='k')
    ring_2 = ax.plot_surface(x_2, y_2, z_2, color='red', rstride=5, cstride=5, edgecolor='k')

    # Imposta i limiti e le etichette del grafico
    range = [-2*R_1, 2*R_1]
    ax.set_xlim(range)
    ax.set_ylim(range)
    ax.set_zlim(range)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.view_init(elev=elev, azim=azim)


    ani = FuncAnimation(fig, partial(update, ax=ax, sol=sol, x_1=x_a, x_2=x_b, range=range), frames=np.arange(0, t_points), interval=50)

    # Mostra il grafico
    plt.show()
