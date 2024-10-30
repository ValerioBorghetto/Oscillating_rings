from vpython import ring, vector, rate, scene, box, cylinder, color
import numpy as np
import matplotlib.pyplot as plt


def update_ring(ring_object, angle):
    #Update the rotation of the ring based on the given angle.
    ring_object.axis = vector(np.cos(angle), np.sin(angle), 0)

def simulation( R_1, R_2, sol, t_points):
    # Create the scene
    theta1_array=sol.y[0]
    theta2_array=sol.y[2]
    
    scene.title = "Oscillalting rings simulation"
    scene.width = 960
    scene.height = 720

    # Set up the camera
    scene.camera.pos = vector(2*R_1, 2*R_1, R_1)  # Position the camera above and slightly behind
    scene.camera.axis = vector(-1, -1, -0.5)  # Point the camera towards the origin
    scene.camera.up = vector(0, 0, 1)  # Set the up direction

    # Create pavement
    pavement_z = -1.5 * R_1
    pavement = box(pos=vector(0, 0, pavement_z), size=vector(5*R_1, 5 *R_1, 0.1), color=vector(0.5, 0.5, 0.5))

    #Create the raws connecting the rings
    cylinder(pos=vector(0, 0, -R_1), axis=vector(0, 0, (R_1-R_2)), color=color.white, radius=R_2/50)
    cylinder(pos=vector(0, 0, pavement_z), axis=vector(0, 0, (-pavement_z-R_1)), color=color.white, radius=R_2/50)

    # Create the rings 
    ring1 = ring(pos=vector(0, 0, 0), radius=R_1, thickness=R_1/20, color=vector(1, 0, 0), axis=vector(1, 0, 0))
    ring2 = ring(pos=vector(0, 0, 0), radius=R_2, thickness=R_2/20, color=vector(0, 0, 1), axis=vector(1, 0, 0))

    # Animation loop
    for i in range(t_points):
        # Update each ring's rotation independently
        if i < len(theta1_array):
            update_ring(ring1, theta1_array[i])
        if i < len(theta2_array):
            update_ring(ring2, theta2_array[i])

        rate(60)  # Controls the animation speed


#plotting both phase space and theta(t) with matplotlib
def plot_results(sol, t_min, t_max, t_points):
    dt=(t_max-t_min)/t_points
    time_array = [t_min + i * dt for i in range(t_points)]
        
    fig, axs = plt.subplots(2, 2, figsize=(10, 8)) 
    #1
    axs[0, 0].plot(sol.y[0], sol.y[2], color='blue')
    axs[0, 0].set_title("Theta_a phase space")
    axs[0, 0].set(xlabel=' Theta_a (rad)', ylabel= ' Momentum_a ')
    #2
    axs[0, 1].plot(sol.y[1], sol.y[3], color='green')
    axs[0, 1].set_title("Theta_b phase space")
    axs[0, 1].set(xlabel=' Theta_b (rad)', ylabel= ' Momentum_b ')
    #3
    axs[1, 0].plot(time_array, sol.y[0], color='red')
    axs[1, 0].set_title("Theta_a over time")
    axs[1, 0].set(xlabel=' Time (s)', ylabel= ' Theta_a (rad)')
    #4
    axs[1, 1].plot(time_array, sol.y[3], color='purple')
    axs[1, 1].set_title("Theta_b over time")
    axs[1, 1].set(xlabel=' Time (s)', ylabel= ' Theta_b (rad)')

    plt.tight_layout()
    plt.show()


#plotting just theta(t)
def plot_theta(name, theta, t_min, t_max, t_points):  
    dt=(t_max-t_min)/t_points
    time_array = [t_min + i * dt for i in range(t_points)]
    while (t < t_max):
        time_array.append(t)
        t = t + dt

    fig, ax = plt.subplots()
    ax.plot(time_array, theta)

    ax.set(xlabel='time (s)', ylabel= name + ' (rad)',) 
    ax.grid()

    plt.show()

#plotting just phasespace
def plot_phasespace(name, theta_array, momentum_array):
    fig, ax = plt.subplots()
    ax.plot(theta_array, momentum_array)

    ax.set(xlabel=name + ' (rad)', ylabel= 'momentum ') 
    ax.grid()

    plt.show()

