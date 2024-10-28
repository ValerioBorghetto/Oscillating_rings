from vpython import ring, vector, rate, scene, box, cylinder, color
import numpy as np
import matplotlib.pyplot as plt


def update_ring(ring_object, angle):
    #Update the rotation of the ring based on the given angle.
    ring_object.axis = vector(np.cos(angle), np.sin(angle), 0)  # Update the axis direction

def simulation( R_1, R_2, sol, t_points):
    # Create the scene
    theta1_array=sol.y[0]
    theta2_array=sol.y[2]
    scene.title = "Oscilalting rings simulation"
    scene.width = 800
    scene.height = 600

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
    ring1 = ring(pos=vector(0, 0, 0), radius=R_1, thickness=R_1/10, color=vector(1, 0, 0), axis=vector(1, 0, 0))
    ring2 = ring(pos=vector(0, 0, 0), radius=R_2, thickness=R_2/10, color=vector(0, 0, 1), axis=vector(1, 0, 0))

    # Animation loop
    for i in range(t_points):
        # Update each ring's rotation independently
        if i < len(theta1_array):
            update_ring(ring1, theta1_array[i])
        if i < len(theta2_array):
            update_ring(ring2, theta2_array[i])

        rate(60)  # Controls the animation speed

#plotting theta(t)
def plot_theta(theta, t_min, t_max, t_points):  
    dt=(t_max-t_min)/t_points
    t=t_min + dt
    time_array=[]
    while (t < t_max):
        time_array.append(t)
        t = t + dt

    fig, ax = plt.subplots()
    ax.plot(time_array, theta)

    ax.set(xlabel='time (s)', ylabel='theta (rad)',)
    ax.grid()

    plt.show()








