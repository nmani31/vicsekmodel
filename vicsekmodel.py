# Written in Collaboration with Dr. Weijie Chen (Alumnus of Brown University)

import numpy as np
import matplotlib.pyplot as plt

# Fixed values:

R = 37  # confinement radius
v = 40  # bacteria velocity
t = 0  # initial time
T = 3  # final time
dt = 0.1  # time step

# Variables:

N = 600  # number of bacteria
eta = 0.4  # noise coefficient
r = 6  # effective radius


def get_random_particles(radius, N):
    particles = np.zeros((N, 2))
    i = 0
    while i < N:
        # Generate the random point
        x = np.random.uniform(-radius, radius)
        y = np.random.uniform(-radius, radius)
        # Check that it is inside the circle
        if np.sqrt(x ** 2 + y ** 2) < radius:
            particles[i, 0] = x
            particles[i, 1] = y
            i = i + 1
    return particles


def get_neighbors(particles, r, x0, y0):
    neighbors = []

    for j, (x1, y1) in enumerate(particles):
        dist = np.sqrt((x0 - x1) ** 2 + (y0 - y1) ** 2)
        if dist < r:
            neighbors.append(j)

    return neighbors


def get_average(thetas, neighbors):
    n_neighbors = len(neighbors)
    tot_vector = np.zeros(2)
    current_vector = np.zeros(2)

    for index in neighbors:
        theta = thetas[index, 0]
        current_vector[0] = np.cos(theta)
        current_vector[1] = np.sin(theta)
        tot_vector = tot_vector + current_vector

    avg = tot_vector / n_neighbors

    return avg


def plot_vectors(coords, thetas):
    # Generate random color for every particle
    colors = ["b", "g", "y", "m", "c", "pink", "purple", "seagreen",
              "salmon", "orange", "paleturquoise", "midnightblue",
              "crimson", "lavender"]

    for i, (x, y) in enumerate(coords):
        c = colors[i % len(colors)]

        # Plot point
        plt.scatter(x, y, color='b', marker=".")

        # Plot tail
        v = np.zeros(2)
        v[0] = np.cos(thetas[i])
        v[1] = np.sin(thetas[i])
        x1 = x + (2.5 * v[0])
        y1 = y + (2.5 * v[1])
        plt.plot([x, x1], [y, y1], color=c)

    return


def calculate_phi(particles, thetas):
    phi = 0.
    n = len(thetas)
    for i, (x, y) in enumerate(particles):
        mol = np.sqrt(x ** 2 + y ** 2)
        norm_vector = np.array([-y / mol, x / mol])
        vel_vector = np.array([np.cos(thetas[i]), np.sin(thetas[i])])
        phi += np.abs(np.dot(norm_vector, vel_vector))

    return (phi / n - 2 / np.pi) / (1 - 2 / np.pi)


# core program

particles = get_random_particles(R, N)
thetas = np.random.uniform(-np.pi, np.pi, size=(N, 1))

while t < T:

    for i, (x, y) in enumerate(particles):
        rand_angle = np.random.uniform(-np.pi, np.pi)
        ran_x = np.cos(rand_angle)
        ran_y = np.sin(rand_angle)
        ran_vector = np.array([ran_x, ran_y])
        noise = eta * ran_vector

        neighbors = get_neighbors(particles, r, x, y)
        avg = get_average(thetas, neighbors)
        orient_vector = avg + noise

        temp = particles[i] + dt * v * orient_vector

        # Circular boundary condition

        a = temp[0]
        b = temp[1]

        if (a ** 2 + b ** 2) < 1300:
            particles[i] += dt * v * orient_vector
            thetas[i] = np.arctan2(orient_vector[1], orient_vector[0])

        else:
            x1 = particles[i, 0]
            y1 = particles[i, 1]
            beta = np.arctan2(y1, x1)
            thetas[i] = 2 * beta - thetas[i] + np.pi
    t = t + dt

print(u'\u03A6 = %.3f.' % calculate_phi(particles, thetas))
plot_vectors(particles, thetas)
circle2 = plt.Circle((0, 0), radius=R, color='b', fill=False)
ax = plt.gca()
ax.set_xlim((-53.84, 53.84))
ax.set_ylim((-40, 40))
ax.add_artist(circle2)
plt.savefig('both.png')
plt.show()

