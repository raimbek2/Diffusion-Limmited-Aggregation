import numpy as np
import random as rnd
from matplotlib import pyplot as plt
from time import time

t0 = time()
# Constants
N = 501  # grid side
s = (N, N)
grid = np.zeros(s)
P_nn = 1
R_max = 20

# Seed
Seed_x = round(N / 2)
Seed_y = round(N / 2)
grid[Seed_x][Seed_y] = 1
# Cluster
Cluster = []
# print('Seed position:', round(N / 2), round(N / 2))


# Clear the grid
def grid_clear():
    for i in range(N):
        for j in range(N):
            grid[i][j] = 0
    grid[Seed_x][Seed_y] = 1


# Particles class
class Particle:
    x = 0
    y = 0
    bond = False

    def __init__(self, distance):
        Rad = R_max + 5 - distance
        # Random position for the particle on the circle of R_max + 5 radius
        x_cor = rnd.randint(-Rad, Rad)
        if rnd.random() < 0.5:
            y_cor = - round(np.sqrt(Rad ** 2 - x_cor ** 2))
        else:
            y_cor = + round(np.sqrt(Rad ** 2 - x_cor ** 2))
        x_cor += Seed_x
        y_cor += Seed_y
        self.x = x_cor
        self.y = y_cor
        grid[x_cor][y_cor] = 1

    def __del__(self):
        grid[self.x][self.y] = 0

    # Perform a random jump to the next position
    def jump(self):
        grid[self.x][self.y] = 0
        temp = rnd.random()
        x_disp = 0
        y_disp = 0
        if temp < 0.25:
            x_disp += 1
        elif temp < 0.5:
            x_disp -= 1
        elif temp < 0.75:
            y_disp += 1
        else:
            y_disp -= 1
        # Make sure that a particles don't get into each other.
        if grid[self.x + x_disp][self.y + y_disp] != 1:
            self.x += x_disp
            self.y += y_disp
        grid[self.x][self.y] = 1

    # Check the particle if it has a cluster in the neighborhood
    # or if it went out of range
    def check_in(self):
        if grid[self.x + 1][self.y] or grid[self.x - 1][self.y] \
                or grid[self.x][self.y - 1] or grid[self.x][self.y + 1]:
            if rnd.random() < P_nn:
                self.bond = True
                Cluster.append(self)
                return 1
        return 0

    def check_out(self):
        # More than 3 R_max
        if (self.x - Seed_x) ** 2 + (self.y - Seed_y) ** 2 > 9 * R_max ** 2:
            return 1
        return 0


# Main process
def simulate(n, pn):
    global R_max
    global P_nn
    P_nn = pn
    for i in range(n):
        # Offset from R_max + 5
        offset = 0
        # Create a prototype particle at the R_max + 5
        prototype1 = Particle(offset)
        # Perform jumps until the particle gets to the cluster
        while prototype1.check_in() != 1:
            # print('Postion: ',prototype1.x,prototype1.y)
            # If a particle went out of bounds, delete, create new
            if prototype1.check_out() == 1:
                # Delete old particle
                # print('Deleting')
                # Make it closer to the seed
                if offset < R_max - 5:
                    offset += 1
                del prototype1
                # Create new on the circle
                prototype1 = Particle(offset)
            # Jump
            prototype1.jump()
        # print(len(Cluster),"Captured particle at:", prototype1.x, prototype1.y)
        #         print('R_max = ',R_max)
        # In order to optimize the computational time, I decided that we can
        if (i + 20) % 100 == 0:
            R_max = 2 * round(np.sqrt(i))


# Create a 10000 points simulation
num_points = 500  # you can chage this and make sure my program works for lower number of particles

tot_cells_at = np.zeros(51)
# Number of points at a distance x
for i in range(1, 50):
    tot_cells_at[i] = 4 * i


# Distance between two particles
def distance(par1, par2):
    return abs(par2.x - par1.x) + abs(par2.y - par1.y)


tt = time()


# BFS algorithm
def compute(cluster):
    # Density-density correlation
    C_all = np.zeros(51)
    for idx in range(num_points):
        X = cluster[idx].x
        Y = cluster[idx].y
        processed = {(X, Y)}
        queue = [(X, Y)]
        for level in range(1, 50):
            new_queue = []
            for point in queue:
                x_coor = point[0]
                y_coor = point[1]
                if not (x_coor - 1, y_coor) in processed:
                    processed.add((x_coor - 1, y_coor))
                    new_queue.append((x_coor - 1, y_coor))
                if not (x_coor, y_coor - 1) in processed:
                    processed.add((x_coor, y_coor - 1))
                    new_queue.append((x_coor, y_coor - 1))
                if not (x_coor + 1, y_coor) in processed:
                    processed.add((x_coor + 1, y_coor))
                    new_queue.append((x_coor + 1, y_coor))
                if not (x_coor, y_coor + 1) in processed:
                    processed.add((x_coor, y_coor + 1))
                    new_queue.append((x_coor, y_coor + 1))
            queue = new_queue

        for point in cluster:
            if (point.x, point.y) in processed:
                dist = distance(point, cluster[idx])
                C_all[dist] += 1
    # Normalization
    for i in range(1, 50):
        C_all[i] /= tot_cells_at[i] * num_points
    temp = np.zeros(49)
    for i in range(1, 50):
        temp[i-1] = C_all[i]
    return temp


# calculate the avarage for 20 clusters
all_clusters = []

# For P_nn = 1
for i in range(20):
    c = []
    print('P_n=1 ', i ,'# run')
    Cluster.clear()
    grid_clear()
    R_max = 20
    simulate(num_points, 1)
    for j in Cluster:
        c.append(j)
    all_clusters.append(c)
all_c = []
for i in range(20):
    all_c.append(compute(all_clusters[i]))
result_c = np.zeros(49)
for i in range(49):
    result_c[i] = np.mean(np.transpose(all_c)[i])
lnr = np.zeros(49)
lnc = np.zeros(49)
for var in range(1, 50):
    lnr[var-1] = np.log(var)
for var in range(49):
    lnc[var] = np.log(result_c[var])


#Plot the c for P_nn =1
plt.figure()
plt.plot(lnr, lnc)
plt.title('Logarithms relation of Density and radius at $P_{nn} = 1$')
plt.xlabel('$\ln{r}$')
plt.ylabel('$\ln(C(r))$')
plt.show()

# For P_nn = 0.3
all_c.clear()
all_clusters.clear()


for i in range(5):
    c = []
    print('P_n=1 ', i ,'# run')
    Cluster.clear()
    grid_clear()
    R_max = 20
    simulate(num_points, 1)
    for j in Cluster:
        c.append(j)
    all_clusters.append(c)
for i in range(5):
    all_c.append(compute(all_clusters[i]))
result_c = np.zeros(49)
for i in range(49):
    result_c[i] = np.mean(np.transpose(all_c)[i])
lnc = np.zeros(49)
for var in range(49):
    lnc[var] = np.log(result_c[var])

# Plot the c for P_nn =0.3
plt.figure()
plt.plot(lnr, lnc)
plt.title('Logarithms relation of Density and radius at $P_{nn} = 0.3$')
plt.xlabel('$\ln{r}$')
plt.ylabel('$\ln(C(r))$')
plt.show()