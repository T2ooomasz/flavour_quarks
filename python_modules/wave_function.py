
from python_modules.potential import distance_sqr_from_center

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

#
# Wave function psi
#
def psi_in_dot(psi, i, Layout_base, site_elements, r_0):    # calculate percentage of wave function in a quantum dot (in potential range r_0)
    psi_test = np.reshape(psi[i], (-1, site_elements)).flatten()

    psi_inside = 0
    psi_outside = 0

    for i in range(len(Layout_base["realLocation"])):
        elem = Layout_base["realLocation"][i]
        dist_sqr = distance_sqr_from_center(elem[0], Layout_base)
        if np.sqrt(dist_sqr) < r_0:
            psi_inside += abs(psi_test[i])
        else:
            psi_outside += abs(psi_test[i])

    return 100*psi_inside/(psi_inside+psi_outside)

def plot_psi_percentage_range(psi, n, Layout_base, site_elements, r_0):    # for n first wave fuction plot percentage of wave function in a quantum dot
    Y = np.zeros(n)
    for i in range(n):
        Y[i] = psi_in_dot(psi, i, Layout_base, site_elements, r_0)
    plt.plot(Y)
    plt.xlabel("psi id")
    plt.ylabel("perfent of psi in QD")
    title = """percentage of wave funtion in QD for first {n}"""
    plt.title(title.format(n=n))
    plt.show()

def plot_psi_on_grid(psi, i, Layout_base, site_size, site_elements, r_0):   # plot i-th wave function on the grid 
    center = site_size/2
    psi_test = np.reshape(psi[i], (-1, site_elements))
    circle = plt.Circle((center, center), r_0, color='black', fill=False)
    percentage = psi_in_dot(psi, i, Layout_base, site_elements, r_0)

    X = np.zeros(len(Layout_base["realLocation"]))
    Y = np.zeros(len(Layout_base["realLocation"]))
    for k in range(len(Layout_base["realLocation"])):
        elem = Layout_base["realLocation"][k]
        X[k] = elem[1][0]
        Y[k] = elem[1][1]

    plt.scatter(X, Y, c=psi_test.flatten(), linewidths=1)
    plt.gca().set_aspect('equal', 'box')
    plt.gca().add_patch(circle)
    plt.colorbar()
    plt.xlabel("x [A]")
    plt.ylabel("y [A]")
    title = """psi {i} # {percentage:.5f}% in QD"""
    plt.title(title.format(i=i, percentage=percentage))
    plt.show()

def psi_in_dot_1D(psi, i, site_size, site_elements, r_0):
    center = site_size/2
    x = np.linspace(0, site_size, site_elements)
    psi_test = np.reshape(psi[i], (-1, site_elements))[int(site_elements/2)]
    plt.plot(x, psi_test)
    plt.xlabel("x [A]")
    plt.ylabel("y [A]")
    plt.vlines(center-r_0, np.amin(psi_test), np.amax(psi_test), colors="black", linestyles="--")
    plt.vlines(center+r_0, np.amin(psi_test), np.amax(psi_test), colors="black", linestyles="--")
    title = """psi {i} profile through the center"""
    plt.title(title.format(i=i))
    plt.show()

def psi_n_3D(psi, n, site_size, site_elements):
    x = np.linspace(0,site_size,site_elements)
    xx, yy = np.meshgrid(x, x)

    #R = ((a*xx)**2 + (a*yy)**2)

    #f = lambda R, V_max, r_0: V_n(R, V_max, r_0)
    #Potencial = np.array(list(map(f, R, np.full(len(R), V_max), np.full(len(R), r_0))))
    psi_test = np.reshape(psi[n], (-1, site_elements))
    # show all Potencial
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    surf = ax.plot_surface(xx, yy, psi_test, cmap=cm.coolwarm,
                        linewidth=0, antialiased=False)
    title = """psi {i}"""
    plt.title(title.format(i=n))
    plt.show()
    