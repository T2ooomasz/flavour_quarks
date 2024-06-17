import python_modules.base as base
import python_modules.layout as layout

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


#
# Potential
#
def V(i, Layout, V_max, r_0):
    dist_sqr = distance_sqr_from_center(i, Layout)
    return V_n_direct(dist_sqr, V_max, r_0)

def V_n_direct(dist_sqr, V_max, r_0):
    dist_sqr = abs(dist_sqr)
    if np.sqrt(dist_sqr) >= r_0:
        return 0
    else:
        return V_max/(r_0*r_0) * dist_sqr - V_max

def distance_sqr_from_center(elem1, Layout):
    site_size = Layout["realLocation"][len(Layout["realLocation"])-1][1][0]
    center = site_size/2
    #for p_elem in Layout["realLocation"]:
    #    if elem1 == p_elem[0]:
    #        x1, y1 = p_elem[1][0], p_elem[1][1]
    x1, y1 = Layout["realLocation"][elem1-1][1][0], Layout["realLocation"][elem1-1][1][1]
    return (x1 - center)**2 + (y1 - center)**2

def plot_potential_on_grid(Layout_base, site_size, site_elements, V_max, r_0):
    a = base.a_calc(site_size, site_elements)
    Layout = layout.initialize_layout(site_elements, a)

    Potencial = np.zeros_like(Layout["2Dneighbours"], dtype=float)
    for i in range(site_elements):
        for j in range(site_elements):
            Potencial[j,i] = V((i*site_elements+j+1), Layout, V_max, r_0)

    X = np.zeros(len(Layout_base["realLocation"]))
    Y = np.zeros(len(Layout_base["realLocation"]))
    for i in range(len(Layout_base["realLocation"])):
        elem = Layout_base["realLocation"][i]
        X[i] = elem[1][0]
        Y[i] = elem[1][1]

    plt.scatter(X, Y, c='blue', linewidths=1)
    plt.gca().set_aspect('equal', 'box')
    plt.imshow(Potencial, cmap=cm.plasma, extent = [0, site_size, 0, site_size])
    plt.colorbar()
    plt.xlabel("x [A]")
    plt.ylabel("y [A]")
    plt.title("Better use plot_potential_on_grid_heat() function")
    plt.show()

def plot_potential_on_grid_heat(Layout_base, site_size, site_elements, V_max, r_0):
    a = base.a_calc(site_size, site_elements)
    Layout = layout.initialize_layout(site_elements, a)
    center = site_size/2
    circle = plt.Circle((center, center), r_0, color='black', fill=False)

    Potencial = np.zeros_like(Layout["2Dneighbours"], dtype=float)
    for i in range(site_elements):
        for j in range(site_elements):
            Potencial[j,i] = V((i*site_elements+j+1), Layout, V_max, r_0)

    X = np.zeros(len(Layout_base["realLocation"]))
    Y = np.zeros(len(Layout_base["realLocation"]))
    for i in range(len(Layout_base["realLocation"])):
        elem = Layout_base["realLocation"][i]
        X[i] = elem[1][0]
        Y[i] = elem[1][1]

    plt.scatter(X, Y, c=Potencial, linewidths=1)
    plt.gca().set_aspect('equal', 'box')
    #plt.plot(center, center, "o", color="black")
    plt.gca().add_patch(circle)
    plt.colorbar()
    plt.xlabel("x [A]")
    plt.ylabel("y [A]")
    title = """Potencial V={V_max}eV on the grid"""
    plt.title(title.format(V_max=V_max))
    plt.show()
