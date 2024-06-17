# required libraries
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import math
import scipy
import itertools
import unittest
import pandas as pd
import ast


#
# Fundamental values of the grid
#
def t_calc(
        a,
        m_eff,
        m_0,
        hbar,
        _1eV
):
    return 1 / (2 * m_eff * a ** 2) * (hbar ** 2) / (m_0 * 10e-20 * _1eV)   # [eV]

def a_calc(
        site_size,
        site_elements
):
    return site_size / (site_elements - 1) # [Angstrom] - distance between elements 

#
# Dispersion 2D
#
def dispersion_2D(
    site_elements,
    t,
    a,
    hbar,
    m_eff, 
    m_0,
    _1eV,
    E0,
    dEs):
    pi = np.pi
    kx = np.linspace(-pi/a,pi/a,site_elements)
    ky = np.linspace(-pi/a,pi/a,site_elements)
    kxx, kyy = np.meshgrid(kx, ky)
    dEs = 0

    # kp
    f = lambda kx, ky: (hbar ** 2) * (kx ** 2 + ky ** 2) / (2 * m_eff * m_0 * _1eV * 10e-20) + E0
    ek_kp = np.array(list(map(f, kxx, kyy)))


    # cos
    ek_cos = dEs - 2 * t * (np.cos(kxx * a) + np.cos(kyy * a))

    # scale the ek_cos so it will start at 0
    f = lambda x: x - delta
    delta = np.amin(ek_cos)
    ek_cos = np.array(list(map(f, ek_cos)))

    # plot dispersion
    plt.plot(kx,ek_cos[0],'o', mfc = 'none', mec = 'blue', color="blue", label="cos()")
    for i in range(1,site_elements):
        plt.plot(kx,ek_cos[i],'o', mfc = 'none', mec = 'blue', color="blue")
    plt.plot(kx,ek_kp[0],'.', color="red", label="kp")
    for i in range(1,site_elements):    
        plt.plot(kx,ek_kp[i],'.', color="red")
    plt.xlabel("k (1/A)")
    plt.ylabel("E (eV)")
    plt.hlines(8*t, -pi/a, pi/a, colors="black", linestyles="--")
    plt.hlines(0, -pi/a, pi/a, colors="black", linestyles="--")
    plt.vlines(-pi/a, 0, np.amax(ek_kp), colors="black", linestyles="--")
    plt.vlines(pi/a, 0, np.amax(ek_kp), colors="black", linestyles="--")
    plt.legend()
    plt.title("Dispersion 2D")
    plt.show()


#
# Layout
#
def initialize_layout(site_elements, a):
    Layout = {"2Dneighbours": np.zeros((site_elements, site_elements), dtype=int),
              "realLocation": [[[]]]}
    k = 1
    for i in range(site_elements):
        for j in range(site_elements):
            Layout["2Dneighbours"][j][i] = k
            Layout["realLocation"].append([k,[a*i,a*j]])
            k += 1
    Layout["realLocation"] = Layout["realLocation"][1:]

    return Layout

def plot_layout(Layout):
    X = np.zeros(len(Layout["realLocation"]))
    Y = np.zeros(len(Layout["realLocation"]))
    for i in range(len(Layout["realLocation"])):
        elem = Layout["realLocation"][i]
        X[i] = elem[1][0]
        Y[i] = elem[1][1]
    plt.scatter(X, Y, c='blue', linewidths=1)
    plt.gca().set_aspect('equal', 'box')
    plt.title("grid")
    plt.show()   



#
# Hamiltonian
#
def initialize_hamiltonian(Layout, t, PBC:bool, V_max, r_0):
    site_elements = Layout["2Dneighbours"].shape[0]
    elements = site_elements**2
    Hamiltonian = np.zeros((elements,elements))
    
    for i in range(site_elements):
        for j in range(site_elements):
            elem = Layout["2Dneighbours"][i,j] - 1

            next_x, prev_x, next_y, prev_y = neighbours(Layout["2Dneighbours"], i, j, PBC)

            if next_x in range(elements):
                Hamiltonian[next_x,elem] = -t
                Hamiltonian[elem,next_x] = -t

            if prev_x in range(elements):    
                Hamiltonian[prev_x,elem] = -t
                Hamiltonian[elem,prev_x] = -t

            if next_y in range(elements):
                Hamiltonian[next_y,elem] = -t
                Hamiltonian[elem,next_y] = -t
            
            if prev_y in range(elements):
                Hamiltonian[prev_y,elem] = -t
                Hamiltonian[elem,prev_y] = -t

    for i in range(elements):

       Hamiltonian[i][i] = 4*t + V(i+1, Layout, V_max, r_0)
       
    return Hamiltonian

def neighbours(Layout, i, j, PBC:bool):
    site_elements = Layout.shape[0]
    if PBC:
        #print("PBC True")
        if j in range(1,site_elements-1):
            next_y = j + 1
            prev_y = j - 1
        elif j == 0:
            next_y = j + 1
            prev_y = site_elements-1
        else:
            next_y = 0
            prev_y = j - 1

        if i in range(1,site_elements-1):
            next_x = i + 1
            prev_x = i - 1
        elif i == 0:
            next_x = i + 1
            prev_x = site_elements-1
        else:
            next_x = 0
            prev_x = i - 1
        
        return Layout[next_x,j]-1, Layout[prev_x,j]-1, Layout[i,next_y]-1, Layout[i,prev_y]-1 
    else:
        next_y = j + 1
        prev_y = j - 1
        next_x = i + 1
        prev_x = i - 1

        if prev_y in range(0,site_elements):
            prev_y_elem = Layout[i,prev_y]-1 
        else:
            prev_y_elem = -1
        if next_y in range(0,site_elements):
            next_y_elem = Layout[i,next_y]-1
        else:
            next_y_elem = -1
        if prev_x in range(0,site_elements):
            prev_x_elem =  Layout[prev_x,j]-1
        else:
            prev_x_elem = -1
        if next_x in range(0,site_elements):
            next_x_elem =  Layout[next_x,j]-1
        else:
            next_x_elem = -1

        return next_x_elem, prev_x_elem, next_y_elem, prev_y_elem 

def diagonalize_hamiltonian(Hamiltonian):
    #E,psiT = np.linalg.eigh(Hamiltonian) # This computes the eigen values and eigenvectors
    #psi = np.transpose(psiT)   # We take the transpose of psiT to the wavefunction vectors can accessed as psi[n]
    #plt.plot(E, 'o')
    #plt.show()
    return np.linalg.eigvalsh(Hamiltonian)

def diagonalize_hamiltonian_with_psi(Hamiltonian):
    E,psiT = np.linalg.eigh(Hamiltonian) # This computes the eigen values and eigenvectors
    psi = np.transpose(psiT)   # We take the transpose of psiT to the wavefunction vectors can accessed as psi[n]
    #E, psi = scipy.linalg.eigh(Hamiltonian)
    #plt.plot(E, 'o')
    #plt.show()
    return E, psi

def show_hamiltonian(Hamiltonian):
    Hamiltonian = np.flip(Hamiltonian,axis=1)
    ran = math.floor(math.sqrt(Hamiltonian.size))
    fig, ax = plt.subplots()
    pcm = ax.imshow(Hamiltonian, cmap='hot', extent=[0, ran, 0, ran])
    ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    ax.invert_yaxis()
    ax.xaxis.set_label_position('top')
    fig.colorbar(pcm, ax=ax)
    plt.show()

def plot_eigenvalues_hamiltonian_whole(diagonalizated_hamiltonain, t):
    plt.plot(diagonalizated_hamiltonain, '.')
    # bottom = np.amin(diagonalizated_hamiltonain)
    plt.hlines(8*t, 0, len(diagonalizated_hamiltonain), colors="red", linestyles="--")
    plt.hlines(0, 0, len(diagonalizated_hamiltonain), colors="red", linestyles="--")
    plt.xlabel("eigenvalue index")
    plt.ylabel("E (eV)")
    plt.title("Eigenvalues of Hamiltonian")
    plt.show()

def plot_eigenvalues_hamiltonian(diagonalizated_hamiltonain):
    plt.plot(diagonalizated_hamiltonain, 'o')
    plt.xlabel("eigenvalue index")
    plt.ylabel("E (eV)")
    plt.title("Eigenvalues of Hamiltonian")
    plt.show()

def plot_multiple_eigenvalues_hamiltonian(diagonalizated_hamiltonain_list, label_list, plot_first):
    if len(label_list) < len(diagonalizated_hamiltonain_list):
        label_list = ["Analitical"] + list(label_list)
    i=0
    for diagH in diagonalizated_hamiltonain_list:
        plt.plot(diagH[:plot_first], 'o', label=str(label_list[i]))
        i+=1
    plt.xlabel("eigenvalue index")
    plt.ylabel("E (eV)")
    plt.legend()
    plt.title("Comparison of different results.")
    plt.show()

def eigenvalues_hamiltonian(
        site_size,
        site_elements,
        PBC:bool,
        V_max, 
        r_0,
        m_eff,
        m_0,
        hbar,
        _1eV
):
    a = a_calc(site_size=site_size, site_elements=site_elements)
    t = t_calc(a=a, m_eff=m_eff, m_0=m_0, hbar=hbar, _1eV=_1eV)
    Layout = initialize_layout(site_elements, a)
    Hamiltonian = initialize_hamiltonian(Layout, t, PBC, V_max, r_0)
    return diagonalize_hamiltonian(Hamiltonian)

def comparison_eigenvalues(
        site_size,
        elements_list:list,
        PBC:bool,
        V_max, 
        r_0,
        m_eff,
        m_0,
        hbar,
        _1eV
):
    omega = omega_calc(V_max, r_0, m_eff, m_0,  _1eV)
    ham_list = [analitical_energies(11, omega, V_max, hbar, _1eV)]
    idx = 0 
    for i in elements_list:
        ham_list.append(eigenvalues_hamiltonian(
                site_size=site_size,
                site_elements=i,
                PBC=PBC,
                V_max=V_max, 
                r_0=r_0,
                m_eff=m_eff,
                m_0=m_0,
                hbar=hbar,
                _1eV=_1eV
        ))
        idx+=1
    return ham_list

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
    a = a_calc(site_size, site_elements)
    Layout = initialize_layout(site_elements, a)

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
    a = a_calc(site_size, site_elements)
    Layout = initialize_layout(site_elements, a)
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

#
# Analitical energies
#
def omega_calc(V, r, m_eff, m_0, _1eV): # in [J]
    return 4* np.sqrt((2*V*_1eV)/(m_eff*m_0*(r*10e-10)**2)) # 

def analitical_state(n_1, n_2, omega, V, hbar, _1eV):
    state_Jule = hbar * omega  * (n_1 + n_2 + 1) - V * _1eV
    state_eV = hbar * omega  * (n_1 + n_2 + 1) / _1eV - V
    state_eV_2 = state_Jule / _1eV
    # print(state_Jule, state_eV, state_eV_2)
    # hbar = np.pi * hbar / _1eV
    return state_eV # hbar * omega  * (n_1 + n_2 + 1) - V

def analitical_energies(N, omega, V, hbar, _1eV):
    n_x = N #math.floor(math.sqrt(N))
    #fill_more = N - n_x**2
    n = np.arange(0, n_x, 1)
    
    x = list(itertools.product(n, n))
    Energies = np.zeros(len(x))
    i = 0
    for elem in x:
        Energies[i] = analitical_state(elem[0], elem[1], omega, V, hbar, _1eV)
        i+=1
    return np.sort(Energies)

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

class TestMethods(unittest.TestCase):

    '''def test_df(self):
        data = {
        "calories": [420, 380, 390],
        "duration": [50, 40, 45],
        "a":[[100,101, 1.0],
        [500,200, 2.512562814070352],
        [3000,1000, 3.003003003003003]]
        }
        df = pd.DataFrame(data)
        df.to_csv('data1.csv', index=False)
        df_read = pd.read_csv('data1.csv')
        for elem in df_read.loc[:,"a"]: # data1, data2, value
            data1, data2, value= ast.literal_eval(elem)
            print(data1, data2, value, a_calc(data1, data2))'''


    def test_a(self):
        df_read = pd.read_csv('a_data.csv')
        for elem in df_read.loc[:,"a"]: # data1, data2, value
            data1, data2, value= ast.literal_eval(elem)
            self.assertEqual(value, a_calc(data1, data2))

    def test_t(self):
        self.assertEqual(t_calc(a_calc(100,101), m_eff, m_0, hbar, _1eV), 5.6865404718729256)
        self.assertEqual(t_calc(a_calc(500,200), m_eff, m_0, hbar, _1eV), 0.9007707569065586)
        self.assertEqual(t_calc(a_calc(3000,1000), m_eff, m_0, hbar, _1eV), 0.6305747863855169)



if __name__ == '__main__':
    # definitions of fundamental constant used in calculations
    site_size = 100
    site_elements = 101
    PBC = False
    V_max = 1 
    r_0 = 9
    m_eff = 0.067
    m_0 = scipy.constants.physical_constants["atomic unit of mass"][0]
    hbar = scipy.constants.hbar
    _1eV = scipy.constants.physical_constants["electron volt"][0]
    a = a_calc(site_size, site_elements)
    t = t_calc(a, m_eff, m_0, hbar, _1eV)
    a_2 = a_calc(500,200)
    a_3 = a_calc(3000,1000)
    t = t_calc(a, m_eff, m_0, hbar, _1eV)
    t_2 = t_calc(a_2, m_eff, m_0, hbar, _1eV)
    t_3 = t_calc(a_3, m_eff, m_0, hbar, _1eV)
    unittest.main()