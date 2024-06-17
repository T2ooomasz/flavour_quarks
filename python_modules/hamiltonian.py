import python_modules.base as base
import python_modules.layout as layout
import python_modules.potential as potential

import numpy as np
import math
import matplotlib.pyplot as plt
import itertools

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

       Hamiltonian[i][i] = 4*t + potential.V(i+1, Layout, V_max, r_0)
       
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
        label_list = ["Analitical"] + label_list
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
    a = base.a_calc(site_size=site_size, site_elements=site_elements)
    t = base.t_calc(a=a, m_eff=m_eff, m_0=m_0, hbar=hbar, _1eV=_1eV)
    Layout = layout.initialize_layout(site_elements, a)
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
# Analitical energies
#
def omega_calc(V, r, m_eff, m_0, _1eV): # in [J]
    return np.sqrt((2*V*_1eV)/(m_eff*m_0*(r*1e-10)**2)) # 

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
