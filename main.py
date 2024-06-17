import python_modules.initial_values as init
import python_modules.base as base
import python_modules.layout as layout
import python_modules.potential as potential
import python_modules.hamiltonian as hamiltonian
import python_modules.wave_function as psi

import scipy

def main():
    # definitions of fundamental constant used in calculations
    QD = init.QD_first
    site_size = QD["site_size"]
    site_elements = QD["site_elements"]
    PBC = QD["PBC"]
    V_max = QD["V_max"] 
    r_0 = QD["r_0"]
    m_eff =QD["m_eff"]
    m_0 = scipy.constants.physical_constants["atomic unit of mass"][0]
    hbar = scipy.constants.hbar
    _1eV = scipy.constants.physical_constants["electron volt"][0]
    a = base.a_calc(site_size, site_elements)
    t = base.t_calc(a, m_eff, m_0, hbar, _1eV)

    L = layout.initialize_layout(site_elements,a)
    potential.plot_potential_on_grid_heat(L, site_size, site_elements, V_max, r_0)
    H = hamiltonian.initialize_hamiltonian(L, t, PBC, V_max, r_0)
    diagH, psi_ = hamiltonian.diagonalize_hamiltonian_with_psi(H)
    #hamiltonian.plot_eigenvalues_hamiltonian_whole(diagH, t)
    #hamiltonian.plot_eigenvalues_hamiltonian(diagH[:40])
    psi.psi_n_3D(psi_, 0, site_size, site_elements)

if __name__ == "__main__":
    main()