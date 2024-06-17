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
    return (hbar*hbar) / (2*m_eff*m_0*a*a*1e-20*_1eV) # [eV]

def a_calc(
        site_size,
        site_elements
):
    return site_size / (site_elements - 1) # [Angstrom] - distance between elements 
