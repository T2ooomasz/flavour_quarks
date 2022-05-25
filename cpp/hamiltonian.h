#include "layout.h"

#ifndef hamiltonian_H
#define hamiltonian_H

class Hamiltonian {
    public:
        int **H;
        
        Hamiltonian(Layout X);
        ~Hamiltonian();
};

#endif