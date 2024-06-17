#include "layout.h"
#include<iostream>

Layout::Layout(int x){
    layout_size_x = x;
    layout_size_y = x;
    initializeLayout();
}

Layout::Layout(int x, int y) {
    layout_size_x = x;
    layout_size_y = y;
    initializeLayout();
}

Layout::~Layout() {
    for(int i = 0; i < layout_size_y; ++i) 
        delete [] layout[i];
    delete [] layout;
}

void Layout::initializeLayout() {
    layout = new int*[layout_size_x];
    for(int i = 0; i < layout_size_x; i++)
        layout[i] = new int[layout_size_y];
    fillLayout();
}

void Layout::fillLayout(){
    int c = 0;
    for (int i = 0; i < layout_size_x; i++) {
        for (int j = 0; j < layout_size_y; j++) {
            layout[i][j] = ++c;
        }
    }
}

void Layout::showLayout(){
    for (int i = 0; i < layout_size_x; i++) {
        for (int j = 0; j < layout_size_y; j++) {
            std::cout << layout[i][j] << " ";
        }
        std::cout << std::endl;
    }
}