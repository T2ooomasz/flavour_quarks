#ifndef layout_H
#define layout_H

class Layout {
    public:
        int layout_size_x;
        int layout_size_y;
        int **layout;

        Layout(int x);
        Layout(int x, int y);
        ~Layout();
        void initializeLayout();
        void fillLayout();
        void showLayout();
};

#endif