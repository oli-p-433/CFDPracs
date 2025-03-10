#include "fluid.H"
#include "operators.H"

class boundary{
    public:
        void reflectiveLeftBC(fluid&);
        void reflectiveRightBC(fluid&);
        void reflectiveTopBC(fluid&);
        void reflectiveBottomBC(fluid&);

        void transmissiveLeftBC(fluid&);
        void transmissiveRightBC(fluid&);
        void transmissiveTopBC(fluid&);
        void transmissiveBottomBC(fluid&);

        void updateBottomLeftCorner(fluid&,bool,bool);
        void updateTopRightCorner(fluid&,bool,bool);
        void updateBottomRightCorner(fluid&,bool,bool);
        void updateTopLeftCorner(fluid&,bool,bool);

        boundary(int,int,int);

    private:
        int nCellsX;
        int nCellsY;
        int nGhost;
};
