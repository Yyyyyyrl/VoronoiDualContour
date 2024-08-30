#include "debug.h"


void print_cell(Delaunay::Cell c)
{
    using namespace std;

    cerr << "Cell: [";
    for (int i = 0; i < 4; i++)
    {
        cerr << "(" << c.vertex(i)->point() << ")";
        if (i < 3)
        {
            cerr << ",";
        }
    }
    cerr << "]" << endl;
}

void print_facet(Facet f)
{
    int iFacet = f.second;
    int d1, d2, d3;
    d1 = (iFacet + 1) % 4;
    d2 = (iFacet + 2) % 4;
    d3 = (iFacet + 3) % 4;
    Cell_handle c = f.first;
    std::cout << "Facet: " << c->vertex(d1)->point() << ", " << c->vertex(d2)->point() << ", " << c->vertex(d3)->point() << std::endl;
    std::cout << "ifacet: " << iFacet << std::endl;
}


