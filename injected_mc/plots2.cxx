#include <iostream>
using namespace std;
#include "style.h"

void plots2()
{
    TFile *f = new TFile("f0(1710)_gen.root", "read");
    if (f->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    
}