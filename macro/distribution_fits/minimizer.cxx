#include <iostream>
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

using namespace std;

//------------------------------------------------------
// Function to minimize
//------------------------------------------------------

double MyFunction(const double *par)
{
    double a = par[0];
    double b = par[1];

    return pow(a-3.0,2) + pow(b+2.0,2);
}

//------------------------------------------------------

void minimizer()
{
    // Create Minuit2 minimizer using Migrad algorithm
    ROOT::Math::Minimizer *min =
        ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");

    // Tell ROOT which function to minimize
    ROOT::Math::Functor f(&MyFunction,2);
    min->SetFunction(f);

    // Set parameters
    // name , initial value , step size

    min->SetVariable(0,"a",0.0,0.01);
    min->SetVariable(1,"b",0.0,0.01);

    // Perform minimization
    min->Minimize();

    //--------------------------------------------------
    // Results
    //--------------------------------------------------

    const double *x = min->X();

    cout<<"Best fit parameters\n";
    cout<<"a = "<<x[0]<<endl;
    cout<<"b = "<<x[1]<<endl;

    cout<<endl;

    cout<<"Errors"<<endl;
    cout<<"Error(a) = "<<min->Errors()[0]<<endl;
    cout<<"Error(b) = "<<min->Errors()[1]<<endl;

    cout<<endl;

    cout<<"Minimum function value = "
        <<min->MinValue()<<endl;

    delete min;
}