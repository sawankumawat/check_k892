#include <iostream>
using namespace std;

float f1270Mass = 1275.4;
float f1270MassErr = 0.8;
float f1270Width = 185.8;
float f1270WidthErr = 2.8;

float a1320Mass = 1318.2;
float a1320MassErr = 0.6;
float a1320Width = 107;
float a1320WidthErr = 5;

float f1525Mass = 1517.3;
float f1525MassErr = 2.4;
float f1525Width = 84.4;
float f1525WidthErr = 2.7;

float f1710Mass = 1733;
float f1710MassErr = 8;
float f1710Width = 150;
float f1710WidthErr = 12;

void calculate_sigma()
{
    float pdg_value = f1270Mass;
    float pdg_error = f1270MassErr;
    float value, valueErr;
    cin >> value >> valueErr;

    // cout<< "value = " << value << " valueErr = " << valueErr << endl;
    double sigma = abs(value - pdg_value) / sqrt(valueErr * valueErr + pdg_error * pdg_error);

    cout << Form("sigma = %.2f", sigma) << endl;

    // f1270 mass sigma 14.7, 3.1, 3.7
    // f1270 width sigma 2.9, 0.6, 10.3
    // a1320 mass sigma 2.9, 1.1, 7.4
    // a1320 width sigma 10.7, 0.4, 0.1
    // f1525 mass sigma 0.4, 0.4, 2.2
    // f1525 width sigma 0.1, 0.8, 0.6
    // f1710 mass sigma 2.4, 2.3, 0.9
    // f1710 width sigma 0.1, 0.1, 3.1
}
