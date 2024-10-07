#include <iostream>
#include <tuple>
#include <vector>
#include <algorithm>
#include "../src/fitfunc.h"
using namespace std;

int check()
{
    // // //checking the usage of tuple
    // vector<tuple<int, int, int, int>> check;

<<<<<<< HEAD
    int arr1[5] = {1, 2, 3, 4, 5};      // First element of the tuple
    int arr2[5] = {10, 20, 30, 40, 50}; // Second element of the tuple
    int arr3[5] = {9, 8, 8, 2, 6};      // Third element of the tuple
    int arr4[5] = {5, 10, 50, 20, 25};  // Fourth element of the tuple
=======
    // int arr1[5] = {1, 2, 3, 4, 5};      // First element of the tuple
    // int arr2[5] = {10, 20, 30, 40, 50}; // Second element of the tuple
    // int arr3[5] = {9, 8, 8, 2, 6};      // Third element of the tuple
    // int arr4[5] = {5, 10, 15, 20, 25};  // Fourth element of the tuple
>>>>>>> a9da7fb (added the coherent BW fit with the real and imaginary parts)

    // // Create tuples using elements at corresponding positions in the arrays
    // for (int i = 0; i < 5; i++)
    // {
    //     check.push_back(make_tuple(arr1[i], arr2[i], arr3[i], arr4[i]));
    // }

<<<<<<< HEAD
    // sort in asceding order w.r.t to the third array i.e. chi2/NDF
    sort(check.begin(), check.end(),
         [](const auto &a, const auto &b)
         {
             return std::get<3>(a) < std::get<3>(b);
         });
=======
    // // sort in asceding order w.r.t to the third array i.e. chi2/NDF
    // sort(check.begin(), check.end(),
    //      [](const auto &a, const auto &b)
    //      {
    //          return std::get<2>(a) < std::get<2>(b);
    //      });
>>>>>>> a9da7fb (added the coherent BW fit with the real and imaginary parts)

    // for (int i = 0; i < 5; i++)
    // {
    //     float best_ipar1 = std::get<0>(check[i]);
    //     float best_ipar2 = std::get<1>(check[i]);
    //     float best_ipar3 = std::get<2>(check[i]);
    //     float best_chi2ndf = std::get<3>(check[i]);
    //     cout << "ipar1: " << best_ipar1 << ",  ipar2: " << best_ipar2 << ",  ipar3: " << best_ipar3 << ",  chi2/NDF: " << best_chi2ndf << endl;
    // }

    //checking the fitting function expol
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    TF1 *expol = new TF1("expol", exponential_bkg, 0, 10, 3);
    expol->SetParameter(0, -1000); // This tells the curvature. Positive and negative values has opposite curvatures. Magnitude of this parameters scales the overall curve
    expol->SetParameter(1, 0);
    expol->SetParameter(2, -0.5); // This parameter tells the shape. Positive value means exponential decay and negative value means exponential rise. Larger the value, faster the decay or rise. 0 means flat line.
    expol->Draw();

    return 0;
}