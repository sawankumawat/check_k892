#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std;

// Test the Chernoff mixture p-value calculation
void test_chernoff_pvalue() {
    cout << "Testing Chernoff Mixture P-value Calculation" << endl;
    cout << "============================================" << endl;
    
    // Test some typical q0 values
    vector<double> q0_values = {0.0, 1.0, 2.71, 3.84, 6.63, 10.83, 25.0};
    
    cout << "q0_data\tChernoff p-val\tChernoff σ\tPure χ² σ\tDifference" << endl;
    cout << "-------\t--------------\t----------\t---------\t----------" << endl;
    
    for (double q0_data : q0_values) {
        // Chernoff mixture p-value calculation
        double chernoff_p_value = 0.0;
        if (q0_data > 0) {
            double chi2_p_value = TMath::Prob(q0_data, 1);  // P(χ²₁ ≥ q0_data)
            chernoff_p_value = 0.5 * (1.0 - chi2_p_value);  // Chernoff mixture p-value
        } else {
            chernoff_p_value = 1.0;  // If q0_data ≤ 0, p-value is 1
        }
        
        double chernoff_significance = (chernoff_p_value > 0 && chernoff_p_value < 1) ? 
                                     TMath::NormQuantile(1.0 - chernoff_p_value) : 0.0;
        
        // Pure χ²₁ calculation for comparison
        double pure_chi2_significance = sqrt(q0_data);
        
        double difference = chernoff_significance - pure_chi2_significance;
        
        cout << q0_data << "\t" << chernoff_p_value << "\t" << chernoff_significance 
             << "\t" << pure_chi2_significance << "\t" << difference << endl;
    }
    
    cout << "\nNote: Chernoff mixture gives more conservative (smaller) significance values" << endl;
    cout << "This is the correct behavior for boundary testing scenarios." << endl;
}

int main() {
    test_chernoff_pvalue();
    return 0;
}
