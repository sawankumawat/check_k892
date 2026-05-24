#ifndef HISTOGRAM_OPERATIONS_H
#define HISTOGRAM_OPERATIONS_H

#include <TH1.h>
#include <vector>

class HistogramOperations
{
public:
    HistogramOperations(); // Constructor declaration
    // std::vector<TH1 *> CalculateRatio(TH1D *hist1, std::vector<TH1 *> &variationHistograms);
    TH1D *CalculateRatio(TH1D *hist1, TH1 *variationHistograms);
    TH1D *RelativeUncertainty(TH1D *hist1, std::vector<TH1 *> &variationHistograms);
    TH1D *sigma(std::vector<TH1D *> &variationHistograms);
    TH1D *smooth(TH1D *hist1, int n = 1);
    TH1D *barlowcheck(TH1D *hist1, TH1 *hist2, bool &checkbar);
};

#endif
