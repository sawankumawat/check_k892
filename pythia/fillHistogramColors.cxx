
void fillHistogramColors()
{
    // Create histogram
    TH1F *h = new TH1F("h", "Example", 100, 0, 100);

    // Fill histogram (example: flat distribution)
    for (int i = 0; i < 10000; i++)
        h->Fill(gRandom->Uniform(0, 100));

    // Total entries
    double total = h->Integral();

    // Number of regions (10% each → 10 regions)
    const int nRegions = 10;

    // Colors
    int colors[nRegions] = {kRed, kOrange, kYellow, kGreen, kCyan,
                            kBlue, kViolet, kMagenta, kPink, kGray};

    // Clone histograms
    TH1F *hRegions[nRegions];
    for (int i = 0; i < nRegions; i++)
    {
        hRegions[i] = (TH1F *)h->Clone(Form("h_%d", i));
        hRegions[i]->Reset();
        hRegions[i]->SetFillColor(colors[i]);
        hRegions[i]->SetLineColor(kBlack);
    }

    // Assign bins to percentile regions
    double cumulative = 0;
    int currentRegion = 0;

    for (int bin = 1; bin <= h->GetNbinsX(); bin++)
    {
        double content = h->GetBinContent(bin);
        cumulative += content;

        double frac = cumulative / total;

        int region = frac * nRegions;
        if (region >= nRegions)
            region = nRegions - 1;

        hRegions[region]->SetBinContent(bin, content);
    }

    // Draw
    TCanvas *c = new TCanvas("c", "c", 800, 600);

    h->SetLineColor(kBlack);
    h->Draw("HIST");

    for (int i = 0; i < nRegions; i++)
    {
        hRegions[i]->Draw("HIST SAME");
    }
}