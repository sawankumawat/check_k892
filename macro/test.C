#include <TCanvas.h>
#include <TGraph.h>
#include <TPad.h>

void test() {
    // Example data for the first graph
    const int nPoints1 = 5;
    double x1[nPoints1] = {0, 1, 2, 3, 4};
    double y1[nPoints1] = {0, 2, 4, 6, 8};

    // Example data for the second graph
    const int nPoints2 = 5;
    double x2[nPoints2] = {0, 1, 2, 3, 4};
    double y2[nPoints2] = {8, 6, 4, 2, 0};

    // Create the first graph
    TGraph* graph1 = new TGraph(nPoints1, x1, y1);

    // Create the second graph
    TGraph* graph2 = new TGraph(nPoints2, x2, y2);

    // Create a canvas to hold the two graphs
    TCanvas* canvas = new TCanvas("canvas", "Two Graphs", 800, 600);

    // Create two pads on the canvas: one for each graph
    TPad* pad1 = new TPad("pad1", "pad1", 0, 0.5, 1, 1);
    TPad* pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.5);

    // Divide the canvas into two pads vertically (one below the other)
    pad1->SetBottomMargin(0); // Remove the bottom margin for the first pad
    pad1->Draw();
    pad2->SetTopMargin(0);    // Remove the top margin for the second pad
    pad2->Draw();

    // Set the first pad as active
    pad1->cd();

    // Draw the first graph on the first pad
    graph1->SetTitle("Graph 1;X-axis;Y-axis");
    graph1->Draw("alp");

    // Set the second pad as active
    pad2->cd();

    // Draw the second graph on the second pad
    graph2->SetTitle("Graph 2;X-axis;Y-axis");
    graph2->Draw("alp");

    // Update the canvas
    // canvas->Update();
}

