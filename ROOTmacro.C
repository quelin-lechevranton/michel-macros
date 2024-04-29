#include <TGraph.h>
#include <TCanvas.h>


void macro() {
    auto graph = new TGraph();
    graph->SetTitle("test");
    graph->GetXaxis()->SetName("X");
    graph->SetPoint(1,2,3);
    graph->SetPoint(2,4,5);
    auto canvas = new TCanvas("c","TEST");
    canvas->cd();
    graph->Draw("AP");
}