#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <ROOT/RDataFrame.hxx>
#include <string>
#include <vector>

std::string filename="pdvd_Mitrees_out.root",
            micheltreename="checks/Michel",
            eventtreename="checks/Event";

struct Binning {
    int n;
    double min, max;
};

std::vector<Binning> binChan = {{1175, 1900, 3075},
                                {1175, 4975, 6150},
                                {1175, 8045, 9220},
                                {1175, 11115, 12290}}; 

Binning binTick = {2000, 0, 6000};

ROOT::RDataFrame michel(micheltreename, filename),
                 event(eventtreename, filename);

TFile* f = new TFile("test4_cpp.root", "RECREATE");

size_t n_event = 10;
size_t n_section = 4;
size_t n_channel = 12288;

size_t GetSection(size_t ch) {
    return size_t(4.*ch / n_channel);
}

void test4() {

    std::vector<TCanvas*> c(n_event);

    size_t i_event = 0;
    event.Range(n_event).Foreach( [&c, &i_event](
        size_t nhit, 
        std::vector<float> &ticks, 
        std::vector<size_t> &chans, 
        std::vector<float> &adcs, 
        size_t nmichel, 
        std::vector<size_t> &imichel
    ){
        std::vector<TH2F*> h2(n_section);

        for (size_t s=0; s<n_section; s++) {
            h2[s] = new TH2F(
                Form("h2_%ld_%ld", i_event, s),
                Form("Event %ld, Section %ld", i_event, s),
                binChan[s].n, binChan[s].min, binChan[s].max,
                binTick.n, binTick.min, binTick.max
            );
        }

        for (size_t i=0; i<nhit; i++) {
            h2[GetSection(chans[i])]->Fill(chans[i], ticks[i], adcs[i]);
        }


        std::vector<TGraph*> g(nmichel);
        int i_michel = 0;
        michel.Range(imichel.front(), imichel.back()+1).Foreach( [&g, &i_michel](
            int michel_nhit, 
            std::vector<float> &michel_ticks, 
            std::vector<int> &michel_chans
        ) {
            g[i_michel] = new TGraph();

            for (int i=0; i<michel_nhit; i++) {
                g[i_michel]->AddPoint(michel_ticks[i], michel_chans[i]);
            }
            i_michel++;
        }, {"MCPNHit", "MCPHitTick", "MCPHitChan"} );

        c[i_event] = new TCanvas(
            Form("c%ld", i_event),
            Form("Event %ld", i_event)
        );

        c[i_event]->Divide(2, 2);
        for (size_t s=0; s<n_section; s++) {
            auto pad = c[i_event]->cd(s+1);

            pad->DrawFrame(
                binChan[s].min, binTick.min,
                binChan[s].max, binTick.max,
                Form("Section %ld, W plane;Tick;Channel",s)
            );
            pad->SetLogz();

            // std::cout << "evt#" << i_event << " sec#" << s << " nhit=" << h2[s]->GetEntries() << std::endl;

            pad->cd();
            h2[s]->SetMinimum(0.1);
            h2[s]->Draw("colz same");
            // pad->Update();

            for (size_t i=0; i<nmichel; i++) {
                if (g[i]->GetN() == 0) continue;
                g[i]->SetEditable(kFALSE);
                g[i]->SetMarkerStyle(kFullDiamond);
                g[i]->SetMarkerSize(2.5);
                g[i]->SetMarkerColor(kPink-5);
                g[i]->Draw("same P");
            }
        }

        c[i_event]->Write();

        i_event++;

    }, {"NHit", "HitTick", "HitChan", "HitADC", "NMichel", "iMichel"} );

    f->Close();
}

