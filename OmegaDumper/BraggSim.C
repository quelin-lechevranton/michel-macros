#include "OmegaLight_tools.h"

/*
 * muon Bragg peak analysis preliminary for TrueMichel_v2_1.C
*/

const size_t n_file=3;
const vector<string> filelist = omega::ReadFileList(n_file,"list/jeremy.list");

const bool v = false;

const omega::Limits det = omega::cryostat;
// const omega::Limits det = omega::fiducial;

const omega::Binning bAvg = {30,0,30};
const omega::Binning bBragg = {100,0,10}; //#avg dEdx

size_t n_dep_bragg=50; // 0.03 cm/pt
const double dEdx_min_ratio = 1;
const size_t n_peak_tail=2;


void BraggSim(size_t i=0) {
    clock_t start_time=clock();

    // TH1D* hAvgdEdx = new TH1D ("hAvgdEdx","Average dEdx on last deposits;MeV/cm;#",bAvg.n,bAvg.min,bAvg.max);

    size_t nBragg=4;
    vector<TH1D*> hBragg(nBragg);
    vector<string> col = {"Bragg","Treshol","Unnormalized","Treshold Unnormalized"};
    for (int i=0; i<nBragg; i++) {
        stringstream t; t << "h" << i;
        stringstream n; n << col[i] << ";#AvgdEdx;#";
        hBragg[i] = new TH1D(
            t.str().c_str(), n.str().c_str(),
            bBragg.n, bBragg.min, bBragg.max
        );
    }

    vector<vector<size_t>> res(filelist.size());

    size_t n_mu=0;
    size_t n_mu_in=0;
    size_t n_mu_stop=0;

    size_t i_file=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << i_file+1 << "/" << n_file << ": " << filename << "\e[0m" << endl;

        omega::Truth T(filename.c_str());

        for (size_t i_evt=0; i_evt < T.N; i_evt++) {

            if (!v) cout << "Event#" << i_evt+1 << "/" << T.N << "\r" << flush;
            if (v) cout << "Event#" << i_evt+1 << "\r";

            T.GetEvt(i_evt);

            res[i_file].push_back(0);

            for (size_t i_prt=0; i_prt < T.Prt.N; i_prt++) {
                n_dep_bragg=50;

                T.GetPrt(i_prt);

                if (T.Prt.Pdg!=13 && T.Prt.Pdg!=-13) continue;

                T.GetPrtDep(i_prt);
                // if (T.Dep.N < n_dep_bragg+n_peak_tail) continue;
                if (T.Dep.N <= n_dep_bragg) continue;

                n_mu++;

                if (!omega::IsInside(
                    T.Dep.X,
                    T.Dep.Y,
                    T.Dep.Z,
                    det
                )) continue;
                n_mu_in++;


                double avg_body_dEdx=0;
                for (size_t i_dep=0; i_dep < T.Dep.N-n_dep_bragg; i_dep++) {
                    avg_body_dEdx += *T.Dep.E[i_dep] / 0.03;
                }
                avg_body_dEdx /= T.Dep.N - n_dep_bragg;

                const double dEdx_min = avg_body_dEdx*dEdx_min_ratio;

                double bragg_int=0;
                double bragg_treshold_int=0;
                size_t n_bragg_int=0;
                for (size_t i_dep=T.Dep.N-n_dep_bragg; i_dep < T.Dep.N; i_dep++) {
                    double dEdx = *T.Dep.E[i_dep] / 0.03;

                    bragg_int += dEdx;
                    if (dEdx < dEdx_min) continue;
                    bragg_treshold_int += dEdx;
                    n_bragg_int++;
                }

                // n_bragg_int=1;
                // n_dep_bragg=1;
                
                hBragg[0]->Fill(bragg_int / avg_body_dEdx / n_dep_bragg);
                hBragg[1]->Fill(bragg_treshold_int / avg_body_dEdx / n_bragg_int);
                hBragg[2]->Fill(bragg_int / 2 / n_dep_bragg);
                hBragg[3]->Fill(bragg_treshold_int / 2 / n_bragg_int);

                // ofstream("braggsim_results.txt")
                if (bragg_int/n_dep_bragg > 6) n_mu_stop++;

                res[i_file][i_evt] += 1;

                // double avg_dEdx=0;
                // for (size_t i_dep=n_dep-n_dep_bragg-n_peak_tail; i_dep < n_dep-n_peak_tail; i_dep++) {
                //     double E = (*T.DepE)[i_prt][i_dep];
                //     avg_dEdx += E/0.03/n_dep_bragg;
                // } //end deposit loop
                // hAvgdEdx->Fill(avg_dEdx);
                // cout << "\t\t\tavg_dEdx=" << avg_dEdx << endl;
            } //end particle loop
            if (v) cout << endl;
        } //end event loop
        i_file++;
    } //end file loop
    cout << endl;

    TCanvas* c1 = new TCanvas("c1","BraggSim");
    c1->Divide(2,2);
    for (int i=0; i<nBragg; i++) {
        c1->cd(i+1);
        hBragg[i]->Draw("hist");
    }

    cout << "n_mu / n_mu_in / n_mu_stop: " << n_mu << " / " << n_mu_in << " / " << n_mu_stop << endl;

    // ofstream f("braggsim_res.txt");
    // for (const vector<size_t> r : res) {
    //     for (const size_t s : r) {
    //         f << s << " ";
    //     }
    //     f << "\n";
    // }
    // f.close();

    cout << "total time of execution: " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}

void PlotBragg(size_t i_file=1, size_t i_evt=1) {
    clock_t start_time=clock();

    TGraph* gdEdx = new TGraph();
    gdEdx->SetTitle("muon Energy Loss;Range? (cm);dEdx (MeV/cm)");    

    string filename = filelist[i_file-1];
    
    cout << "\e[3mOpening file #" << i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

    omega::Truth T(filename.c_str());
    T.GetEvt(i_evt-1);

    for (size_t i_prt=0; i_prt < T.Prt.N; i_prt++) {

        T.GetPrt(i_prt);

        if (T.Prt.Pdg!=13 && T.Prt.Pdg!=-13) continue;

        T.GetPrtDep(i_prt);

        if (!omega::IsInside(
            T.Dep.X,
            T.Dep.Y,
            T.Dep.Z,
            det
        )) continue;
        if (T.Dep.N<n_dep_bragg+n_peak_tail) continue;

        for (size_t i_dep=T.Dep.N-n_dep_bragg-n_peak_tail; i_dep < T.Dep.N; i_dep++) {
            double E = *T.Dep.E[i_dep];
            double X = i_dep*0.03;
            gdEdx->AddPoint(X,E/0.03);
        } //end deposit loop
    } //end particle loop
    cout << endl;

    TCanvas* c2 = new TCanvas("c2","BraggSim");
    c2->cd();
    gdEdx->Draw();

    cout << "total time of execution: " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}