#include "YAD_tools.h"

/*
 * muon Bragg peak analysis preliminary for TrueMichel_v2_1.C
*/

vector<vector<double>> detlim={{-350, 350}, {-350,350}, {0,300}}; //cm
// vector<vector<double>> detlim={{-320, 350}, {-317,317}, {20,280}}; //cm

const size_t n_bragg_integration=50;
const size_t n_bragg_tail=1;

const size_t n_file=3;
const vector<string> filelist = yad::readFileList(n_file,"list/ijclab.list");

void BraggCalo(size_t i=0) {
    clock_t start_time=clock();

    TH1D* hAvgdEdx = new TH1D ("hAvgdEdx","Average dEdx on last deposits;MeV/cm;#",30,0,30);

    size_t i_file=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        // yad::Truth T(filename.c_str());
        yad::Reco R(filename.c_str());

        size_t n_evt = R.GetEntries();
        for (size_t i_evt=0; i_evt < n_evt; i_evt++) {

            cout << "Event#" << i_evt+1 << "/" << n_evt << "\r" << flush;

            // T.GetEntry(i_evt);
            R.GetEntry(i_evt);
            

            for (size_t i_trk=0; i_trk < R.NTrk; i_trk++) {

                // int pdg = T.PrtPdg->at(i_prt);
                size_t n_tpt = R.TrkNPt->at(i_trk);
                size_t n_cal = R.TrkCalNPt->at(i_trk);

                // if (pdg!=13 && pdg!=-13) continue;
                // if (!yad::isInside(
                //     R.TrkPtX->at(i_trk),
                //     R.TrkPtY->at(i_trk),
                //     R.TrkPtZ->at(i_trk),
                //     detlim[0][0], detlim[0][1],
                //     detlim[1][0], detlim[1][1],
                //     detlim[2][0], detlim[2][1]
                // )) continue;
                if (n_cal<n_bragg_integration+n_bragg_tail) continue;

                double avg_dEdx=0;
                for (size_t i_cal=n_cal-n_bragg_integration-n_bragg_tail; i_cal < n_cal-n_bragg_tail; i_cal++) {
                    double dEdx = (*R.TrkCaldEdx)[i_trk][i_cal];
                    avg_dEdx += dEdx/n_bragg_integration;
                } //end deposit loop
                hAvgdEdx->Fill(avg_dEdx);
                // cout << "\t\t\tavg_dEdx=" << avg_dEdx << endl;
            } //end particle loop
        } //end event loop
    } //end file loop
    cout << endl;

    TCanvas* c1 = new TCanvas("c1","BraggCalo");
    c1->cd();
    hAvgdEdx->Draw("hist");

    cout << "total time of execution: " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}

void PlotBragg(size_t i_file=1, size_t i_evt=1) {
    clock_t start_time=clock();

    TGraph* gdEdx = new TGraph();
    gdEdx->SetTitle("muon Energy Loss;Range (cm);dEdx (MeV/cm)");    

    // TGraph* gResRange = new TGraph();
    // gResRange->SetTitle(";Calo point number;Residual Range (cm)");    

    string filename = filelist[i_file-1];
    
    cout << "\e[3mOpening file #" << i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

    // yad::Truth T(filename.c_str());
    yad::Reco R(filename.c_str());
    
    // T.GetEntry(i_evt-1);
    R.GetEntry(i_evt-1);

    for (size_t i_trk=0; i_trk < R.NTrk; i_trk++) {

        // int pdg = T.PrtPdg->at(i_prt);
        size_t n_tpt = R.TrkNPt->at(i_trk);
        size_t n_cal = R.TrkCalNPt->at(i_trk);

        // if (pdg!=13 && pdg!=-13) continue;
        if (!yad::isInside(
            R.TrkPtX->at(i_trk),
            R.TrkPtY->at(i_trk),
            R.TrkPtZ->at(i_trk),
            detlim[0][0], detlim[0][1],
            detlim[1][0], detlim[1][1],
            detlim[2][0], detlim[2][1]
        )) continue;
        if (n_cal<n_bragg_integration+n_bragg_tail) continue;

        for (size_t i_cal=n_bragg_tail; i_cal < n_bragg_integration + n_bragg_tail; i_cal++) {
            double dEdx = (*R.TrkCaldEdx)[i_trk][i_cal];
            double X = (*R.TrkCalRange)[i_trk]-(*R.TrkCalResRange)[i_trk][i_cal];
            gdEdx->AddPoint(X,dEdx);
        } //end deposit loop

        // for (size_t i_cal=0; i_cal < n_cal; i_cal++) {
        //     double ResRange = (*R.TrkCalResRange)[i_trk][i_cal];
        //     gResRange->AddPoint(i_cal,ResRange);
        // } //end deposit loop
    } //end particle loop
    cout << endl;

    TCanvas* c2 = new TCanvas("c2","BraggCalo");
    c2->cd();
    gdEdx->Draw();

    // TCanvas* c3 = new TCanvas("c3","BraggCalo");
    // c3->cd();
    // gResRange->Draw();

    cout << "total time of execution: " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}