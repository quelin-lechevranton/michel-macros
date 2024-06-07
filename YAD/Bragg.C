#include "YAD_tools.h"

/*
 * muon Bragg peak analysis preliminary for TrueMichel_v3.C
*/

vector<vector<double>> detlim={{-320, 350}, {-317,317}, {20,280}}; //cm

const size_t n_bragg_integration=15;
const size_t n_bragg_tail=2;

const size_t n_file=30;
const vector<string> filelist = yad::readFileList(n_file,"list/jeremy.list");

void Bragg(size_t i=0) {
    clock_t start_time=clock();

    TH1D* hAvgdEdx = new TH1D ("hAvgdEdx","Average dEdx on last deposits;MeV/cm;#",30,0,30);

    size_t i_file=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        yad::Truth T(filename.c_str());

        size_t n_evt = T.GetEntries();
        for (size_t i_evt=0; i_evt < n_evt; i_evt++) {

            cout << "Event#" << i_evt+1 << "/" << n_evt << "\r" << flush;

            T.GetEntry(i_evt);

            for (size_t i_prt=0; i_prt < T.NPrt; i_prt++) {

                int pdg = T.PrtPdg->at(i_prt);
                size_t n_ppt = T.PrtNPt->at(i_prt);
                size_t n_dep = T.PrtNDep->at(i_prt);

                if (pdg!=13 && pdg!=-13) continue;
                if (!yad::isInside(
                    T.DepX->at(i_prt),
                    T.DepY->at(i_prt),
                    T.DepZ->at(i_prt),
                    detlim[0][0], detlim[0][1],
                    detlim[1][0], detlim[1][1],
                    detlim[2][0], detlim[2][1]
                )) continue;
                if (n_dep<n_bragg_integration+n_bragg_tail) continue;

                // for (size_t i_ppt=0; i_ppt < T.PrtNPt->at(i_prt); i_ppt++) {

                // } //end particlepoint loop

                double avg_dEdx=0;
                for (size_t i_dep=n_dep-n_bragg_integration-n_bragg_tail; i_dep < n_dep-n_bragg_tail; i_dep++) {
                    double E = (*T.DepE)[i_prt][i_dep];
                    avg_dEdx += E/0.03/n_bragg_integration;
                } //end deposit loop
                hAvgdEdx->Fill(avg_dEdx);
                // cout << "\t\t\tavg_dEdx=" << avg_dEdx << endl;
            } //end particle loop
        } //end event loop
    } //end file loop
    cout << endl;

    TCanvas* c1 = new TCanvas("c1","Bragg");
    c1->cd();
    hAvgdEdx->Draw("hist");

    cout << "total time of execution: " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}

void PlotBragg(size_t i_file=1, size_t i_evt=0) {
    clock_t start_time=clock();

    TGraph* gdEdx = new TGraph();
    // gdEdx->SetTitle("muon Energy Loss;Deposit Number;dEdx (Mev/cm)");    

    string filename = filelist[i_file-1];
    
    cout << "\e[3mOpening file #" << i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

    yad::Truth T(filename.c_str());
    T.GetEntry(i_evt-1);

    for (size_t i_prt=0; i_prt < T.NPrt; i_prt++) {

        int pdg = T.PrtPdg->at(i_prt);
        size_t n_ppt = T.PrtNPt->at(i_prt);
        size_t n_dep = T.PrtNDep->at(i_prt);

        if (pdg!=13 && pdg!=-13) continue;
        if (!yad::isInside(T.DepX->at(i_prt),T.DepY->at(i_prt),T.DepZ->at(i_prt))) continue;
        if (n_dep<n_bragg_integration+n_bragg_tail) continue;

        for (size_t i_dep=n_dep-n_bragg_integration-n_bragg_tail; i_dep < n_dep-n_bragg_tail; i_dep++) {
            double E = (*T.DepE)[i_prt][i_dep];
            gdEdx->AddPoint(i_dep,E);
        } //end deposit loop
    } //end particle loop
    cout << endl;

    TCanvas* c1 = new TCanvas("c1","Bragg");
    c1->cd();
    gdEdx->Draw();

    cout << "total time of execution: " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}