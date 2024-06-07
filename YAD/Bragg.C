#include "YAD_tools.h"

const size_t n_bragg_integration=15;
const size_t n_bragg_tail=2;

const size_t n_file=1;
const vector<string> filelist = yad::readFileList(n_file,"list/jeremy.list");

void Bragg(size_t i=0) {
    clock_t start_time=clock();

    TGraph* gdEdx = new TGraph();
    // gdEdx->SetTitle("muon Energy Loss;Deposit Number;dEdx (Mev/cm)");    

    size_t i_file=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        yad::Truth T(filename.c_str());

        size_t n_evt = T.GetEntries();
        for (size_t i_evt=0; i_evt < n_evt; i_evt++) {
        // size_t i_evt=i; {

            cout << "Event#" << i_evt+1 << "/" << n_evt << "\r" << flush;

            T.GetEntry(i_evt);

            for (size_t i_prt=0; i_prt < T.NPrt; i_prt++) {

                int pdg = T.PrtPdg->at(i_prt);
                size_t n_ppt = T.PrtNPt->at(i_prt);
                size_t n_dep = T.PrtNDep->at(i_prt);

                if (pdg!=13 && pdg!=-13) continue;
                if (!yad::isInside(T.DepX->at(i_prt),T.DepY->at(i_prt),T.DepZ->at(i_prt))) continue;
                if (n_dep<n_bragg_integration+n_bragg_tail) continue;

                // cout << i_evt << endl;

                // for (size_t i_ppt=0; i_ppt < T.PrtNPt->at(i_prt); i_ppt++) {

                // } //end particlepoint loop

                double avg_dEdx=0;
                for (size_t i_dep=n_dep-n_bragg_integration-n_bragg_tail; i_dep < n_dep-n_bragg_tail; i_dep++) {
                    // gdEdx->AddPoint(n_dep-i_dep,(*T.DepE)[i_prt][i_dep]/0.03);
                    double E = (*T.DepE)[i_prt][i_dep];
                    avg_dEdx += E/0.03/n_bragg_integration;
                    gdEdx->AddPoint(i_dep,E);
                } //end deposit loop
                cout << endl << "avg_dEdx=" << avg_dEdx << endl;
            } //end particle loop
        } //end event loop
    } //end file loop
    cout << endl;

    // TCanvas* c1 = new TCanvas("c1","Bragg");
    // c1->cd();
    // gdEdx->Draw();

    cout << "total time of execution: " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}