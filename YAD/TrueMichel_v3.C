#include "YAD_tools.h"

vector<vector<double>> detlim={{-320, 350}, {-317,317}, {20,280}}; //cm

const size_t n_file=30;
const vector<string> filelist = yad::readFileList(n_file,"list/jeremy.list");


void TrueMichel_v3() {
    clock_t start_time=clock();

    TColor color;

    size_t  N_evt=0,
            N_prt=0;
    size_t i_file=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        yad::Truth T(filename.c_str());

        size_t n_evt = T.GetEntries();
        N_evt+=n_evt;
        for (size_t i_evt=0; i_evt < n_evt; i_evt++) {

            cout << "Event#" << i_evt+1 << "/" << n_evt << "\r" << flush;

            T.GetEntry(i_evt);

            // N_prt+=T.NPrt;

            for (size_t i_prt=0; i_prt < T.NPrt; i_prt++) {
                if (T.PrtMomID->at(i_prt)!=-1) continue;
                if (T.PrtPdg->at(i_prt)!=13 && T.PrtPdg->at(i_prt)!=-13) continue;
                if (yad::isInside(
                    T.DepX->at(i_prt),
                    T.DepY->at(i_prt),
                    T.DepZ->at(i_prt),
                    detlim[0][0], detlim[0][1],
                    detlim[1][0], detlim[1][1],
                    detlim[2][0], detlim[2][1]
                )) {
                    N_mu_inside++;
                }
            }

            for (size_t i_prt=0; i_prt < T.NPrt; i_prt++) {

                size_t n_el_dep = T.PrtNDep->at(i_prt);
            } //end particle loop
        } //end event loop
    } //end file loop

    TCanvas* c1 = new TCanvas("c1","TrueMichel_v3");
    c1->cd();

    cout << N_evt << " events treated in " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}