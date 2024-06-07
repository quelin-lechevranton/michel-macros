#include "YAD_tools.h"

const size_t n_file=1;
const vector<string> filelist = yad::readFileList(n_file,"list/jeremy.list");

void PrtDep_Pt() {
    clock_t start_time=clock();

    /* graph declarations */
    
    TH1D* hPrtDist = new TH1D("hPrtDist","Distance between two PrtPt;distance (cm);#",100,0,1);
    double avg_prt_distance=0;
    TH1D* hDepDist = new TH1D("hDepDist","Distance between two PrtDep;distance (cm);#",100,0,1);
    double avg_dep_distance=0;
    TH1D* hn_depPt = new TH1D("hn_depPt","ratio n_dep/NPt;ratio;#",100,0,100);

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
                if (pdg!=13 && pdg!=-13) continue;

                size_t n_ppt = T.PrtNPt->at(i_prt);
                size_t n_dep = T.Prtn_dep->at(i_prt);

                hn_depPt->Fill((double) n_dep/n_ppt);

                double X0 = (*T.PrtX)[i_prt][0];
                double Y0 = (*T.PrtY)[i_prt][0];
                double Z0 = (*T.PrtZ)[i_prt][0];
                for (size_t i_ppt=1; i_ppt < n_ppt; i_ppt++) {
                    double X = (*T.PrtX)[i_prt][i_ppt];
                    double Y = (*T.PrtY)[i_prt][i_ppt];
                    double Z = (*T.PrtZ)[i_prt][i_ppt];

                    hPrtDist->Fill(yad::distance(X0,Y0,Z0,X,Y,Z));
                    avg_prt_distance+=yad::distance(X0,Y0,Z0,X,Y,Z);

                    X0=X; Y0=Y; Z0=Z;
                } //end particlepoint loop
                if (n_ppt>1) avg_prt_distance/=(n_ppt-1);


                X0 = (*T.DepX)[i_dep][0];
                Y0 = (*T.DepY)[i_dep][0];
                Z0 = (*T.DepZ)[i_dep][0];
                for (size_t i_dep=0; i_dep < n_dep; i_dep++) {
                    double X = (*T.DepX)[i_prt][i_dep];
                    double Y = (*T.DepY)[i_prt][i_dep];
                    double Z = (*T.DepZ)[i_prt][i_dep];

                    hDepDist->Fill(yad::distance(X0,Y0,Z0,X,Y,Z));
                    avg_dep_distance+=yad::distance(X0,Y0,Z0,X,Y,Z);

                    X0=X; Y0=Y; Z0=Z;

                } //end deposit loop
                if (n_dep>1) avg_dep_distance/=(n_dep-1);
            } //end particle loop

        } //end event loop
    } //end file loop
    cout << endl;

    cout << avg_prt_distance << endl;
    TCanvas* c1 = new TCanvas("c1","PrtDep_Pt");
    c1->Divide(2,1);
    c1->cd(1);
    hPrtDist->Draw("hist");
    c1->cd(2);
    hDepDist->Draw("hist");
    // hn_depPt->Draw("hist");

    cout << "total time of execution: " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}