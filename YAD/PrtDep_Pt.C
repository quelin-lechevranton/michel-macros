#include "YAD_tools.h"

const size_t n_file=1;
const vector<string> filelist = yad::readFileList(n_file,"list/jeremy.list");

void PrtDep_Pt() {
    clock_t start_time=clock();

    /* graph declarations */
    
    TH1D* hPrtDist = new TH1D("hPrtDist","Distance between two PrtPt;distance (cm);#",100,0,1);
    double avg_prt_distance=0;
    TH1D* hNDepPt = new TH1D("hNDepPt","ratio NDep/NPt;ratio;#",100,0,100);

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

                size_t nppt = T.PrtNPt->at(i_prt);
                size_t ndep = T.PrtNDep->at(i_prt);

                hNDepPt->Fill((double) ndep/nppt);

                double X0 = (*T.PrtX)[i_prt][0];
                double Y0 = (*T.PrtY)[i_prt][0];
                double Z0 = (*T.PrtZ)[i_prt][0];
                for (size_t i_ppt=1; i_ppt < nppt; i_ppt++) {
                    double X = (*T.PrtX)[i_prt][i_ppt];
                    double Y = (*T.PrtY)[i_prt][i_ppt];
                    double Z = (*T.PrtZ)[i_prt][i_ppt];

                    hPrtDist->Fill(yad::distance(X0,Y0,Z0,X,Y,Z));
                    avg_prt_distance+=yad::distance(X0,Y0,Z0,X,Y,Z);

                    X0=X; Y0=Y; Z0=Z;
                } //end particlepoint loop
                if (nppt>1) avg_prt_distance/=(nppt-1);


                // for (size_t i_dep=0; i_dep < T.PrtNDep->at(i_prt); i_dep++) {

                // } //end deposit loop
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
    hNDepPt->Draw("hist");

    cout << "total time of execution: " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}