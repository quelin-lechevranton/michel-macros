#include "YAD_tools.h"

/*
 * muon Bragg peak analysis preliminary for TrueMichel_v2_1.C
*/

vector<vector<double>> detlim={{-350, 350}, {-350,350}, {0,300}}; //cm
// vector<vector<double>> detlim={{-320, 350}, {-317,317}, {20,280}}; //cm

const size_t n_bragg_integration=50;
const size_t n_bragg_tail=1;

const size_t n_file=1;
const vector<string> filelist = yad::readFileList(n_file,"list/ijclab.list");

void BraggHit(size_t i=0) {
    clock_t start_time=clock();

    // TH1D* hAvgdEdx = new TH1D ("hAvgdEdx","??Average dEdx on last deposits;MeV/cm;#",30,0,30);

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
            

            for (size_t i_pfp=0; i_pfp < R.NPfp; i_pfp++) {

                if (R.PfpTrkID->at(i_pfp)<0) continue;
                size_t i_trk = R.PfpTrkID->at(i_pfp);

                size_t n_spt = R.PfpNSpt->at(i_pfp);

                cout << "\t\tpfp#" << i_pfp << ":#spt" << n_spt << endl;

                // if (!yad::isInside(
                //     R.SptX->at(i_pfp),
                //     R.SptY->at(i_pfp),
                //     R.SptZ->at(i_pfp),
                //     detlim[0][0], detlim[0][1],
                //     detlim[1][0], detlim[1][1],
                //     detlim[2][0], detlim[2][1]
                // )) continue;
                // if (n_spt<n_bragg_integration+n_bragg_tail) continue;

                // // double avg_dEdx=0;
                // double Range=0;
                // double X0=(*R.SptX)[i_pfp][0];
                // double Y0=(*R.SptY)[i_pfp][0];
                // double Z0=(*R.SptZ)[i_pfp][0];
                // for (size_t i_spt=1; i_spt < n_spt; i_spt++) {
                //     double X=(*R.SptX)[i_pfp][i_spt];
                //     double Y=(*R.SptY)[i_pfp][i_spt];
                //     double Z=(*R.SptZ)[i_pfp][i_spt];

                //     double localRange=yad::distance(X0,Y0,Z0,X,Y,Z);
                //     Range += localRange;
                //     X0=X; Y0=Y; Z0=Z;

                //     // for (size_t i_hit=0; i_hit < n_hit; i_hit++) {

                //     //     if ((*R.HitSptID)[i_pfp][i_hit] != i_spt) continue;
                //     //     if ((*R.HitPlane)[i_pfp][i_hit] != 2) continue;
                //     // }
                // }

                // double TrkLength=R.TrkLength->at(i_trk);
                // cout << "\t\t\tTrkLen: " << TrkLength << " vs. SptRange: " << Range << endl;
                // if (TMath::Abs(TrkLength-Range)/TrkLength > .1) {
                //     cout << "\t\t\t\e[91mstrange spt\e[0m" << endl;
                //     continue;
                // }

                // double avg_nHit_perSpt=0;
                // for (size_t i_spt=0; i_spt < n_spt; i_spt++) {

                //     double nHit=0;
                //     for (size_t i_hit=0; i_hit < n_hit; i_hit++) {

                //         if ((*R.HitSptID)[i_pfp][i_hit] != i_spt) continue;
                //         if ((*R.Hit))
                //         nHit++;

                //     }
                //     avg_nHit_perSpt+=nHit;
                // }
                // avg_nHit_perSpt/=n_spt;
                // cout << "\t\t\tavg_nHit_perSpt: " << avg_nHit_perSpt << endl;

                // hAvgdEdx->Fill(avg_dEdx);
            } //end particle loop
        } //end event loop
    } //end file loop
    // cout << endl;

    // TCanvas* c1 = new TCanvas("c1","BraggHit");
    // c1->cd();
    // hAvgdEdx->Draw("hist");

    cout << "total time of execution: " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}

// void PlotBragg(size_t i_file=1, size_t i_evt=1) {
//     clock_t start_time=clock();

//     TGraph* gSptX = new TGraph();
//     gSptX->SetTitle(";i_spt;SptX (cm)");    
//     TGraph* gSptY = new TGraph();
//     gSptY->SetTitle(";i_spt;SptY (cm)");    
//     TGraph* gSptZ = new TGraph();
//     gSptZ->SetTitle(";i_spt;SptZ (cm)");    

//     TGraph* gSumADC = new TGraph();
//     gSumADC->SetTitle(";Range (cm);SumADC");

//     TGraph* gLocalRange = new TGraph();
//     gLocalRange->SetTitle(";Range (cm);Local Range (cm)");

//     TGraph* gdSdx = new TGraph();
//     gdSdx->SetTitle(";Range (cm);SumADC/Local Range (/cm)");

//     string filename = filelist[i_file-1];

//     // yad::Truth T(filename.c_str());
//     yad::Reco R(filename.c_str());
    
//     // T.GetEntry(i_evt-1);
//     R.GetEntry(i_evt-1);

//     for (size_t i_pfp=0; i_pfp < R.NPfp; i_pfp++) {

//         if (R.PfpTrkID->at(i_pfp)<0) continue;
//         size_t i_trk = R.PfpTrkID->at(i_pfp);

//         size_t n_spt = R.PfpNSpt->at(i_pfp);
//         size_t n_hit = R.PfpNHit->at(i_pfp);

//         // if (!yad::isInside(
//         //     R.SptX->at(i_pfp),
//         //     R.SptY->at(i_pfp),
//         //     R.SptZ->at(i_pfp),
//         //     detlim[0][0], detlim[0][1],
//         //     detlim[1][0], detlim[1][1],
//         //     detlim[2][0], detlim[2][1]
//         // )) continue;

//         double Range=0;
//         double X0=(*R.SptX)[i_pfp][0];
//         double Y0=(*R.SptY)[i_pfp][0];
//         double Z0=(*R.SptZ)[i_pfp][0];
//         for (size_t i_spt=1; i_spt < n_spt; i_spt++) {
//             double X=(*R.SptX)[i_pfp][i_spt];
//             double Y=(*R.SptY)[i_pfp][i_spt];
//             double Z=(*R.SptZ)[i_pfp][i_spt];

//             double localRange=yad::distance(X0,Y0,Z0,X,Y,Z);
//             // if (localRange < 0.5) continue;
//             Range += localRange;
//             X0=X; Y0=Y; Z0=Z;

//             gSptX->AddPoint(i_spt,X);
//             gSptY->AddPoint(i_spt,Y);
//             gSptZ->AddPoint(i_spt,Z);

//             for (size_t i_hit=0; i_hit < n_hit; i_hit++) {

//                 if ((*R.HitSptID)[i_pfp][i_hit] != i_spt) continue;
//                 if ((*R.HitPlane)[i_pfp][i_hit] != 2) continue;

//                 double SumADC=(*R.HitSumADC)[i_pfp][i_hit];
//                 gSumADC->AddPoint(Range,SumADC);
//                 gLocalRange->AddPoint(Range,localRange);
//                 gdSdx->AddPoint(Range,SumADC/localRange);
//             }
//         }

//     } //end particle loop
//     cout << endl;

//     TCanvas* c2 = new TCanvas("c2","BraggHit");
//     c2->Divide(2,2);
//     c2->cd(1);
//     gSptX->Draw();
//     c2->cd(2);
//     gSptY->Draw();
//     c2->cd(3);
//     gSptZ->Draw();

//     TCanvas* c3 = new TCanvas("c3","BraggHit");
//     c3->Divide(2,2);
//     c3->cd(1);
//     gSumADC->Draw();
//     c3->cd(2);
//     gLocalRange->Draw();
//     c3->cd(3);
//     gdSdx->Draw();

//     cout << "total time of execution: " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
// }