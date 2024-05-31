#include "YAD_tools.h"

//analysis paramenters
vector<vector<double>> detlim={{-320, 350}, {-317,317}, {20,280}}; //cm
double thres_trk_len=0.; //cm
size_t n_last_clu=20;
double thres_dEdx=4; //Mev/cm

const size_t n_file=3;
const vector<string> filelist = yad::readFileList(n_file,"list/ijclab.list");

void Clustorimetry() {
    clock_t start_time=clock();

    TH2D* hdEdx = new TH2D("hdEdx",";Residual Range (cm);dE/dx (MeV/cm)",150,0,300,50,0,6);

    TH1D* hTrkLen = new TH1D("hTrkLen",";Track Length (cm);#",200,0,800);
    TH1D* hLastdEdx = new TH1D("hLastdEdx",";dE/dx (Mev/cm);#",50,0,6);

    size_t N_trk=0;
    size_t N_trk_sel=0;
    double avg_tpt=0;
    double avg_clu=0;
    double avg_len=0;

    size_t i_file=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        yad::Truth T(filename.c_str());
        yad::Reco R(filename.c_str());

        size_t n_evt = R.GetEntries();
        for (size_t i_evt=0; i_evt < n_evt; i_evt++) {

            cout << "Event#" << i_evt+1 << "/" << n_evt << "\r" << flush;

            T.GetEntry(i_evt);
            R.GetEntry(i_evt);

            N_trk+=R.NTrk;

            // for (size_t i_trk=0; i_trk < R.NTrk; i_trk++) {
            for (size_t i_pfp=0; i_pfp < R.NPfp; i_pfp++) {

                if (R.PfpTrkID->at(i_pfp) == -1) continue;
                size_t i_trk = R.PfpTrkID->at(i_pfp);

                if(!yad::isInside(
                    R.SptX->at(i_pfp),
                    R.SptY->at(i_pfp),
                    R.SptZ->at(i_pfp),
                    detlim[0][0], detlim[0][1],
                    detlim[1][0], detlim[1][1],
                    detlim[2][0], detlim[2][1]
                )) continue;

                size_t n_collection_clu=0;
                cout << "\tPfpNClu: " << R.PfpNClu->at(i_pfp) << endl;
                for (size_t i_clu=0; i_clu < R.PfpNClu->at(i_pfp); i_clu++) {
                    cout << "\tCluPlane: " << (*R.CluPlane)[i_pfp][i_clu] << endl;
                    if ((*R.CluPlane)[i_pfp][i_clu]==0) n_collection_clu++;
                } //end clu loop
                cout << "nColClu: " << n_collection_clu << " vs. nPfpSpt: " << R.PfpNSpt->at(i_pfp) << endl;

                // hTrkLen->Fill(R.TrkLength->at(i_trk));
                // if(R.TrkLength->at(i_trk) < thres_trk_len) continue;

                // double avg_dEdx=0;
                
                // for (size_t i_clu=R.PfpNClu->at(i_pfp)-n_last_clu; i_clu < R.PfpNClu->at(i_pfp); i_clu++) {
                //     avg_dEdx+=(*R.TrkCaldEdx)[i_trk][i_clu]/n_last_clu;
                // } //end cal loop
                // hLastdEdx->Fill(avg_dEdx);
                // if (avg_dEdx < thres_dEdx) continue;

                // N_trk_sel++;
                // avg_clu+=R.PfpNClu->at(i_pfp);
                // avg_len+=R.TrkLength->at(i_trk);

                // for (size_t i_clu=0; i_clu < R.TrkCalNPt->at(i_trk); i_clu++) {
                //     hdEdx->Fill((*R.TrkCalResRange)[i_trk][i_clu],(*R.TrkCaldEdx)[i_trk][i_clu]);
                // } //end cal loop
            } //end trk loop
        } //end event loop
    } //end file loop
    cout << endl;

    // TCanvas* c1 = new TCanvas("c1","Clustorimetry");
    // c1->cd();
    // c1->Divide(2,2);
    // c1->cd(1);
    // gPad->SetLogz();
    // hdEdx->Draw("colZ");
    // c1->cd(2);
    // gPad->SetLogy();
    // hTrkLen->Draw("hist");
    // c1->cd(3);
    // gPad->SetLogy();
    // hLastdEdx->Draw("hist");
    
    // avg_clu /= N_trk_sel;
    // avg_len /= N_trk_sel;
    // cout << N_trk <<  " tracks" << endl;
    // cout << N_trk_sel << " tracks selected" << endl;
    // cout << avg_clu << " calorimetry point per track on average" << endl;
    // cout << avg_len << " cm track length on average" << endl;
    cout << "total time of execution: " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}
