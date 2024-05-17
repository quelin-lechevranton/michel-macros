#include <vector>
#include <string>
#include <iostream>
#include "TH1D.h"

using namespace std;

const vector<string> filelist = {
    "/silver/DUNE/quelin-lechevranton/out/protodunevd_10_muon_500MeV_Jdumped.root",
    "/silver/DUNE/quelin-lechevranton/out/protodunevd_100_muon_800MeV_Jdumped.root",
    "/silver/DUNE/quelin-lechevranton/out/pdvd_1k_muon_2GeV_allangles_Jdumped.root",
    "/silver/DUNE/quelin-lechevranton/out/pdvd_1k_muon_1500MeV_Jdumped.root",
    "/silver/DUNE/quelin-lechevranton/out/pdvd_100_muon_1GeV_Jdumped.root"
};

void Clusters() {

    int N_evt=0;

    int n_clu=0;
    int n_pfp=0;

    vector<int> *PfpPdg=nullptr,
                *PfpID=nullptr,
                *CluPfpID=nullptr;

    vector<double>  *CluPlane=nullptr,
                    *CluNHits=nullptr,
                    *CluSum=nullptr,
                    *CluInt=nullptr,
                    *CluWidth=nullptr;    

    // Hist ==================================================
    vector<TH1D*> h(2);
    int n_bin=30, x_min=0, x_max=2000;

    h[0] = new TH1D("hSum/Width","",n_bin,x_min,x_max);
    h[0]->SetLineColor(kRed+1);
    h[0]->SetLineWidth(2);

    h[1] = new TH1D("hInt/Width","",n_bin,x_min,x_max);
    h[1]->SetLineColor(kBlue-3);
    h[1]->SetLineWidth(2);

    for (int i_file=0; i_file < filelist.size(); i_file++) {
        string filename = filelist[i_file];

        cout << "\e[3mOpening file #" << i_file+1 << ": " << filename << "\e[0m" << endl;

        TFile file(filename.c_str());
        TTree *Reco = file.Get<TTree>("JDumper/Reco");
        // TTree *Truth = file.Get<TTree>("JDumper/Truth");
        int n_evt=Reco->GetEntries();
        N_evt+=n_evt;

        Reco->SetBranchAddress("fPfpID",            &PfpID);
        Reco->SetBranchAddress("fNClusters",        &n_clu);
        Reco->SetBranchAddress("fNPFParticles",     &n_pfp);
        Reco->SetBranchAddress("fPfpPdg",           &PfpPdg);

        Reco->SetBranchAddress("fCluPfpID",         &CluPfpID);
        Reco->SetBranchAddress("fCluPlaneID",       &CluPlane);
        Reco->SetBranchAddress("fCluNHits",         &CluNHits);
        Reco->SetBranchAddress("fCluSumADC",        &CluSum);
        Reco->SetBranchAddress("fCluIntFit",        &CluInt);
        Reco->SetBranchAddress("fCluWidth",         &CluWidth);

        // Truth->SetBranchAddress("", &);

        double Sum,Int,Width;

        for (int i_evt=0; i_evt < n_evt; i_evt++) {

            Reco->GetEntry(i_evt);
            // Truth->GetEntry(i_evt);

            for(int i_pfp=0; i_pfp < n_pfp ; i_pfp++) {

                if(PfpPdg->at(i_pfp)!=13) continue;

                for(int i_clu=0; i_clu < n_clu; i_clu++) {
                    
                    // if(CluPlane->at(i_pfp)[i_clu]!=0) continue;

                    if (CluPfpID->at(i_clu)!=PfpID->at(i_pfp)) {continue;}

                    Sum=CluSum->at(i_clu);
                    Int=CluInt->at(i_clu);
                    Width=CluWidth->at(i_clu);

                    h[0]->Fill(Sum/Width);
                    h[1]->Fill(Int/Width);
                }

            }
        
        } //end of event loop

        file.Close();
    
    } //end of file loop

    auto canvas = new TCanvas("c1","muon dQ/dx on collection");
    canvas->cd();
    h[1]->Draw("hist");
    h[0]->Draw("samehist");
    // Reco->Draw("pfpCluSummedADC/pfpCluWidth","pfpPdgCode==13 && pfpCluPlane==0");
    // canvas->SaveAs("Cluster.root");
    // canvas->SaveAs("Cluster.pdf");
}