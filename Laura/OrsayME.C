/*
Usage:
Authors:
- Thibaut Houdy, 4/11/22
*/

#include <vector>
#include <cmath>
#include <iostream>
#include "TGraph.h"
#include <string>
#include "TVector3.h"
#include "TPolyMarker3D.h"
#include "TPolyMarker.h"
#include "TPolyLine.h"
#include "TPolyLine3D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"

#include "TTree.h"


// #include "TH3F.h"
// #include "TH2F.h"
// #include "TLine.h"
// #include "TH1F.h"
// #include "TFile.h"
// #include "TStyle.h"
// #include "TTree.h"
// #include <iostream>
// #include <fstream>
// #include "TH1D.h"
// #include "TH2D.h"
// #include "TCanvas.h"
// #include <string>
// #include <map>
// #include "TF1.h"
// #include "TGraph2D.h"
// #include <TPolyMarker.h>
// #include <TAttMarker.h>
// #include "TGraph.h"
// #include "TGraphErrors.h"
// #include "TMath.h"
// #include "TGaxis.h"
// #include "TPRegexp.h"
// #include <iomanip>
// 
// #include <sstream>

using namespace std;

bool b_Simulation = true;

void GetCharacteristics(string in_root_file="PDVD_cosmics_LauraP_dumped.root", string out_root_file="no.root"){

    vector<double> fPFPTrackStartDirectionX,fPFPTrackStartDirectionY,fPFPTrackStartDirectionZ,fPFPTrackVertexDirectionX,fPFPTrackVertexDirectionY,fPFPTrackVertexDirectionZ;
    vector<int> *fPFPNClusters =0;
    vector<vector<double>> *fPFPCluPlane =0    ;
    vector<vector<double>> *fPFPCluView  =0    ;
    vector<vector<double>> *fPFPCluNHits  =0   ;
    vector<vector<double>> *fPFPCluSummedADC  =0;
    vector<vector<double>> *fPFPCluIntegral  =0;
    vector<vector<double>> *fPFPCluWidth  =0;
    
    TBranch *bvpx = nullptr;

    unsigned int fNPFParticles = 0;
    
    TFile file1(in_root_file);
    TTree* tree_reco = (TTree*) file1.Get("LauraPDumper/Reco");    
    
    tree_reco->SetBranchAddress("pfpTrackStartDirectionX",&fPFPTrackStartDirectionX);
    tree_reco->SetBranchAddress("pfpNClusters",             &fPFPNClusters             );
    tree_reco->SetBranchAddress("pfpCluNHits",              &fPFPCluNHits              );
    tree_reco->SetBranchAddress("pfpCluPlane",    &fPFPCluPlane);
    tree_reco->SetBranchAddress("pfpCluView",     &fPFPCluView);
    tree_reco->SetBranchAddress("pfpCluIntegral", &fPFPCluIntegral);
    tree_reco->SetBranchAddress("pfpCluSummedADC", &fPFPCluSummedADC);
    tree_reco->SetBranchAddress("nPFParticles", &fNPFParticles);
   // tree_reco->SetBranchAddress("pfpCluWidth", &fPFPCluWidth);

    Int_t nentries = tree_reco->GetEntries();
    cout<<nentries<<endl;

    for (Int_t iev=0; iev < tree_reco->GetEntries(); ++iev) { //Loop over the events
        cout<<"starting " <<iev<<endl;
        tree_reco->GetEntry(iev);
//         fPFPCluPlane[NumPart][Cluster] = [0,1,2]
        Int_t nclust = fPFPNClusters->size();
        Int_t ntest = fPFPCluPlane->size();
        Int_t ntest2 = fPFPCluPlane->at(0).size();
        cout<<"Treating event number: "<<iev<<"  "<<fNPFParticles<<"  "<< nclust<<"  "<<ntest<<"  "<<ntest2<<endl;
//         for(Int_t ipart=0; ipart<fNPFParticles;ipart++){
//                 if(ipart>0) continue;
//                 try{
//                     cout<<ipart<<"  "<<(fPFPCluPlane.at(ipart)).at(0)<<endl;
//                   //  cout<<"  "<<(fPFPCluPlane.at(ipart)).size()<<endl;
//                 }
//                 catch(const std::exception& e){
//                      std::cout << e.what(); 
//                 }
                
//         }
        
    }
     file1.Close();
     return;
}
/*
void Drawing3D_HitsAndTracks(TCanvas* canvas_3D, vector<vector<float>> l_vertexSpacePoints){
  //Drawing 3D tracks

  TH3F *axes = new TH3F("axes", "3D deposits distribution; X; Y; Z", 1, -400, 400, 1, -400, 400, 1, -10, 400);
  axes->SetDirectory(0);
  canvas_3D->cd();
  axes->Draw();

  TPolyMarker3D* otherPoly = new TPolyMarker3D(l_vertexSpacePoints.size());
  int point=0;
  for (auto const& sp : l_vertexSpacePoints) {
    otherPoly->SetPoint(point, sp[0], sp[1], sp[2]);
    point++;
  }
  otherPoly->SetMarkerStyle(8);
  otherPoly->SetMarkerSize(1);    
  otherPoly->SetMarkerColor(kBlue);

  canvas_3D->cd();
  otherPoly->Draw();   
  return;
}*/


void OrsayME(){
  GetCharacteristics("/silver/DUNE/quelin-lechevranton/out/PDVD_cosmics_LauraP_dumped.root", "st2_dis.root");
  return;
}
