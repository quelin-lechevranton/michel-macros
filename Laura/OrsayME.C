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
#include "TH3F.h"
#include "TH2F.h"
#include "TLine.h"

using namespace std;

bool b_Simulation = true;

void GetCharacteristics(string in_root_file="PDVD_cosmics_LauraP_dumped.root", string out_root_file="no.root"){

    
    bool *fPFPIsTrack = nullptr;
    vector<int> *fPFPNClusters = nullptr;
    vector<int> *fPFPCluSummedADC = nullptr;

    TBranch *bvpx = nullptr;

    unsigned int fNPFParticles = 0;
    
    TFile file1("PDVD_cosmics_LauraP_dumped.root");
    TTree* tree_reco = (TTree*) file1.Get("LauraPDumper/Reco");
    Int_t nentries = tree_reco->GetEntries();
    cout<<nentries<<endl;

    for (Int_t iev=0; iev < tree_reco->GetEntries(); ++iev) { //Loop over the events
        tree_reco->GetEntry(iev);
        cout<<"Treating event number: "<<iev<<endl;
    }
    file1.Close();
    Int_t selectedpixel = 1;
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
  GetCharacteristics("PDVD_cosmics_LauraP_dumped.root", "st2_dis.root");

}
