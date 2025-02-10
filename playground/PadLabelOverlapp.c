void PadLabelOverlapp(void)
{
   TCanvas* c = new TCanvas("c","c",500,500);
   Double_t eps = 0.005;
   TPad* p1 = new TPad("p1","p1",0.1,0.5,0.9,0.9,0); p1->Draw();
   TPad* p2 = new TPad("p1","p1",0.1,0.1,0.9,0.5+eps,0); p2->Draw();
   p1->SetBottomMargin(0);
   p2->SetTopMargin(0);   

   TH1F* h1 = new TH1F("h1","",100,-2.5,2.5);  
   TH1F* h2 = new TH1F("h2","",100,-2.5,2.5);  
   h1->Fill(0.,12000.);
   h2->Fill(0.,1.5);   
   p1->cd(); h1->Draw(); 
   p2->cd(); h2->Draw();
   c->cd();
}
