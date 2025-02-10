void BinLabels()
{
   const Int_t nx = 12;
   const Int_t ny = 20;
   const char *months[nx] = {"January","February","March","April","May","June","July",
                      "August","September","October","November","December"};
   const char *people[ny] = {"Jean","Pierre","Marie","Odile","Sebastien","Fons","Rene",
      "Nicolas","Xavier","Greg","Bjarne","Anton","Otto","Eddy","Peter","Pasha",
      "Philippe","Suzanne","Jeff","Valery"};
   TCanvas *c1 = new TCanvas("c1","demo bin labels",10,10,800,800);
   c1->SetGrid();
   c1->SetLeftMargin(0.15);
   c1->SetBottomMargin(0.15);
   TH2F *h = new TH2F("h","test",nx,0,nx,ny,0,ny);
   Int_t i;
   for (i=0;i<5000;i++) {
      h->Fill(gRandom->Gaus(0.5*nx,0.2*nx), gRandom->Gaus(0.5*ny,0.2*ny));
   }
   h->SetStats(0);
   h->GetXaxis()->SetLabelOffset(99);
   h->GetYaxis()->SetLabelOffset(99);
   h->Draw("text");
   // draw labels along X
   Float_t x, y;
   y = gPad->GetUymin() - 0.2*h->GetYaxis()->GetBinWidth(1);
   TText t;
   t.SetTextAngle(60);
   t.SetTextSize(0.02);
   t.SetTextAlign(33);
   for (i=0;i<nx;i++) {
      x = h->GetXaxis()->GetBinCenter(i+1);
      t.DrawText(x,y,months[i]);
   }
   // draw labels along y
   x = gPad->GetUxmin() - 0.1*h->GetXaxis()->GetBinWidth(1);
   t.SetTextAlign(32);
   t.SetTextAngle(0);
   for (i=0;i<ny;i++) {
      y = h->GetYaxis()->GetBinCenter(i+1);
      t.DrawText(x,y,people[i]);
   }
}
