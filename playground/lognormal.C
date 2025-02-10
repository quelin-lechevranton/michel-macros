


void lognormal(){
  
// TF1* m_logNormal = new TF1("logNormal","[2]/(TMath::Log([1])*2.5066)*TMath::Exp(-(TMath::Log(x)-TMath::Log([0]))*(TMath::Log(x)-TMath::Log([0]))/(2*TMath::Log([1])*TMath::Log([1])))/x", 0, 10000);
  TF1* m_logNormal = new TF1("logNormal","[2]*ROOT::Math::lognormal_pdf(x, log([0]), log([1]) )",0,1000);

  m_logNormal->SetParameter(2, 70);
  m_logNormal->SetParameter(0, 15);
  m_logNormal->SetParameter(1, 3.9);
  m_logNormal->SetNpx(999999);

  //copy as the first TF1 will be fitted.
//  TF1* m_logNormal2 = new TF1("logNormal2","[2]/(TMath::Log([1])*2.5066)*TMath::Exp(-(TMath::Log(x)-TMath::Log([0]))*(TMath::Log(x)-TMath::Log([0]))/(2*TMath::Log([1])*TMath::Log([1])))/x", 0, 10000);
  TF1* m_logNormal2 = new TF1("logNormal2","[2]*ROOT::Math::lognormal_pdf(x, log([0]), log([1]) )",0,1000);
  m_logNormal2->SetParameter(2, 70);
  m_logNormal2->SetParameter(0, 15);
  m_logNormal2->SetParameter(1, 3.9);
  m_logNormal2->SetNpx(999999);

  const int nbins=200;
  Double_t xbins[nbins+1];
  double dx = 4./nbins;
  double l10 = TMath::Log(10);
  for (int i=0;i<=nbins;i++) {
     xbins[i] = 0.1*TMath::Exp(l10*i*dx);
  }
  TH1I* trueclusterIntensities = new TH1I("trueclusterintensities", "", nbins, xbins);
  for(int i=0; i < 165000; i++){
    trueclusterIntensities->Fill(m_logNormal->GetRandom());
  }

  trueclusterIntensities->Sumw2();
  trueclusterIntensities->Scale(1,"width");
  trueclusterIntensities->Draw();

 TCanvas* c1=  new TCanvas();
  c1->cd();
  c1->SetLogx();
  trueclusterIntensities->Fit(m_logNormal,"");
  //draw TF1 from which random values are generated.

  m_logNormal2->SetLineColor(kBlue);
  m_logNormal2->SetParameter(2, m_logNormal->GetParameter(2));
  m_logNormal2->Draw("same");

  
}
