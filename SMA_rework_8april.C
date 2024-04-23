/*
  G4 criteria selection for Michel electrons
  13/07/23 -> Camilo Andres-Cortes
  03/01/24 -> Thibaut Houdy
  08/04/24 -> Jeremy Quelin Lechevranton
  Particle PDGcode:
  https://pdg.lbl.gov/2019/reviews/rpp2018-rev-monte-carlo-numbering.pdf
*/




R__ADD_INCLUDE_PATH("gallery/Event.h")
R__ADD_INCLUDE_PATH("lardataobj/Simulation/SimEnergyDeposit.h")
R__ADD_INCLUDE_PATH("larcoreobj/SimpleTypesAndConstants/geo_types.h")
R__ADD_INCLUDE_PATH("larcoreobj/SimpleTypesAndConstants/geo_vectors.h")
R__ADD_INCLUDE_PATH("nusimdata/v1_27_01/include/nusimdata/SimulationBase")
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "tools.h"

/************************************ MuonAnalysis ************************************/

vector<pair<TVector3, int>> MuonAnalysis(int n_files, int i_first_event, int i_last_event, string tag_energy_deposit, string tag_MCTruth, string listname, bool is_muon)
{
  //Particle options
  vector<string> pdgnames = {"e^{-}", "e^{+}", "#mu^{-}", "#mu^{+}", "p^{+}"};
  vector<int> pdgcodes    = {11, -11, 13, -13, 2212};
  Color_t colors[6]       = {kBlue+2, kAzure+8, kOrange+7, kRed+2, kGreen+2, kMagenta+2};
  
  //To distinguish between muon/antimuon generated
  int pdg;
  string particle_name;
  string generated;
  if(is_muon) {pdg = pdgcodes[2]; particle_name = pdgnames[2]; generated = "muon";}
  if(!is_muon){pdg = pdgcodes[3]; particle_name = pdgnames[3]; generated = "antimuon";}

  vector<string> ytitle_muon = {("Deposited "+particle_name+" Energy (MeV)").c_str(), ("Remaining "+particle_name+" Energy (MeV)").c_str()};
  vector<TGraph*> muon_graph(2);

  for(int i = 0; i < muon_graph.size(); i++){
    muon_graph[i] = new TGraph();
    muon_graph[i]->SetName(Form("muon_graph_%d", i));
    muon_graph[i]->SetMarkerColor(colors[i+3]);
    muon_graph[i]->SetLineColor(colors[i+3]);
    muon_graph[i]->SetMarkerStyle(5);
    muon_graph[i]->SetMarkerSize(0.5);
    muon_graph[i]->GetXaxis()->SetTitle("Track Length (cm)");
    muon_graph[i]->GetYaxis()->SetTitle(ytitle_muon[i].c_str());
  }

  //All tags that you need
  art::InputTag InputTag_energy_deposit(tag_energy_deposit);
  art::InputTag InputTag_MCTruth(tag_MCTruth);  
  //Search your files
  string file_list          = listname+".list";
  vector<string> filenames = ReadFileList(n_files, file_list);

  int i_muon_total_depo = 0;    //to account for all points from all events
  int i_event = 0;
  vector<pair<TVector3, int>> muon_analysis(i_last_event-(i_first_event-1));

  //Go for events!
  for(gallery::Event event(filenames); !event.atEnd(); event.next()){
    //Physical quantities
    vector<double> array_muon_residual_track;
    vector<double> array_muon_track_length;       //hit position from start hit
    vector<double> array_muon_depo_energy;     //energy deposition in each hit
    vector<double> array_muon_remaining_energy;      //remaining energy in each hit

    if(event.eventAuxiliary().event() < i_first_event) continue;
    if(event.eventAuxiliary().event() > i_last_event) break; 

    cout << "Event: " << event.eventAuxiliary().event() << endl;
      
    //Retrieve list of energy deposits per event
    auto const depo_list = event.getValidHandle<vector<sim::SimEnergyDeposit>>(InputTag_energy_deposit);
    size_t n_depo        = depo_list->size();

    //Retrieve generated information   
    auto const montecarlo        = event.getValidHandle<vector<simb::MCTruth>>(InputTag_MCTruth);
    const simb::MCTruth& gen     = montecarlo->at(0);
    const simb::MCParticle& muon = gen.GetParticle(0);

    double muon_remaining_energy = (1000)*muon.E();   //in MeV
    cout << "Muon/AntiMuon Energy: " << muon_remaining_energy << " MeV" << endl;
    double muon_depo_energy_tot   = 0;
    int    i_muon_last_hit = -1;

    geo::Point_t muon_first_hit = depo_list->at(0).MidPoint();
    geo::Point_t muon_previous_hit = muon_first_hit;
    double muon_track_length = 0.;

    //Depo loop for muon/antimuon information
    for(size_t i_depo = 0; i_depo < n_depo; i_depo++){

      const sim::SimEnergyDeposit& depo = depo_list->at(i_depo);
      
      //Retrieve position of energy deposit
      geo::Point_t hit = depo.MidPoint();
      //####################################################################
      if(depo.PdgCode() == pdg) {

        muon_track_length += sqrt( pow(muon_previous_hit.X() - hit.X(), 2) + pow(muon_previous_hit.Y() - hit.Y(), 2) + pow(muon_previous_hit.Z() - hit.Z(), 2) );
    
        muon_remaining_energy -= depo.Energy();
        muon_depo_energy_tot  += depo.Energy();
        
        //Energy deposit and position for muon/antimuon
        array_muon_track_length.push_back(abs(muon_track_length));
        array_muon_depo_energy.push_back(depo.Energy());
        array_muon_remaining_energy.push_back(muon_remaining_energy);

        i_muon_last_hit = i_depo;              //save last hit
        muon_previous_hit = hit;
      } //end of muon/anitmuon condition
    }   //end of depo loop 
    
    if(array_muon_track_length.size() == 0){cout << "There is no muon/antimuon information here!" << endl; i_event++; continue;}    //go back to begining if you don't have muon information
    
    geo::Point_t muon_last_hit = depo_list->at(i_muon_last_hit).MidPoint();

    cout << "Muon/AntiMuon Total Deposit Energy: " << muon_depo_energy_tot << " MeV" << endl;
    cout << "Muon/AntiMuon Last Hit: " << muon_last_hit << endl;


    //Direction vector  
    TVector3 muon_track_global_vector;
    muon_track_global_vector(0) = muon_first_hit.X()-muon_last_hit.X();
    muon_track_global_vector(1) = muon_first_hit.Y()-muon_last_hit.Y();
    muon_track_global_vector(2) = muon_first_hit.Z()-muon_last_hit.Z();
    muon_track_global_vector    = muon_track_global_vector.Unit();


    int n_track_length = array_muon_track_length.size();
    for (int i=0 ; i<n_track_length ; i++ ){
      array_muon_residual_track.push_back(array_muon_track_length[n_track_length-i-1]);
    }
  
    for (int i = 0; i < array_muon_track_length.size(); i++){   
        muon_graph[0]->SetPoint(i_muon_total_depo, array_muon_residual_track[i], array_muon_depo_energy[i]);
        muon_graph[1]->SetPoint(i_muon_total_depo, array_muon_track_length[i], array_muon_remaining_energy[i]);
        i_muon_total_depo++;
    }   //end of muon/antimuon information loop

    //pair<vector<TVector3>, geo::Point_t> p(muon_track_global_vector, muon_last_hit);
    pair<TVector3, int> p(muon_track_global_vector, i_muon_last_hit);
    muon_analysis[i_event] = p;
    i_event++;
  } //end of event loop

  printf("####################################");

  //Plot muon/antimuon information
  TCanvas* canvas_muon_depo_energy = new TCanvas("canvas_muon_depo_energy", "Muon/AntiMuon Energy Deposit");
  canvas_muon_depo_energy->cd();
  muon_graph[0]->Draw("AP");
  canvas_muon_depo_energy->Update();

  TCanvas* canvas_muon_remaining_energy = new TCanvas("canvas_muon_remaining_energy", "Muon/AntiMuon Remaining Energy");
  canvas_muon_remaining_energy->cd();
  muon_graph[1]->Draw("AP");
  canvas_muon_remaining_energy->Update();

  //Save it!
  canvas_muon_depo_energy->SaveAs(("output/"+generated+"_depo_energy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".pdf").c_str()); 
  canvas_muon_depo_energy->SaveAs(("output/"+generated+"_depo_energy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".root").c_str());

  canvas_muon_remaining_energy->SaveAs(("output/"+generated+"_remaining_energy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".pdf").c_str()); 
  canvas_muon_remaining_energy->SaveAs(("output/"+generated+"_remaining_energy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".root").c_str());

  return muon_analysis;
} // end of MuonAnalysis

/************************************ ElectronAnalysis ************************************/


void ElectronAnalysis(int n_files, int i_first_event, int i_last_event, string tag_energy_deposit, string tag_MCTruth, string tag_MCParticle, string listname, vector<pair<TVector3, int>> muon_analysis, bool muon)
{
  //X-axis
  vector<string> xtitle = {"#theta (degrees)", "Energy (MeV)", "Energy (MeV)", "Energy (MeV)", "Energy (MeV)"};
  vector<int> bins      = {180, 90, 90, 90, 90};
  vector<float> xmin    = {-10, 0, 0, 0, 0};
  vector<float> xmax    = {200, 360, 100, 60, 60};
  //Y-axis
  vector<string> ytitle_elec = {"Total Deposited e^{-}/e^{+} Energy (MeV)", "Completeness"};
  //Particle options
  vector<string> pdgnames = {"e^{-}", "e^{+}", "#mu^{-}", "#mu^{+}", "p^{+}"};
  vector<int> pdgcodes    = {11, -11, 13, -13, 2212};
  Color_t colors[7]       = {kBlue+2, kAzure+8, kOrange+7, kRed+2, kGreen-3, kMagenta+2, kBlue-6};

  //To distinguish between muon/antimuon generated
  int pdg;
  string particle_name;
  string generated;
  if(muon) {pdg = pdgcodes[2]; particle_name = pdgnames[2]; generated = "muon";}
  if(!muon){pdg = pdgcodes[3]; particle_name = pdgnames[3]; generated = "antimuon";}

  vector<TH1D*> elec_histograms(5);
  for(int i = 0; i < elec_histograms.size(); i++){
    elec_histograms[i] = new TH1D(Form("elec_histograms_%d", i), Form(";%s;Counts", xtitle[i].c_str()), bins[i], xmin[i], xmax[i]);
    if(muon) {elec_histograms[i]->SetLineColor(colors[6]);  elec_histograms[i]->SetFillColor(colors[6]);}
    if(!muon){elec_histograms[i]->SetLineColor(colors[1]);  elec_histograms[i]->SetFillColor(colors[1]);} 
  }
  elec_histograms[0]->SetLineColor(colors[5]);  elec_histograms[0]->SetFillColor(colors[5]);

  vector<TGraph*> gene_sel(2);
  for(int i = 0; i < gene_sel.size(); i++){
    gene_sel[i] = new TGraph();
    gene_sel[i]->SetName(Form("gene_sel_%d", i));
    if(muon) {gene_sel[i]->SetMarkerColor(colors[6]);  gene_sel[i]->SetLineColor(colors[6]); gene_sel[i]->SetFillColor(colors[6]);}
    if(!muon){gene_sel[i]->SetMarkerColor(colors[1]);  gene_sel[i]->SetLineColor(colors[1]); gene_sel[i]->SetFillColor(colors[1]);} 
    gene_sel[i]->SetMarkerStyle(3);
    gene_sel[i]->SetMarkerSize(0.25);
    gene_sel[i]->GetXaxis()->SetTitle("#theta (degrees)");
    gene_sel[i]->GetYaxis()->SetTitle(ytitle_elec[i].c_str());
  }

  TH1D* hposi = new TH1D("", Form("Decay MC e^{+};%s;Counts", xtitle[4].c_str()), bins[4], xmin[4], xmax[4]);
  hposi->SetLineColor(colors[4]);  hposi->SetFillColor(colors[4]);

  TH1D* hdiff = new TH1D("", Form("E_{e^{+}, MC} - E_{e^{+}, Selected};%s;Counts", xtitle[4].c_str()), bins[4], -xmax[4], xmax[4]);
  hdiff->SetLineColor(colors[3]);  //hdiff->SetFillColor(colors[3]);

  //All tags that you need
  art::InputTag InputTag_energy_deposit(tag_energy_deposit);
  art::InputTag InputTag_MCTruth(tag_MCTruth);
  //art::InputTag InputTag_MCParticle(tag_MCParticle); 
  //Search your files
  string file_list          = listname+".list";
  vector<string> filenames = ReadFileList(n_files, file_list);

  vector<double> Posi_Energy(i_last_event-(i_first_event-1));           //to store each decay positron energy from mc
  const int n_bin = 180;
  vector<float> elec_energy_ratio_angular_distrubtion(n_bin);
  vector<float> completeness(n_bin);
  //int elec_pc = 0;                                          //to account for all points from all events
  int i_event    = 0;
  
  //Go for events!
  for(gallery::Event event(filenames); !event.atEnd(); event.next()){
    //Physical quantities
    // vector<double> elecDist;    //distance of electron hit from last muon hit
    // vector<double> elecDepoE;   //deposited energy of selected electrons

    if(event.eventAuxiliary().event() < i_first_event) continue;
    if(event.eventAuxiliary().event() > i_last_event) break; 
      
    //Retrieve list of energy deposits per event
    auto const depo_list = event.getValidHandle<vector<sim::SimEnergyDeposit>>(InputTag_energy_deposit);
    size_t n_depo       = depo_list->size();

    //Decay mc positron information ~ only for antimuon events
    /*
    if(!muon){
      auto const posicarlo = event.getValidHandle<vector<simb::MCParticle>>(InputTag_MCParticle);
      size_t nposi         = posicarlo->size();

      for(size_t pos_i = 0; pos_i < nposi; pos_i++){
        const simb::MCParticle& positron = posicarlo->at(pos_i);
        if(positron.PdgCode() == -11 && positron.Process() == "Decay"){Posi_Energy[i_event] = (1000)*positron.E(); hposi->Fill(Posi_Energy[i_event]);}
      }    
    } */ 

    if(muon_analysis[i_event].second == 0){cout << "There is no muon information here!" << endl; i_event++; continue;}    //go back to begining if you don't have muon information

    //Information from muon/anitmuon analysis     
    TVector3 muon_track_global_vector = muon_analysis[i_event].first;      
    int i_muon_last_hit = muon_analysis[i_event].second;

    geo::Point_t muon_last_hit = depo_list->at(i_muon_last_hit).MidPoint();

    //Electron information
    double distance_max  = 25.0;        //max. distance from last muon hit (in cm) ~ electron analysis
    double angle_barycentre_depo;               //cone opening angle
    double angle_track_depo;
    int angle_bin; 
    double elec_depo_energy_tot   = 0;
    double elec_depo_energy_mask  = 0; 
    double elec_depo_energy_sphere   = 0;
    double elec_depo_energy_cone = 0;       

    TVector3 depo_weighted_point;             //weighted electron hits
    TVector3 elec_barycenter;         //coordinates of barycenter shower
    TVector3 mu_to_barycenter_vector;            //axis cone vector
    //TVector3 r_ElecHit;         //all electron hits vector
    TVector3 mu_to_depo_vector;
    double mu_to_depo_distance;
    
    //Compute barycenter
    for(size_t i_depo = 0; i_depo < n_depo; i_depo++){

      const sim::SimEnergyDeposit& depo = depo_list->at(i_depo);
      
      if(abs(depo.PdgCode()) == 11){  

        //Retrieve position of energy deposit
        geo::Point_t depo_hit = depo.MidPoint();

        mu_to_depo_vector(0) = depo_hit.X()-muon_last_hit.X();
        mu_to_depo_vector(1) = depo_hit.Y()-muon_last_hit.Y();
        mu_to_depo_vector(2) = depo_hit.Z()-muon_last_hit.Z();

        //Distance and angle of electron hit from muon last hit
        mu_to_depo_distance = mu_to_depo_vector.Mag();
        angle_track_depo  = muon_track_global_vector.Angle(mu_to_depo_vector)*TMath::RadToDeg(); 

        if(abs(angle_track_depo) >= 20.0){     //muon/antimuon track masking 

          if(mu_to_depo_distance <= distance_max){

            elec_depo_energy_sphere += depo.Energy();

            depo_weighted_point(0)  += depo.Energy()*depo_hit.X();
            depo_weighted_point(1)  += depo.Energy()*depo_hit.Y();
            depo_weighted_point(2)  += depo.Energy()*depo_hit.Z();

          } //end of sphere condition
        }   //end of track masking condition
      }     //end of electron condition
    }       //end of depo loop


    if(elec_depo_energy_sphere <= 0){cout << "There is no Michel electron information here!" << endl; i_event++; continue;}    //go back to begining if you don't have Michel electron information

    //Fill barycenter shower vector
    for(int i = 0; i < 3; i++){
      elec_barycenter(i) = depo_weighted_point(i)/elec_depo_energy_sphere;
    }

    cout << "Barycenter: (" << elec_barycenter.X() << "," << elec_barycenter.Y() << "," << elec_barycenter.Z() << ")" << endl;

    //Find normalized axis cone vector  
    mu_to_barycenter_vector(0) = elec_barycenter.X()-muon_last_hit.X();
    mu_to_barycenter_vector(1) = elec_barycenter.Y()-muon_last_hit.Y();
    mu_to_barycenter_vector(2) = elec_barycenter.Z()-muon_last_hit.Z();
    mu_to_barycenter_vector    = mu_to_barycenter_vector.Unit();

    vector<float> elec_energy_angular_distribution(n_bin);  

    //Depo loop for michel electron selection
    for(size_t i_depo = 0; i_depo < n_depo; i_depo++){
      
      const sim::SimEnergyDeposit& depo = depo_list->at(i_depo);
      
      if(abs(depo.PdgCode()) == 11){        

        elec_depo_energy_tot += depo.Energy();  

        //Retrieve position of energy deposit
        geo::Point_t depo_hit = depo.MidPoint();
        
        mu_to_depo_vector(0) = depo_hit.X()-muon_last_hit.X();
        mu_to_depo_vector(1) = depo_hit.Y()-muon_last_hit.Y();
        mu_to_depo_vector(2) = depo_hit.Z()-muon_last_hit.Z();


        //Distance and angle of electron hit from muon last hit
        mu_to_depo_distance = mu_to_depo_vector.Mag();
        angle_track_depo  = muon_track_global_vector.Angle(mu_to_depo_vector)*TMath::RadToDeg(); 


        if(abs(angle_track_depo) >= 20.0){    //muon/antimuon track masking condition

          elec_depo_energy_mask += depo.Energy();

          angle_barycentre_depo = mu_to_barycenter_vector.Angle(mu_to_depo_vector)*TMath::RadToDeg();                //angle between axis and hit in degrees
          angle_bin = floor(((angle_barycentre_depo*n_bin)/180));

          if(mu_to_depo_distance <= distance_max){   //containment sphere condition

            elec_histograms[0]->Fill(angle_barycentre_depo);

            //Fill angle distribution for sphere events 
            elec_energy_angular_distribution[angle_bin] += depo.Energy();                      
            
            if(angle_barycentre_depo <= 70){ // selection cone condition
              elec_depo_energy_cone += depo.Energy();
            } //end of cone condition     
          }   //end of containtment sphere condition
        }     //end of track masking condition
      }       //end of electron hit condition
    }         //end of depo loop


    elec_histograms[1]->Fill(elec_depo_energy_tot);
    elec_histograms[2]->Fill(elec_depo_energy_mask);
    elec_histograms[3]->Fill(elec_depo_energy_sphere);
    elec_histograms[4]->Fill(elec_depo_energy_cone);
    

    double elec_energy_ratio = 0.0;
    for(int i = 0; i < n_bin; i++){
      elec_energy_ratio += elec_energy_angular_distribution[i]/elec_depo_energy_sphere;              //portion of e energy inside cs + tm zone

      elec_energy_ratio_angular_distrubtion[i] = elec_energy_angular_distribution[i]/elec_depo_energy_sphere;      

      completeness[i] += elec_energy_ratio;  
    }

    hdiff->Fill(Posi_Energy[i_event]-elec_depo_energy_cone);          //energy difference between the mc and the selected one
    i_event++;   
  } //end of event loop

  cout << "# Events: " << i_event << endl;
  
  for(int i = 0; i < n_bin; i++){
    elec_energy_ratio_angular_distrubtion[i] /= i_event;
    completeness[i] /= i_event;
    gene_sel[0]->SetPoint(i, (i*180)/n_bin, elec_energy_ratio_angular_distrubtion[i]);
    gene_sel[1]->SetPoint(i, (i*180)/n_bin, completeness[i]);
  }

  //Get some numbers ?????????????????????????????????????????????????????????????????????????????????/
  // double y_me;

  // for(int i = 0; i < gene_sel[1]->GetN(); ++i){
  //   double x, y;
  //   gene_sel[1]->GetPoint(i, x, y);
  //   if(x == 70){y_me = y;}
  // }
  // cout << "Completeness (theta = 70Â°): " << y_me << endl;
  
  //Plot electron information
  TCanvas* ctheta = new TCanvas("ctheta", "Selected Theta Distribution");
  ctheta->cd();
  ctheta->SetLogy();
  elec_histograms[0]->SetStats(0);
  elec_histograms[0]->Draw(); 
  ctheta->Update();

  TCanvas* cene = new TCanvas("cene", "Total Energy Distribution");
  cene->cd();
  cene->SetLogy();  
  elec_histograms[1]->SetStats(0);
  elec_histograms[1]->Draw(); 
  cene->Update();

  TCanvas* cene_2 = new TCanvas("cene_2", "Muon Tracking Mask");
  cene_2->cd();
  cene_2->SetLogy();
  elec_histograms[2]->SetStats(0);
  elec_histograms[2]->Draw(); 
  cene_2->Update();

  TCanvas* cene_3 = new TCanvas("cene_3", "Muon Tracking Mask and Containing Sphere");
  cene_3->cd();
  cene_3->SetLogy();
  elec_histograms[3]->SetStats(0);
  elec_histograms[3]->Draw();; 
  cene_3->Update();

  TCanvas* cene_4 = new TCanvas("cene_4", "Muon Tracking Mask, Containing Sphere and Selection Cone");
  cene_4->cd();
  cene_4->SetLogy();
  elec_histograms[4]->SetStats(0);
  elec_histograms[4]->Draw(); 
  cene_4->Update();

  TCanvas* celec_sel = new TCanvas("celec_sel", "Selected Electron Information");
  celec_sel->Divide(1,2);
  for(int i = 0; i < gene_sel.size(); i++){
    celec_sel->cd(i+1); gene_sel[i]->Draw("AB");
  } 
  celec_sel->Update();

  //Decay mc positron information
  if(!muon){
    THStack* stack = new THStack("stack", "");
    stack->Add(hposi);
    stack->Add(elec_histograms[4]);

    TCanvas* cposi = new TCanvas("cposi", "Decay MC e^{+}");
    cposi->cd();
    cposi->SetLogy();
    stack->Draw("hist");

    //Legend
    TLegend* legend = new TLegend(0.125, 0.8, 0.225, 0.88);
    legend->AddEntry(hposi, " Decay MC Positron", "l");
    legend->AddEntry(elec_histograms[4], " Selected e^{+}/e^{-}", "l");
    legend->SetBorderSize(0);
    legend->Draw();

    TCanvas* cdiff = new TCanvas("cdiff", "Positron Energy Diff");
    cdiff->cd();
    cdiff->SetLogy();
    hdiff->SetStats(0);
    hdiff->Draw(); 
    cdiff->Update();

    cposi->SaveAs(("output/images/simmichelanalysis/"+generated+"_posi_energy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".pdf").c_str()); 
    cposi->SaveAs(("output/images/simmichelanalysis/"+generated+"_posi_energy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".root").c_str());

    cdiff->SaveAs(("output/images/simmichelanalysis/"+generated+"_diff_posi_energy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".pdf").c_str()); 
    cdiff->SaveAs(("output/images/simmichelanalysis/"+generated+"_diff_posi_energy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".root").c_str());
  }

  //Save it!
  ctheta->SaveAs(("output/images/simmichelanalysis/"+generated+"_theta_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".pdf").c_str()); 
  ctheta->SaveAs(("output/images/simmichelanalysis/"+generated+"_theta_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".root").c_str());

  cene->SaveAs(("output/images/simmichelanalysis/"+generated+"_Eenergy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".pdf").c_str()); 
  cene->SaveAs(("output/images/simmichelanalysis/"+generated+"_Eenergy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".root").c_str());

  cene_2->SaveAs(("output/images/simmichelanalysis/"+generated+"_MTMenergy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".pdf").c_str()); 
  cene_2->SaveAs(("output/images/simmichelanalysis/"+generated+"_MTMenergy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".root").c_str());

  cene_3->SaveAs(("output/images/simmichelanalysis/"+generated+"_CSenergy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".pdf").c_str()); 
  cene_3->SaveAs(("output/images/simmichelanalysis/"+generated+"_CSenergy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".root").c_str());

  cene_4->SaveAs(("output/images/simmichelanalysis/"+generated+"_Selenergy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".pdf").c_str()); 
  cene_4->SaveAs(("output/images/simmichelanalysis/"+generated+"_Selenergy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".root").c_str());

  celec_sel->SaveAs(("output/images/simmichelanalysis/"+generated+"_sel_elec_energy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".pdf").c_str()); 
  celec_sel->SaveAs(("output/images/simmichelanalysis/"+generated+"_sel_elec_energy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".root").c_str());
} // end of ElectronAnalysis

// bool debug = false;

// ----- M A I N ----- //

int main()
{
  //debug = Debug;
  //SetDebug(debug);
  //MuonAnred_g4_stage2_muon.rootorServicevolTPCActive", "generator::SinglesGen","test", true);


  string tag_energy_deposit = "largeant:LArG4DetectorServicevolTPCActive";
  string tag_MCTruth = "generator::SinglesGen";
  string tag_MCParticle = "";


  vector<pair<TVector3, int>> muon_analysis = MuonAnalysis(1, 1, 8, tag_energy_deposit, tag_MCTruth,"file", 1);

  ElectronAnalysis(1, 1, 8, tag_energy_deposit, tag_MCTruth, tag_MCParticle ,"file",muon_analysis, true);


  //ElectronAnalysis(n_files, i_first_event, i_last_event, tag_energy_deposit, tag_MCTruth, tag_MCParticle, listname, muon_analysis, muon);
  return 0;
}