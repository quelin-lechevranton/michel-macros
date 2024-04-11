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






// TVector3 to_TVector3(auto geotype) {
//   TVector3 vector;
//   vector.SetXYZ( geotype.X() , geotype.Y(), geotype.Z() );
//   return vector;
// }


double geoAngle(geo::Vector_t u, geo::Vector_t v) {
  return TMath::ACos(u.Dot(v)/u.R()/v.R())*TMath::RadToDeg;
}


/************************************ MuonAnalysis ************************************/

vector<pair<geo::Vector_t, int>> MuonAnalysis(int n_files, int i_first_event, int i_last_event, string depotag, string MCTruthtag, string listname, bool is_muon, bool save, bool verbose)
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

  //Graph
  // TGraph* Graph_muon;
  // Graph_muon = new TGraph();
  // Graph_muon->SetName("remaining_energy_per_track_length");
  // Graph_muon->SetMarkerColor(colors[4]);
  // Graph_muon->SetLineColor(colors[4]);
  // Graph_muon->SetMarkerStyle(5);
  // Graph_muon->SetMarkerSize(0.5);
  // Graph_muon->GetXaxis()->SetTitle("Track Length (cm)");
  // Graph_muon->GetYaxis()->SetTitle(("Remaining "+particle_name+" Energy (MeV)").c_str());


  int n_bin_X=450, n_bin_Y=1000; double X_min=0., X_max=450., Y_min=0., Y_max=1000.; 
  auto Histo_muon_remE = new TH2D("Histo_muon_remE","Remaining Energy (MeV);Track Length (cm)",n_bin_X,X_min,X_max,n_bin_Y,Y_min,Y_max);

  n_bin_X=320, n_bin_Y=100; X_min=0., X_max=160., Y_min=0., Y_max=10.; 
  auto Histo_muon_dEdx = new TH2D("Histo_muon_dEdx","dE/dx (MeV/cm);Residual Track Length (cm)",n_bin_X,X_min,X_max,n_bin_Y,Y_min,Y_max);

  //Tag conversion
  art::InputTag depo_tag(depotag);
  art::InputTag MCTruth_tag(MCTruthtag);  


  //Search your files
  string file_list          = listname+".list";
  vector<string> filenames = ReadFileList(n_files, file_list);

  int i_muon_total_depo = 0;    //to account for all points from all events
  int i_event = 0;
  vector<pair<geo::Vector_t, int>> muon_analysis(i_last_event-(i_first_event-1));

  int depoless_count = 0; //number of depoless events

  //Go for events!
  for(gallery::Event ev(filenames); !ev.atEnd(); ev.next() , i_event++ ){
    //Physical quantities 
    vector<double> array_muon_track_length(0);       //hit position from start hit
    vector<double> array_muon_depo_energy(0);        //energy deposition in each hit
    vector<double> array_muon_remaining_energy(0);   //remaining energy in each hit

    if(ev.eventAuxiliary().event() < i_first_event) continue;
    if(ev.eventAuxiliary().event() > i_last_event) break; 

    if (verbose) {cout << "Event: " << ev.eventAuxiliary().event() << endl;}
      
    //Retrieve list of energy deposits per event
    auto const depo_list = ev.getValidHandle<vector<sim::SimEnergyDeposit>>(depo_tag);
    size_t n_depo        = depo_list->size();

    if (verbose && n_depo ==0) {cout << "There is no deposit here!" << endl; depoless_count++; continue;}

    //Retrieve generated information   
    auto const montecarlo        = ev.getValidHandle<vector<simb::MCTruth>>(MCTruth_tag);
    const simb::MCTruth& gen     = montecarlo->at(0);
    const simb::MCParticle& muon = gen.GetParticle(0);

    double muon_remaining_energy = (1000)*muon.E();   //in MeV
    if (verbose) {cout << "Muon/AntiMuon Energy: " << muon_remaining_energy << " MeV" << endl;}
    double muon_depo_energy_tot  = 0;
    int    i_muon_last_hit       = -1;

    geo::Point_t muon_first_hit = depo_list->at(0).MidPoint();
    geo::Point_t muon_previous_hit = muon_first_hit;
    double muon_track_length = 0.;

    //Depo loop for muon/antimuon information
    for(size_t i_depo = 0; i_depo < n_depo; i_depo++){

      const sim::SimEnergyDeposit& depo = depo_list->at(i_depo);
      geo::Point_t hit = depo.MidPoint();

      if(depo.PdgCode() != pdg) {continue;} //muon/anitmuon condition

      muon_track_length += (muon_previous_hit - hit).R();
      muon_remaining_energy -= depo.Energy();
      muon_depo_energy_tot  += depo.Energy();
      
      //Energy deposit and position for muon/antimuon
      array_muon_track_length.push_back( muon_track_length );
      array_muon_depo_energy.push_back( depo.Energy() );
      array_muon_remaining_energy.push_back(muon_remaining_energy);

      i_muon_last_hit = i_depo;              //save last hit
      muon_previous_hit = hit;
    }   //end of depo loop 
    
    if(array_muon_track_length.size() == 0 && verbose){cout << "There is no muon/antimuon information here!" << endl; continue;}    //go back to begining if you don't have muon information
    
    geo::Point_t muon_last_hit = depo_list->at(i_muon_last_hit).MidPoint();

    if (verbose) {
        cout << "Muon/AntiMuon Total Deposit Energy: " << muon_depo_energy_tot << " MeV" << endl;
        cout << "Muon/AntiMuon Last Hit: ( " << muon_last_hit.X() << " , " << muon_last_hit.Y() << " , " << muon_last_hit.Z() << " )" << endl;
    }


    //Direction vector  
    geo::Vector_t muon_track_global_direction = ( muon_first_hit - muon_last_hit ).Unit();

    //Retrieve residual track length
    int n_track_length = array_muon_track_length.size();
    vector<double> array_muon_residual_track(n_track_length);  
    for (int i=0 ; i<n_track_length ; i++ ){
      array_muon_residual_track[i]=array_muon_track_length[n_track_length-1-i];
    }
  
    //Compute dE/dx in MeV/cm
    double cm_per_bin = 0.03;
    vector<double> array_muon_dEdx( array_muon_depo_energy.size() );
    for (int i=0 ; i < array_muon_depo_energy.size() ; i++ ) {
      array_muon_dEdx[i] = 1/cm_per_bin * array_muon_depo_energy[i];
    } 


    //Graph

    for (int i = 0; i < array_muon_track_length.size(); i++ , i_muon_total_depo++ ){   
        Histo_muon_dEdx->Fill(array_muon_residual_track[i], array_muon_dEdx[i]);
        Histo_muon_remE->Fill(array_muon_track_length[i]  , array_muon_remaining_energy[i]);
        // Graph_muon->SetPoint(i_muon_total_depo, array_muon_track_length[i]  , array_muon_remaining_energy[i]);
    }   //end of muon/antimuon information loop


    pair<geo::Vector_t, int> p(muon_track_global_direction, i_muon_last_hit);
    muon_analysis[i_event] = p;
  } //end of event loop

  cout << depoless_count << " events did not have any deposited energy" << endl;

  printf("####################################");

  //Plot muon/antimuon information
  TCanvas* canvas_muon_1 = new TCanvas("canvas_muon_1", "Muon/AntiMuon dE/dx");
  canvas_muon_1->cd();
  Histo_muon_dEdx->Draw("Colz");
  canvas_muon_1->Update();

  TCanvas* canvas_muon_2 = new TCanvas("canvas_muon_2", "Muon/AntiMuon Remaining Energy");
  canvas_muon_2->cd();
  Histo_muon_remE->Draw("Colz");
  canvas_muon_2->Update();

  //Save it!
  if (save) {
    canvas_muon_1->SaveAs(("output/"+generated+"_depo_energy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".pdf").c_str()); 
    canvas_muon_1->SaveAs(("output/"+generated+"_depo_energy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".root").c_str());

    canvas_muon_2->SaveAs(("output/"+generated+"_remaining_energy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".pdf").c_str()); 
    canvas_muon_2->SaveAs(("output/"+generated+"_remaining_energy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".root").c_str());
  }

  return muon_analysis;

} // end of MuonAnalysis

/************************************ ElectronAnalysis ************************************/


void ElectronAnalysis(int n_files, int i_first_event, int i_last_event, string depotag, string MCTruthtag, string MCParticletag, string listname, vector<pair<geo::Vector_t, int>> muon_analysis, bool is_muon, bool save, bool verbose)
{
  //X-axis
  vector<string> xtitle = {"#theta (degrees)", "Energy (MeV)", "Energy (MeV)", "Energy (MeV)", "Energy (MeV)"};
  vector<int> n_bin     = {180, 90, 90, 90, 90};
  vector<float> X_min    = {-10, 0, 0, 0, 0};
  vector<float> X_max    = {200, 360, 100, 60, 60};
  

  //Particle options
  vector<string> pdgnames = {"e^{-}", "e^{+}", "#mu^{-}", "#mu^{+}", "p^{+}"};
  vector<int> pdgcodes    = {11, -11, 13, -13, 2212};
  Color_t colors[7]       = {kBlue+2, kAzure+8, kOrange+7, kRed+2, kGreen-3, kMagenta+2, kBlue-6};

  //To distinguish between muon/antimuon generated
  int pdg;
  string particle_name;
  string generated;
  if(is_muon) {pdg = pdgcodes[2]; particle_name = pdgnames[2]; generated = "muon";}
  if(!is_muon){pdg = pdgcodes[3]; particle_name = pdgnames[3]; generated = "antimuon";}

  vector<TH1D*> elec_histograms(5);
  for(int i = 0; i < elec_histograms.size(); i++){
    elec_histograms[i] = new TH1D(Form("elec_histograms_%d", i), Form(";%s;Counts", xtitle[i].c_str()), n_bin[i], X_min[i], X_max[i]);
    if(is_muon) {elec_histograms[i]->SetLineColor(colors[6]);  elec_histograms[i]->SetFillColor(colors[6]);}
    if(!is_muon){elec_histograms[i]->SetLineColor(colors[1]);  elec_histograms[i]->SetFillColor(colors[1]);} 
  }
  elec_histograms[0]->SetLineColor(colors[5]);  elec_histograms[0]->SetFillColor(colors[5]);

  vector<TGraph*> gene_sel(2);
  for(int i = 0; i < gene_sel.size(); i++){
    gene_sel[i] = new TGraph();
    gene_sel[i]->SetName(Form("gene_sel_%d", i));
    if(is_muon) {gene_sel[i]->SetMarkerColor(colors[6]);  gene_sel[i]->SetLineColor(colors[6]); gene_sel[i]->SetFillColor(colors[6]);}
    if(!is_muon){gene_sel[i]->SetMarkerColor(colors[1]);  gene_sel[i]->SetLineColor(colors[1]); gene_sel[i]->SetFillColor(colors[1]);} 
    gene_sel[i]->SetMarkerStyle(3);
    gene_sel[i]->SetMarkerSize(0.25);
    gene_sel[i]->GetXaxis()->SetTitle("#theta (degrees)");
  }
  gene_sel[0]->GetYaxis()->SetTitle("Electron Energy Ratio Angular (Barycenter - Muon Last Hit - Depo) Distribution");
  gene_sel[1]->GetYaxis()->SetTitle("Completeness");

  TH1D* hposi = new TH1D("", Form("Decay MC e^{+};%s;Counts", xtitle[4].c_str()), n_bin[4], X_min[4], X_max[4]);
  hposi->SetLineColor(colors[4]);  hposi->SetFillColor(colors[4]);

  TH1D* hdiff = new TH1D("", Form("E_{e^{+}, MC} - E_{e^{+}, Selected};%s;Counts", xtitle[4].c_str()), n_bin[4], -X_max[4], X_max[4]);
  hdiff->SetLineColor(colors[3]);  //hdiff->SetFillColor(colors[3]);

  //All tags that you need
  art::InputTag depo_tag(depotag);
  art::InputTag MCTruth_tag(MCTruthtag);
  //art::InputTag MCParticle_tag(MCParticletag); 
  //Search your files
  string file_list          = listname+".list";
  vector<string> filenames = ReadFileList(n_files, file_list);

  vector<double> Posi_Energy(i_last_event-(i_first_event-1));           //to store each decay positron energy from mc
  vector<float> elec_energy_angular_distribution(n_bin[0]);  
  vector<float> elec_energy_ratio_angular_distrubtion(n_bin[0]);
  vector<float> completeness(n_bin[0]);
  //int elec_pc = 0;                                          //to account for all points from all events
  int i_event    = 0;


  int depoless_count =0; //track the number of depoless events
  
  //Go for events!
  for(gallery::Event ev(filenames); !ev.atEnd(); ev.next(), i_event++){
    //Physical quantities
    // vector<double> elecDist;    //distance of electron hit from last muon hit
    // vector<double> elecDepoE;   //deposited energy of selected electrons

    if(ev.eventAuxiliary().event() < i_first_event) continue;
    if(ev.eventAuxiliary().event() > i_last_event) break; 

    if (verbose) {cout << "Event: " << ev.eventAuxiliary().event() << endl;}
      
    //Retrieve list of energy deposits per event
    auto const depo_list = ev.getValidHandle<vector<sim::SimEnergyDeposit>>(depo_tag);
    size_t n_depo       = depo_list->size();

    if (verbose && n_depo ==0) {cout << "There is no deposit here!" << endl;depoless_count++;continue;}

    //Decay mc positron information ~ only for antimuon events
    /*
    if(!muon){
      auto const posicarlo = ev.getValidHandle<vector<simb::MCParticle>>(MCParticle_tag);
      size_t nposi         = posicarlo->size();

      for(size_t pos_i = 0; pos_i < nposi; pos_i++){
        const simb::MCParticle& positron = posicarlo->at(pos_i);
        if(positron.PdgCode() == -11 && positron.Process() == "Decay"){Posi_Energy[i_event] = (1000)*positron.E(); hposi->Fill(Posi_Energy[i_event]);}
      }    
    } */ 

    if(verbose && muon_analysis[i_event].second == 0){cout << "There is no muon information here!" << endl; continue;}    //go back to begining if you don't have muon information

    //Information from muon/anitmuon analysis     
    geo::Vector_t muon_track_global_direction = muon_analysis[i_event].first;      
    int i_muon_last_hit = muon_analysis[i_event].second;

    geo::Point_t muon_last_hit = depo_list->at(i_muon_last_hit).MidPoint();



    //Electron information
    double distance_max  = 25.0;        //max. distance from last muon hit (in cm) ~ electron analysis
    double elec_depo_energy_tot   = 0;
    double elec_depo_energy_mask  = 0; 
    double elec_depo_energy_sphere   = 0;
    double elec_depo_energy_cone = 0;       

    geo::Point_t depo_weighted_point;     //weighted electron hits
    
    //Compute barycenter
    for(size_t i_depo = 0; i_depo < n_depo; i_depo++){

      const sim::SimEnergyDeposit& depo = depo_list->at(i_depo);
      
      if(abs(depo.PdgCode()) == 11){  

        //Retrieve position of energy deposit
        geo::Point_t depo_hit = depo.MidPoint();
        geo::Vector_t mu_to_depo_vector = depo_hit - muon_last_hit;

        //Distance and angle of electron hit from muon last hit
        double mu_to_depo_distance = mu_to_depo_vector.R();
        double angle_track_depo  = geoAngle(muon_track_global_direction,mu_to_depo_vector); 

        if(abs(angle_track_depo) >= 20.0){     //muon/antimuon track masking 

          if(mu_to_depo_distance <= distance_max){

            elec_depo_energy_sphere += depo.Energy();
            depo_weighted_point  += depo.Energy()*depo_hit;

          } //end of sphere condition
        }   //end of track masking condition
      }     //end of electron condition
    }       //end of depo loop


    if(verbose && elec_depo_energy_sphere <= 0){cout << "There is no Michel electron information here!" << endl; continue;}    //go back to begining if you don't have Michel electron information

    //coordinates of barycenter shower
    geo::Point_t elec_barycenter = depo_weighted_point /elec_depo_energy_sphere; 

    if (verbose) {cout << "Barycenter: (" << elec_barycenter.X() << "," << elec_barycenter.Y() << "," << elec_barycenter.Z() << ")" << endl;}

    //axis cone vector
    geo::Vector_t mu_to_barycenter_direction = ( elec_barycenter - muon_last_hit ).Unit(); 

    //Depo loop for michel electron selection
    for(size_t i_depo = 0; i_depo < n_depo; i_depo++){
      
      const sim::SimEnergyDeposit& depo = depo_list->at(i_depo);
      
      if(abs(depo.PdgCode()) == 11){        

        elec_depo_energy_tot += depo.Energy();  

        //Retrieve position of energy deposit
        geo::Point_t depo_hit = depo.MidPoint();
        
        geo::Vector_t mu_to_depo_vector = depo_hit - muon_last_hit;

        //Distance and angle of electron hit from muon last hit
        double mu_to_depo_distance = mu_to_depo_vector.R();
        double angle_track_depo  = geoAngle(muon_track_global_direction , mu_to_depo_vector); 


        if(abs(angle_track_depo) >= 20.0){    //muon/antimuon track masking condition

          elec_depo_energy_mask += depo.Energy();

          if(mu_to_depo_distance <= distance_max){   //containment sphere condition

            //angle between axis and hit in degrees
            double angle_barycentre_depo = geoAngle( mu_to_barycenter_direction , mu_to_depo_vector);                
            elec_histograms[0]->Fill(angle_barycentre_depo);


            //Fill angle distribution for sphere events 
            int angle_bin = floor( n_bin[0]*angle_barycentre_depo / 180);
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
    for(int i = 0; i < n_bin[0]; i++){

      elec_energy_ratio_angular_distrubtion[i] = elec_energy_angular_distribution[i]/elec_depo_energy_sphere;       
      elec_energy_ratio += elec_energy_ratio_angular_distrubtion[i];              //portion of e energy inside cs + tm zone    
      completeness[i] += elec_energy_ratio;  
    }

    hdiff->Fill(Posi_Energy[i_event]-elec_depo_energy_cone);          //energy difference between the mc and the selected one 
  } //end of event loop

  cout << "# Events: " << i_event << endl;
  cout << depoless_count << " events did not have any deposited energy" << endl;
  cout << "##################################" << endl;
  
  for(int i = 0; i < n_bin[0]; i++){
    elec_energy_ratio_angular_distrubtion[i] /= i_event;
    completeness[i] /= i_event; 
    gene_sel[0]->SetPoint(i, (i*180)/n_bin[0], elec_energy_ratio_angular_distrubtion[i]);
    gene_sel[1]->SetPoint(i, (i*180)/n_bin[0], completeness[i]);
  }

  //Get some numbers ?????????????????????????????????????????????????????????????????????????????????/
  // double y_me;

  // for(int i = 0; i < gene_sel[1]->GetN(); ++i){
  //   double x, y;
  //   gene_sel[1]->GetPoint(i, x, y);
  //   if(x == 70){y_me = y;}
  // }
  // cout << "Completeness (theta = 70°): " << y_me << endl;
  
  //Plot electron information
  TCanvas* ctheta = new TCanvas("ctheta", "Electrons Angular (Barycenter - Muon Last Hit - Depo) Distribution");
  ctheta->cd();
  ctheta->SetLogy();
  elec_histograms[0]->SetStats(0);
  elec_histograms[0]->Draw(); 
  ctheta->Update();

  TCanvas* cene = new TCanvas("cene", "Electron Total Energy Distribution");
  cene->cd();
  cene->SetLogy();  
  elec_histograms[1]->SetStats(0);
  elec_histograms[1]->Draw(); 
  cene->Update();

  TCanvas* cene_2 = new TCanvas("cene_2", "e^{-} E dist, Muon Tracking Mask");
  cene_2->cd();
  cene_2->SetLogy();
  elec_histograms[2]->SetStats(0);
  elec_histograms[2]->Draw(); 
  cene_2->Update();

  TCanvas* cene_3 = new TCanvas("cene_3", "e^{-} E dist, Muon Tracking Mask and Containing Sphere");
  cene_3->cd();
  cene_3->SetLogy();
  elec_histograms[3]->SetStats(0);
  elec_histograms[3]->Draw();; 
  cene_3->Update();

  TCanvas* cene_4 = new TCanvas("cene_4", "e^{-} E dist, Muon Tracking Mask, Containing Sphere and Selection Cone");
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
  // if(!is_muon){
  //   THStack* stack = new THStack("stack", "");
  //   stack->Add(hposi);
  //   stack->Add(elec_histograms[4]);

  //   TCanvas* cposi = new TCanvas("cposi", "Decay MC e^{+}");
  //   cposi->cd();
  //   cposi->SetLogy();
  //   stack->Draw("hist");

  //   //Legend
  //   TLegend* legend = new TLegend(0.125, 0.8, 0.225, 0.88);
  //   legend->AddEntry(hposi, " Decay MC Positron", "l");
  //   legend->AddEntry(elec_histograms[4], " Selected e^{+}/e^{-}", "l");
  //   legend->SetBorderSize(0);
  //   legend->Draw();

  //   TCanvas* cdiff = new TCanvas("cdiff", "Positron Energy Diff");
  //   cdiff->cd();
  //   cdiff->SetLogy();
  //   hdiff->SetStats(0);
  //   hdiff->Draw(); 
  //   cdiff->Update();

  //   cposi->SaveAs(("output/images/simmichelanalysis/"+generated+"_posi_energy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".pdf").c_str()); 
  //   cposi->SaveAs(("output/images/simmichelanalysis/"+generated+"_posi_energy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".root").c_str());

  //   cdiff->SaveAs(("output/images/simmichelanalysis/"+generated+"_diff_posi_energy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".pdf").c_str()); 
  //   cdiff->SaveAs(("output/images/simmichelanalysis/"+generated+"_diff_posi_energy_"+to_string(i_first_event)+"_"+to_string(i_last_event)+".root").c_str());
  // }

  //Save it!
  if (save) {
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
  }

} // end of ElectronAnalysis

// bool debug = false;

// ----- M A I N ----- //

int SimMichelAnalysis(int n_files, int i_first_event, int i_last_event, bool save, bool verbose)
{
  //debug = Debug;
  //SetDebug(debug);
  //MuonAnred_g4_stage2_muon.rootorServicevolTPCActive", "generator::SinglesGen","test", true);


  string depotag = "largeant:LArG4DetectorServicevolTPCActive";
  string MCTruthtag = "generator::SinglesGen";
  string MCParticletag = "";


  vector<pair<geo::Vector_t, int>> muon_analysis = MuonAnalysis(n_files, i_first_event, i_last_event, depotag, MCTruthtag,"list", 1, save, verbose);

  ElectronAnalysis(n_files, i_first_event, i_last_event, depotag, MCTruthtag, MCParticletag ,"list",muon_analysis, true, save, verbose);

  return 0;
}
