namespace detinfo {
class DetectorClocksData; }

namespace dune{


  class MichelStudy : public art::EDAnalyzer 
  {
  public:

    explicit MichelStudy(fhicl::ParameterSet const& pset);
    virtual ~MichelStudy();

    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);

  }; 

  //========================================================================
 
  void MichelStudy::analyze( const art::Event& evt)
  {

   reset();  
   
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  art::ServiceHandle<cheat::PhotonBackTrackerService> pbt_serv;

  bool isMC = !evt.isRealData();
 
  geo::GeometryCore const * fGeometry = &*(art::ServiceHandle<geo::Geometry>());

  // *MC truth cosmic generator information
  art::Handle< std::vector<simb::MCTruth> > mctruthcryListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclistcry;
  if (isMC )
  {
    if (evt.getByLabel(fCryGenModuleLabel,mctruthcryListHandle))
    {
      art::fill_ptr_vector(mclistcry, mctruthcryListHandle);
    }
    else
    {
      // If we requested this info but it doesn't exist, don't use it.
      //fSaveCryInfo = false;
      mf::LogError("AnalysisTree") << "Requested CRY information but none exists, hence not saving CRY information.";
    }
  }

  // * MC truth ProtoDUNE beam generator information
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclistgen;
  if (isMC)
  {
    if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
    art::fill_ptr_vector(mclistgen, mctruthListHandle);
  }

  // *MC truth G4 information
  art::Handle< std::vector<simb::MCParticle> > mcparticleListHandle;
  std::vector<art::Ptr<simb::MCParticle> > mcparticlelist;
  if (isMC )
  {
    if (evt.getByLabel(fG4ModuleLabel,mcparticleListHandle))
    {
      art::fill_ptr_vector(mcparticlelist, mcparticleListHandle);
    }
    else
    {
      // If we requested this info but it doesn't exist, don't use it.
      //fSaveCryInfo = false;
      mf::LogError("AnalysisTree") << "Requested G4 information but none exists, hence not saving G4 information.";
    }
  }


  if (isMC)
  {
    // *MC Shower information
    art::Handle< std::vector<sim::MCShower> > mcshowerh;
    evt.getByLabel(fMCShowerModuleLabel, mcshowerh);

    // *MC Track information
    art::Handle< std::vector<sim::MCTrack> > mctrackh;
    evt.getByLabel(fMCTrackModuleLabel, mctrackh);

//  int nMCShowers = 0;
//  if (mcshowerh.isValid())
//  nMCShowers = mcshowerh->size();

//  int nMCTracks = 0;
//  if (mctrackh.isValid())
//  nMCTracks = mctrackh->size();

    art::Ptr<simb::MCTruth> mctruthcry;
    int nCryPrimaries = 0, nCryOrigin = 0;
    mctruthcry = mclistcry[0];
    nCryPrimaries = mctruthcry->NParticles();
    nCryOrigin =   mctruthcry->Origin();

    art::Ptr<simb::MCTruth> mctruthgen;
    int nGenPrimaries = 0, nGenOrigin = -99;
    mctruthgen = mclistgen[0];
    nGenPrimaries = mctruthgen->NParticles();
    nGenOrigin =   mctruthgen->Origin();
    
    cout<<"nCryPrimaries "<<nCryPrimaries<<" nCryOrigin "<<nCryOrigin<<"\n";
  
    bool isCosmics = false;
  
    //find origin
    std::vector< art::Handle< std::vector<simb::MCTruth> > > allmclists;
    evt.getManyByType(allmclists);
    for(size_t mcl = 0; mcl < allmclists.size(); ++mcl)
    {
      art::Handle< std::vector<simb::MCTruth> > mclistHandle = allmclists[mcl];
      for(size_t m = 0; m < mclistHandle->size(); ++m)
      {
	art::Ptr<simb::MCTruth> mct(mclistHandle, m);
	isCosmics = false;
	if (mct->Origin() == simb::kCosmicRay) isCosmics = true;
	cout<<"mct->Origin() "<<mct->Origin()<<" isCosmics "<<isCosmics<<"\n";
      }
    }
     
    //We have one single gen entry, and one cosmic entry per event
    
    // GENIE
    if (!mclistgen.empty())
    {//at least one mc record

      const sim::ParticleList& plist = pi_serv->ParticleList();
      int nGEANTparticles = plist.size();

       cout<< "Expected GEANT particles: "<< nGEANTparticles << " ,GENIE particles: "<< nGenPrimaries <<" nGenOrigin "<<nGenOrigin<< "\n";

    } //if (!mclistgen.empty())
     
  //************Now coming to reco quantities***************//
  
  // tracks
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if(evt.getByLabel(fTrackModuleLabel,trackListHandle)) art::fill_ptr_vector(tracklist, trackListHandle);
    
  // showers
  art::Handle< std::vector<recob::Shower> > showerListHandle;
  std::vector<art::Ptr<recob::Shower> > showerlist;
  if(evt.getByLabel(fShowerModuleLabel,showerListHandle)) art::fill_ptr_vector(showerlist, showerListHandle);

  // PFParticles
  art::Handle< std::vector<recob::PFParticle> > PFPListHandle; 
  std::vector<art::Ptr<recob::PFParticle> > pfplist;
  if(evt.getByLabel(fPFPModuleLabel,PFPListHandle)) art::fill_ptr_vector(pfplist, PFPListHandle);
    
  // hits
  art::Handle< std::vector<recob::Hit> > hitListHandle; // to get information about the hits   
  std::vector<art::Ptr<recob::Hit>> hitlist;
  if(evt.getByLabel(fHitsModuleLabel, hitListHandle)) art::fill_ptr_vector(hitlist, hitListHandle);
        
  // Space Points
  art::Handle<std::vector<recob::SpacePoint>> spacepointListHandle;
  std::vector<art::Ptr<recob::SpacePoint> > spacepointList;
  if (evt.getByLabel(fSpacePointSolverModuleLabel,spacepointListHandle)) art::fill_ptr_vector(spacepointList, spacepointListHandle);

  // clusters
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  std::vector<art::Ptr<recob::Cluster> > clusterlist;
  if (evt.getByLabel(fClusterModuleLabel,clusterListHandle)) art::fill_ptr_vector(clusterlist, clusterListHandle);

  // Internal flashes
  art::Handle< std::vector<recob::OpFlash> > flashListHandleInt;
  std::vector<art::Ptr<recob::OpFlash> > flashlistInt;
  if (evt.getByLabel(fOpFlashModuleLabelInt,flashListHandleInt)) art::fill_ptr_vector(flashlistInt, flashListHandleInt);

  // External flashes
  art::Handle< std::vector<recob::OpFlash> > flashListHandleExt;
  std::vector<art::Ptr<recob::OpFlash> > flashlistExt;
  if (evt.getByLabel(fOpFlashModuleLabelExt,flashListHandleExt)) art::fill_ptr_vector(flashlistExt, flashListHandleExt);

  // Internal Ophits
  art::Handle< std::vector<recob::OpHit> > ophitListHandleInt;
  std::vector<art::Ptr<recob::OpHit> > ophitlistInt;
  if (evt.getByLabel(fOpHitModuleLabelInt,ophitListHandleInt)) art::fill_ptr_vector(ophitlistInt, ophitListHandleInt);

  // External Ophits
  art::Handle< std::vector<recob::OpHit> > ophitListHandleExt;
  std::vector<art::Ptr<recob::OpHit> > ophitlistExt;
  if (evt.getByLabel(fOpHitModuleLabelExt,ophitListHandleExt)) art::fill_ptr_vector(ophitlistExt, ophitListHandleExt);

  // Sim channels
  art::Handle< std::vector<sim::SimChannel> > simchannelListHandle;
  std::vector<art::Ptr<sim::SimChannel> > simchannellist;
  if (evt.getByLabel(fSimChannelModuleLabel,simchannelListHandle)) art::fill_ptr_vector(simchannellist, simchannelListHandle);

  // Associations 
  art::FindMany<recob::Hit>                          fmhit           (trackListHandle, evt, fTrackModuleLabel);
  art::FindMany<recob::Hit>                          fshwrhit        (showerListHandle, evt, fShowerModuleLabel);
  art::FindManyP<recob::Hit>                         fshwrHit        (showerListHandle, evt, fShowerModuleLabel);
  art::FindManyP<recob::Hit, recob::TrackHitMeta>    fmthm           (trackListHandle, evt, fTrackModuleLabel); // to associate tracks and hits
  art::FindManyP<recob::Hit>                         fmtht           (trackListHandle, evt, fTrackModuleLabel); // to associate tracks and hits
  art::FindManyP<recob::Track>                       thass           (hitListHandle, evt, fTrackModuleLabel); //to associate hits with tracks
  art::FindManyP<anab::Calorimetry>                  fmcal           (trackListHandle, evt, fCalorimetryModuleLabel); //to associate tracks with calorimetry
  art::FindManyP<anab::T0>                           trk_t0_assn_v   (PFPListHandle, evt ,fPFPModuleLabel); //to associate PFParticles with T0
  art::FindManyP<anab::T0>                           shwr_t0_assn_v  (PFPListHandle, evt ,fPFPModuleLabel); //to associate PFParticles with T0
  art::FindManyP<recob::PFParticle>                  pfp_trk_assn    (trackListHandle, evt ,fTrackModuleLabel); //to associate tracks with PFParticles
  art::FindManyP<recob::PFParticle>                  pfp_shwr_assn   (showerListHandle, evt ,fShowerModuleLabel); //to associate showers with PFParticles
  art::FindManyP<anab::T0>                           fmT0            (trackListHandle, evt ,"pmtrack"); //to associate tracks with T0
  art::FindManyP<recob::Cluster>                     fmcl            (hitListHandle,evt,fClusterModuleLabel);
  art::FindManyP<recob::Hit>                         hitsFromClusters(clusterListHandle, evt, fClusterModuleLabel);
  art::FindManyP<recob::SpacePoint>                  fmsp            (hitListHandle,evt,fSpacePointSolverModuleLabel);
  art::FindManyP<recob::OpHit>                       fophitflsh      (flashListHandleInt, evt, fOpFlashModuleLabelInt);
  art::FindManyP<recob::Wire>                        wFromHits       (hitListHandle,evt,fHitsModuleLabel);
        
  anab::MVAReader<recob::Hit,4> hitResults(evt, "emtrkmichelid:emtrkmichel" );

  // Run info	
  frun = evt.run();
  fsubrun = evt.subRun();
  fevent = evt.id().event();
  art::Timestamp ts = evt.time();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  fevttime=tts.AsDouble();
     
  // Use '_detp' to find 'efield' and 'temp'
  auto const* _detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  double efield = _detp -> Efield();
  double temp   = _detp -> Temperature();
  // Determine the drift velocity from 'efield' and 'temp'
  fDriftVelocity = _detp -> DriftVelocity(efield,temp);

    UInt_t year=0;
    UInt_t month=0;
    UInt_t day=0;
     
    fyear_month_date=tts.GetDate(kTRUE,0,&year,&month,&day);
     
    UInt_t hour=0;
    UInt_t min=0;
    UInt_t sec=0;
    
    fhour_min_sec=tts.GetTime(kTRUE,0,&hour,&min,&sec);
    
    cout<<"year_month_date: "<<fyear_month_date<<" hour_min_sec: "<<fhour_min_sec<<"\n";
    
    //External flashes info
    cout<<"flashlistExt.size() "<<flashlistExt.size()<<"\n";
    fno_flashes_ext =  0; 
    fext_trigger_time = 0;
    double biggest_PE = 0;    
    for (size_t i = 0; i < flashlistExt.size() && i < kMaxFlashes ; ++i)
    {//loop over flashes
      art::Ptr<recob::OpFlash> pflashExt(flashListHandleExt, i);
      const recob::OpFlash& flashExt = *pflashExt;

      fflash_time_ext[fno_flashes_ext]       = flashExt.Time();
      fflash_pe_ext[fno_flashes_ext]         = flashExt.TotalPE();
      fflash_ycenter_ext[fno_flashes_ext]    = flashExt.YCenter();
      fflash_zcenter_ext[fno_flashes_ext]    = flashExt.ZCenter();
      fflash_ywidth_ext[fno_flashes_ext]     = flashExt.YWidth();
      fflash_zwidth_ext[fno_flashes_ext]     = flashExt.ZWidth();
      fflash_timewidth_ext[fno_flashes_ext]  = flashExt.TimeWidth();
      fflash_abstime_ext[fno_flashes_ext]    = flashExt.AbsTime(); // Time by PMT readout clock (?)
      fflash_frame_ext[fno_flashes_ext]      = flashExt.Frame(); // Frame number 
      
      if(flashExt.TotalPE()>biggest_PE)
      { 
        biggest_PE = flashExt.TotalPE();
        fext_trigger_time = flashExt.Time();
      }

      fno_flashes_ext++;
    }
    
    // loop through all reco'd internal flash times
    double dt_min = 99999.; // us
    
    cout<<"flashlistInt.size() "<<flashlistInt.size()<<"\n";
    fno_flashes_int=0;
    
    for (size_t nn = 0; nn < flashlistInt.size() && nn < kMaxFlashes ; ++nn)
    {
      art::Ptr<recob::OpFlash> pflashInt(flashListHandleInt, nn);
      const recob::OpFlash& flashInt = *pflashInt;
 
      //fall_flash_time_diff[fno_flashes_int]     = 3*(flashInt.Time()-fext_trigger_time)+50 - fseltrkrecotime[fsel_mu]; // Time diff relative to trigger
	  	    
      fflash_time_int[fno_flashes_int]       = flashInt.Time(); // Time relative to trigger
      fflash_pe_int[fno_flashes_int]         = flashInt.TotalPE();
      fflash_ycenter_int[fno_flashes_int]    = flashInt.YCenter(); // Geometric center in y
      fflash_zcenter_int[fno_flashes_int]    = flashInt.ZCenter(); // Geometric center in z
      fflash_ywidth_int[fno_flashes_int]     = flashInt.YWidth(); // Geometric width in y
      fflash_zwidth_int[fno_flashes_int]     = flashInt.ZWidth(); // Geometric width in z
      fflash_timewidth_int[fno_flashes_int]  = flashInt.TimeWidth(); // Width of the flash in time
      fflash_abstime_int[fno_flashes_int]    = flashInt.AbsTime(); // Time by PMT readout clock (?)
      fflash_frame_int[fno_flashes_int]      = flashInt.Frame(); // Frame number 
		  
      fno_flashes_int++;	  

    }// loop over flashes
    

    //********Storing hits of all long tracks, will use later*****************//

    std::vector <double> trkHitsKey;
//    double trkHitsKey[100000] = {-999};
//    int longtrk_hits = 0;
    //size_t total=0;
    // int mmm = 0 ;
    for(size_t oo=0; oo<tracklist.size(); oo++)
    {
      art::Ptr<recob::Track> ptrack_l(trackListHandle, oo);
      const recob::Track& track_l = *ptrack_l;
      auto vhit=fmthm.at(oo);
      auto vmeta=fmthm.data(oo);

//      std::vector<const recob::Hit* > hitsll = fmhit.at(oo);
//      cout<<"vhit.size() "<<vhit.size()<<"\n";	
      if(track_l.Length()> _otherlongtrklen )
      {
	for (size_t zz = 0; zz<vhit.size(); ++zz)
        {  
         // looping over the hits of the selected long track
          if(vhit[zz]->WireID().Plane == 2)
          {
	       trkHitsKey.push_back(vhit[zz].key());
//	       if(vhit[zz].key()==38025)cout<<"vhit[zz].ptime "<<vhit[zz]->PeakTime()<<" charge "<<vhit[zz]->Integral()<<"\n";
//	       trkHitsKey[longtrk_hits] = vhit[zz].key();
//	       longtrk_hits++;

	       
	  }
	}//loop over hits 	    
      }//if(track_l.Length()>75 )
    }//loop over tracks   
//    cout<<"trkHitsKey.size() "<<trkHitsKey.size()<<"\n";
    cout<<"trkHitsKey.size() "<<trkHitsKey.size()<<"\n";



		    //**************The real track loop starts****************//

    for(size_t i=0; i<tracklist.size();++i)
    {
      art::Ptr<recob::Track> ptrack(trackListHandle, i);
      const recob::Track& track = *ptrack;

      //std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(i);

      std::vector<art::Ptr<recob::PFParticle>> pfps=pfp_trk_assn.at(i);

      fall_trks++;
      _fall_trks++;
            
      //--------Look for the true Michel tracks------------//
      ftrueMichel1 = 0;
      float mcd_energy1=0;

      if(isMC)
      {

        // Get true MCParticle associated with recob::Track
        const simb::MCParticle *particleP1 = truthUtil.GetMCParticleFromRecoTrack(track,evt,fTrackModuleLabel);
        if(!particleP1) continue;
        const art::Ptr<simb::MCTruth> mcP1=pi_serv->TrackIdToMCTruth_P(particleP1->TrackId());
        if(!mcP1) continue;
        bool hasElectron1 = false, hasNuMu1 = false, hasNuE1 = false;
        if((abs(particleP1->PdgCode())==13) && (particleP1->NumberDaughters())>0)
	{		    
	  for (int ii=0; ii<(particleP1->NumberDaughters());++ii) 
	  {
            const simb::MCParticle* daughter1 = pi_serv->TrackIdToParticle_P((particleP1->Daughter(ii)));
            if(!daughter1) continue;
            int d_pdg1 = abs(daughter1->PdgCode());
            if (d_pdg1 == 11) hasElectron1 = true;
            else if (d_pdg1 == 14) hasNuMu1 = true;
            else if (d_pdg1 == 12) hasNuE1 = true;
	  }
	}
	if(hasElectron1 && hasNuMu1 && hasNuE1)  
	{
	  if((abs(particleP1->PdgCode())==13) && (particleP1->NumberDaughters())>0)
	  {
	    double daugE1 = 0;
		    
	    for (int ii=0; ii<(particleP1->NumberDaughters());++ii) 
	    {
              const simb::MCParticle* daughter1 = pi_serv->TrackIdToParticle_P((particleP1->Daughter(ii)));
              if(!daughter1) continue;
	      if(abs(daughter1->PdgCode())==11 && daughter1->E()>daugE1 && daughter1->Vx()>fiducialBounds1[0] && daughter1->Vx()<fiducialBounds1[1] &&
		      daughter1->Vy()>fiducialBounds1[2] && daughter1->Vy()<fiducialBounds1[3] && 
		      daughter1->Vz()>fiducialBounds1[4] && daughter1->Vz()<fiducialBounds1[5]) 
	      {     
		daugE1 = daughter1->E();
                ftrueMichel1 = 1;    
		mcd_energy1 = daughter1->E();

	      } 
             }		 		       
	   }
	 }
       } //if(isMC)  
	     
              fpfpsana[fpfpana] = pfps.size();
	      fMichelcountpfpana[fpfpana] = ftrueMichel1;
	      ftrueEpfpana[fpfpana] = mcd_energy1;
	      fpfpana++;


      //at least one PFP association
      if(pfps.size()) 
      { 
        
	fPFP_trks++;
	_fPFP_trks++;
        if(ftrueMichel1==1) _ftruePFP_trks++;


        //**************Considering only T0 tagged tracks******************//
	  	
	double T00 = 0;
	
        //---Only consider T0 tagged tracks-----//

        std::vector<art::Ptr<anab::T0>> t0s=trk_t0_assn_v.at(pfps[0].key());

        ft0sana[ft0ana] = t0s.size();
	fMichelcountt0ana[ft0ana] = ftrueMichel1;
	ftrueEt0ana[ft0ana] = mcd_energy1;
	ft0ana++;

       	if(t0s.size() )
	{ 
	  fT0_trks++;  // T0 tagged tracks per event count
         _fT0_trks++; // total T0 tagged track count
	 if(ftrueMichel1==1) _ftrueT0_trks++;
	 	  
          T00      = double(t0s.at(0)->Time()); // in nsec			

          // Store track's properties
          TVector3 dir_start = track.DirectionAtPoint<TVector3>(track.FirstValidPoint());
          TVector3 dir_end   = track.DirectionAtPoint<TVector3>(track.LastValidPoint());
          TVector3 end       = track.End<TVector3>();
          TVector3 pos(track.LocationAtPoint(track.FirstValidPoint()).X(), track.LocationAtPoint(track.FirstValidPoint()).Y(), track.LocationAtPoint(track.FirstValidPoint()).Z());
//          orderRecoStartEnd(pos, end, dir_start,dir_end);
/*        float fdcosStartX = dcosStart.X();
          float fdcosStartY = dcosStart.Y();
          float fdcosStartZ = dcosStart.Z();
          float fdcosEndX   = dcosEnd.X();
          float fdcosEndY   = dcosEnd.Y();
          float fdcosEndZ   = dcosEnd.Z();
*/        // using the ordered start and end points calculate the angles theta_xz and theta_yz
//        float ftheta_xz = TMath::RadToDeg() * TMath::ATan2(fStartX-fEndX, fStartZ-fEndZ);
//        float ftheta_yz = TMath::RadToDeg() * TMath::ATan2(fStartY-fEndY, fStartZ-fEndZ); 

//        fstartinboundsana[fstartinboundana] = t0s.size();
	fMichelcountstartinboundana[fstartinboundana] = ftrueMichel1;
	ftrueEstartinboundana[fstartinboundana] = mcd_energy1;
	fstartinboundana++;

          if(IsPointInBounds(activeBounds_eff,pos) || IsPointInBounds(activeBounds_eff,end)) // Start or end point within the edges of Active volume bounds.
	  {
	    fstartinbound_trks++;  // start in bounds tracks per event count
            _fstartinbound_trks++; // total start in bounds track count
	    if(ftrueMichel1==1) _ftruestartinbound_trks++;
	    
	    if(IsPointInBounds(activeBounds_eff,end)) 
	    {
	      SwitchEndPoints(pos, end, dir_start,dir_end);
	    }  

            float trkbegx = pos.X();
	    float trkbegy = pos.Y();
	    float trkbegz = pos.Z();
            float trkstopx = end.X();
	    float trkstopy = end.Y();
	    float trkstopz = end.Z();	      
	  
	   //-----------Checking which tracks are cathode crossers-------///
  	   int ccrosser = 0;
	   if((trkbegx*trkstopx) < 0) ccrosser = 1;
	   
	   //////////////////////////////////////////////////////////////

        fendinFVsana[fendinFVana] = IsPointInFV(fiducialBounds,end);
	fMichelcountendinFVana[fendinFVana] = ftrueMichel1;
	ftrueEendinFVana[fendinFVana] = mcd_energy1;
	fendinFVana++;

          if(IsPointInFV(fiducialBounds,end)) // End Point within the fiducial volume.
	    {
	      fendinFV_trks++;  // start in bounds tracks per event count
              _fendinFV_trks++; // total start in bounds track count

                //************************Getting truth info of the tracks to count the true michels**********************//
                ftrueMichel = 0;
	        double daugx = -999, daugy = -999, daugz = -999;
		int mcd_trkid        = -999;
                double mcd_vx         = -999;
                double mcd_vy         = -999;
                double mcd_vz         = -999;
                double mcd_t	     = -999;
                double mcd_endx       = -999;
                double mcd_endy       = -999;
                double mcd_endz       = -999;
                double mcd_endt	     = -999;
                double mcd_px         = -999;
                double  mcd_py        = -999;
                double mcd_pz	      = -999;
                double mcd_momentum   = -999;
                double mcd_energy     = -999;
                double mcd_endpx      = -999;
                double mcd_endpy      = -999;
                double mcd_endpz      = -999;
                double mcd_endenergy  = -999;
	        double mcd_pathlen      = -999;
	        double mcd_vxdrifted	     = -999;
                double mcd_vydrifted         = -999;
                double mcd_vzdrifted	    = -999;
                double mcd_tdrifted	     = -999;
                double mcd_endxdrifted       = -999;
                double mcd_endydrifted	    = -999;
                double mcd_endzdrifted       = -999;
                double mcd_endtdrifted       = -999;
                double mcd_pxdrifted         = -999;
                double mcd_pydrifted	    = -999;
                double mcd_pzdrifted         = -999;
                double mcd_momentumdrifted   = -999;
                double mcd_energydrifted     = -999;
                double mcd_endpxdrifted      = -999;
                double mcd_endpydrifted      = -999;
                double mcd_endpzdrifted      = -999;
                double mcd_endenergydrifted  = -999;
	        double mcd_pathlendrifted    = -999;
	        int mcd_endprocess = -999;
	        double mcd_theta      = -999;
	        double mcd_phi	     = -999;
                int mcd_pdg          = -999;
                int mcd_status_code  = -999;
                double mcd_mass       = -999;
                int mcd_ND         = -999;
                int mcd_mother     = -999;
	        int mcd_origin     = -999;
	        int mcd_process    = -999;
                int mcd_rescatter  = -999;
		
		cout<<"\n";
	        if(isMC)
                {
                  // Get true MCParticle associated with recob::Track
                  const simb::MCParticle *particleP = truthUtil.GetMCParticleFromRecoTrack(track,evt,fTrackModuleLabel);
                  if(!particleP) continue;
                  const art::Ptr<simb::MCTruth> mcP=pi_serv->TrackIdToMCTruth_P(particleP->TrackId());
                  if(!mcP) continue;
                  double distance = 99999/*, distance0 = 99999*/;
	     	                    
//                  unsigned int firstPointInAV = TrueParticleFirstPointInAV(fActiveBounds_eff,*particleP);
//                  unsigned int lastPointInAV = TrueParticleLastPointInAV(fActiveBounds_eff,*particleP);
	     

//                  if((abs(particleP->PdgCode())==13) && (particleP->NumberDaughters())>0)
                  bool hasElectron = false, hasNuMu = false, hasNuE = false;
                  if((abs(particleP->PdgCode())==13) && (particleP->NumberDaughters())>0)
	          {		    
	            for (int ii=0; ii<(particleP->NumberDaughters());++ii) 
	            {
                       const simb::MCParticle* daughter = pi_serv->TrackIdToParticle_P((particleP->Daughter(ii)));
                       if(!daughter) continue;
                       int d_pdg = abs(daughter->PdgCode());
                       if (d_pdg == 11) hasElectron = true;
                       else if (d_pdg == 14) hasNuMu = true;
                       else if (d_pdg == 12) hasNuE = true;
		    }
		  }
		  if(hasElectron && hasNuMu && hasNuE)  
		  {
		      
                  if((abs(particleP->PdgCode())==13) && (particleP->NumberDaughters())>0)
	          {
		    double daugE = 0;
		    
	            for (int ii=0; ii<(particleP->NumberDaughters());++ii) 
	            {
                       const simb::MCParticle* daughter = pi_serv->TrackIdToParticle_P((particleP->Daughter(ii)));
                       if(!daughter) continue;
	                      TLorentzVector dmcstart, dmcend, dmcstartdrifted, dmcenddrifted;
	                      unsigned int dpstarti, dpendi, dpstartdriftedi, dpenddriftedi; //mcparticle indices for starts and ends in tpc or drifted volumes
	                      double dplen = length(*daughter, dmcstart, dmcend, dpstarti, dpendi);
	                      double dplendrifted = driftedLength(*daughter, dmcstartdrifted, dmcenddrifted, dpstartdriftedi, dpenddriftedi);
	                      // bool isActive = plen != 0;
	                      bool disDrifted = dplendrifted!= 0;

	                      distance=sqrt(pow(daughter->Vx()-trkstopx,2)+pow(daughter->Vz()-trkstopz,2)); //only in collection plane

		      cout<<"mother pdg "<<particleP->PdgCode()<<" mother energy "<<particleP->E()<<" true end x, y, z "<<particleP->EndX()<<" "<<particleP->EndY()<<" "<<particleP->EndZ()<<" trackid "<<particleP->TrackId()<<"\n";
		      cout<<"daughter start distance "<<distance<<" x, y, z "<<daughter->Vx()<<" "<<daughter->Vy()<<" "<<daughter->Vz()<<" end x, y, z "<<daughter->EndX()<<" "<<daughter->EndY()<<" "<<daughter->EndZ()<<"\n";
		      cout<<"daughter->PdgCode() "<<daughter->PdgCode()<<" daughter energy "<<daughter->E()<<" daughter process "<<daughter->Process()<<" daughter track id "<<daughter->TrackId()<<"\n";


		      
//		      cout<<"mother pdg "<<particleP->PdgCode()<<", mother particle ID "<<particleP->TrackId()<<", daughter pdg "<<daughter->PdgCode()<<", daughter particle ID  "<<daughter->TrackId()<<", daughter true energy "<<daughter->E()<<"\n";


	              if(distance<=30 && abs(daughter->PdgCode())==11 && daughter->E()>daugE && daughter->Vx()>fiducialBounds1[0] && daughter->Vx()<fiducialBounds1[1] &&
		      daughter->Vy()>fiducialBounds1[2] && daughter->Vy()<fiducialBounds1[3] && 
		      daughter->Vz()>fiducialBounds1[4] && daughter->Vz()<fiducialBounds1[5]) // z positions are slightly different between a true and reco particles 
		      {     
           
                        fdistanceana[fdistana] = distance;
	                fdistana++;
			
//			  distance0 = distance;
			  
			  daugE = daughter->E();

		          ftrueMichel   = 1;
		          daugx = daughter->Vx();
		          daugy = daughter->Vy();
		          daugz = daughter->Vz();
//		        }
//	                if(ftrueMichel == 1)
//		        {
//		          _ftrueMichelcount++; 
//		          break;	
//		        }	     
                    
		              mcd_trkid      = daughter->TrackId();
                              mcd_vx         = daughter->Vx();
                              mcd_vy         = daughter->Vy();
                              mcd_vz         = daughter->Vz();
                              mcd_t          = daughter->T();
                              mcd_endx       = daughter->EndX();
                              mcd_endy       = daughter->EndY();
                              mcd_endz       = daughter->EndZ();
                              mcd_endt	     = daughter->EndT();
                              mcd_px         = daughter->Px();
                              mcd_py         = daughter->Py();
                              mcd_pz         = daughter->Pz();
                              mcd_momentum   = daughter->P();
                              mcd_energy     = daughter->E();
                              mcd_endpx      = daughter->EndPx();
                              mcd_endpy      = daughter->EndPy();
                              mcd_endpz      = daughter->EndPz();
                              mcd_endenergy  = daughter->EndE();
	                      mcd_pathlen       = dplen;
	                      if(disDrifted)
	                      {
                                mcd_vxdrifted         = daughter->Vx();
                                mcd_vydrifted         = daughter->Vy();
                                mcd_vzdrifted         = daughter->Vz();
                                mcd_tdrifted          = daughter->T();
                                mcd_endxdrifted	      = daughter->EndX();
                                mcd_endydrifted	      = daughter->EndY();
                                mcd_endzdrifted       = daughter->EndZ();
                                mcd_endtdrifted       = daughter->EndT();
                                mcd_pxdrifted         = daughter->Px();
                                mcd_pydrifted	      = daughter->Py();
                                mcd_pzdrifted         = daughter->Pz();
                                mcd_momentumdrifted   = daughter->P();
                                mcd_energydrifted     = daughter->E();
                                mcd_endpxdrifted      = daughter->EndPx();
                                mcd_endpydrifted      = daughter->EndPy();
                                mcd_endpzdrifted      = daughter->EndPz();
                                mcd_endenergydrifted  = daughter->EndE();
	                        mcd_pathlendrifted    = dplendrifted;
	                      }
	                      mcd_endprocess = int(daughter->Trajectory().ProcessToKey(daughter->EndProcess()));
	                      mcd_theta      = daughter->Momentum().Theta();
	                      mcd_phi        = daughter->Momentum().Phi();
                              mcd_pdg        = daughter->PdgCode();
                              mcd_status_code= daughter->StatusCode();
                              mcd_mass       = daughter->Mass();
                              mcd_ND         = daughter->NumberDaughters();
                              mcd_mother     = daughter->Mother();
	                      mcd_origin     = mcP->Origin();
	                      mcd_process    = int(daughter->Trajectory().ProcessToKey(daughter->Process()));
                              mcd_rescatter  = daughter->Rescatter();
		     
		              cout<<"daughter->PdgCode() "<<daughter->PdgCode()<<" daughter track id "<<daughter->TrackId()<<"\n";
//	                      cout<<" daughter.Mother() "<<daughter->Mother()<<" vx, vy, vz "<<fmcd_vx<<" "<<fmcd_vy<<" "<<fmcd_vz<<" endx, endy, endz "<<fmcd_endx<<" "<<fmcd_endy<<" "<<fmcd_endz<<"\n";

		      } 
                    }		 		       
	          }
	         }

    	        } //if(isMC)  
	     
		if(ftrueMichel==1) _ftrueendinFV_trks++;

              //***************Consider only cathode crossers******************//

              fccrosserana[fccrossana] = ccrosser;
	      fMichelcountccrosserana[fccrossana] = ftrueMichel;
	      ftrueEccrosserana[fccrossana] = mcd_energy;
	      fccrossana++;
	      
	      if(ccrosser == 1)
	      {

	        fccrosser_trks++;
	        _fccrosser_trks++;
		if(ftrueMichel==1) _ftrueccrosser_trks++;


              //***************Removing electron diverters boundaries******************//
	  
              felecdivstopzana[fstopzana] = trkstopz;
	      fMichelcountelecdivstopzana[fstopzana] = ftrueMichel;
	      ftrueEelecdivstopzana[fstopzana] = mcd_energy;
	      fstopzana++;

	      if((trkstopz<APAnegbound1 || trkstopz>APAposbound1) && (trkstopz<APAnegbound2 || trkstopz>APAposbound2))
	      {
	      
	        fnoelecdiv_bound++;
	        _fnoelecdiv_bound++;
		if(ftrueMichel==1) _ftruenoelecdiv_bound++;


	     /**********************************Broken tracks removal****************************************/
             size_t bcount = 0;
  	     for(size_t k=0;k<tracklist.size();++k)
	     {
	       art::Ptr<recob::Track> ptrack_k(trackListHandle, k);
	       const recob::Track& track_k = *ptrack_k;

	       TVector3 dir_start_k, dir_end_k;
               TVector3 pos_k(track_k.LocationAtPoint(track_k.FirstValidPoint()).X(), track_k.LocationAtPoint(track_k.FirstValidPoint()).Y(), track_k.LocationAtPoint(track_k.FirstValidPoint()).Z());
               TVector3 end_k = track_k.End<TVector3>();
               dir_start_k= track_k.DirectionAtPoint<TVector3>(track_k.FirstValidPoint());
               dir_end_k   = track_k.DirectionAtPoint<TVector3>(track_k.LastValidPoint());

               if(k==i) continue;
	 	       
	       double bwdiststart          = sqrt(pow(pos_k.Y()-end.Y(),2)+pow(pos_k.Z()-end.Z(),2)); //since the non-T0 tagged tracks have wrong x position associated with them
               double cosopeningangleStart = cos(dir_end.Theta())*cos(dir_start_k.Theta())+sin(dir_end.Theta())*sin(dir_start_k.Theta())*cos(dir_start_k.Phi()-dir_end.Phi());
	       double bwdistend            = sqrt(pow(end_k.Y()-end.Y(),2)+pow(end_k.Z()-end.Z(),2));
               double cosopeningangleEnd   = cos(dir_end.Theta())*cos(dir_end_k.Theta())+sin(dir_end.Theta())*sin(dir_end_k.Theta())*cos(dir_end_k.Phi()-dir_end.Phi());

	       if((abs(bwdiststart)<30 && abs(cosopeningangleStart) > 0.97) || (abs(bwdistend)<30 && abs(cosopeningangleEnd) > 0.97)
	        ||(abs(bwdiststart)<50 && abs(cosopeningangleStart) > 0.998) || (abs(bwdistend)<50 && abs(cosopeningangleEnd) > 0.998))
	       {
	         bcount++;
	       } 
	       
/*               if((std::abs(((end_k.Y()-pos_k.Y())/(end_k.Z()-pos_k.Z()))*(trkstopz-pos_k.Z())+pos_k.Y()-trkstopy)<30) &&
                  (std::abs(dir_end.X()*dir_pos_k.X()+dir_end.Y()*dir_pos_k.Y()+dir_end.Z()*dir_pos_k.Z())>0.97 ||
       	           std::abs(dir_end.X()*dir_end_k.X()+dir_end.Y()*dir_end_k.Y()+dir_end.Z()*dir_end_k.Z())>0.97 )) break;
		
	       if((std::abs(((end_k.Y()-pos_k.Y())/(end_k.Z()-pos_k.Z()))*(trkstopz-pos_k.Z())+pos_k.Y()-trkstopy)<50) &&
	          (std::abs(dir_end.X()*dir_pos_k.X()+dir_end.Y()*dir_pos_k.Y()+dir_end.Z()*dir_pos_k.Z())>0.998 ||
       	           std::abs(dir_end.X()*dir_end_k.X()+dir_end.Y()*dir_end_k.Y()+dir_end.Z()*dir_end_k.Z())>0.998)) break;
		
	        count++;
*/	    }	 
      
//              fbrokencountana[fbrokcoana] = tracklist.size()-count;
              fbrokencountana[fbrokcoana] = bcount;
	      fMichelcountbrokencountana[fbrokcoana] = ftrueMichel;
	      ftrueEbrokencountana[fbrokcoana] = mcd_energy;
	      fbrokcoana++;

//	      if(count==tracklist.size()-1) //unbroken tracks
	      if(bcount==0) //unbroken tracks
	      {
                funbroken_trks++;
	        _funbroken_trks++;
	        if(ftrueMichel==1)_ftrueunbroken_trks++;
		

                //------------------track length cut------------------// 
                ftrklengthana[ftrklenana] = track.Length();
	        fMichelcountlenana[ftrklenana] = ftrueMichel;
		ftrueElenana[ftrklenana] = mcd_energy;
	        ftrklenana++;
 
 	        if(track.Length()> _longtrklen ) //the selected track has to have at least X cm long
	        {
	          flongcm_trks++;
	          _flongcm_trks++;
	          if(ftrueMichel==1)_ftruelongcm_trks++;

                  //----------------Get Hits associated with track-----------//
		  int ncolhits = 0;
                  const std::vector<const recob::Hit*> Hits = trackUtil.GetRecoTrackHits(track,evt,fTrackModuleLabel);
                  // Fill vector with Hit Peak Times and store minimum
                  std::vector<double> HitPeakTimes;
                  for (unsigned int hitIndex = 0; hitIndex < Hits.size(); ++hitIndex)  
                  {
                    HitPeakTimes.push_back(Hits[hitIndex]->PeakTime());
                    if(Hits[hitIndex]->WireID().Plane==2) ncolhits++;      
                  }
                  float MinHitPeakTime = *(std::min_element(HitPeakTimes.begin(), HitPeakTimes.end()));
                  float MaxHitPeakTime = *(std::max_element(HitPeakTimes.begin(), HitPeakTimes.end()));
		  
		  HitPeakTimes.clear();
    
                  //--------------Min hit peak time cut--------------------// 
                  fminhitptimeana[fminhitptana] = MinHitPeakTime;
	          fMichelcountminhitptimeana[fminhitptana] = ftrueMichel;
		  ftrueEminhitptimeana[fminhitptana] = mcd_energy;
	          fminhitptana++;
		  
		  if(MinHitPeakTime > _minhitpeakcut)
		  {
	            fminhitpt_trks++;
	            _fminhitpt_trks++;
	            if(ftrueMichel==1)_ftrueminhitpt_trks++;

		    
		    //-----------Max hit peak time cut------------------//
                    fmaxhitptimeana[fmaxhitptana] = MaxHitPeakTime;
	            fMichelcountmaxhitptimeana[fmaxhitptana] = ftrueMichel;
		    ftrueEmaxhitptimeana[fmaxhitptana] = mcd_energy;
	            fmaxhitptana++;
		  
		    if(MaxHitPeakTime < _maxhitpeakcut)
		    {
	              fmaxhitpt_trks++;
	              _fmaxhitpt_trks++;
	              if(ftrueMichel==1)_ftruemaxhitpt_trks++;

/*
                    //--------------collection plane hits cut-------------//
	            ftrkcollhitsana[ftrkcolana] = ncolhits;
	            fMichelcountcollana[ftrkcolana] = ftrueMichel;
                    ftrkcolana++;

	            if(ncolhits >= _longtrkcollhits) // the selected track has to have at least Y coll plane hits
	            {
                      flonghits_trks++;
	              _flonghits_trks++;
	              if(ftrueMichel==1)_ftruelonghits_trks++;
		      
*/

		    
/*		    //-----------Muon end y cut------------------//
                    fmuonendyana[fmuendyana] = trkstopy;
	            fMichelcountmuonendyana[fmuendyana] = ftrueMichel;
	            fmuendyana++;
		  
		    if(trkstopy > _muendycut)
		    {
	              fmuendy_trks++;
	              _fmuendy_trks++;
	              if(ftrueMichel==1)_ftruemuendy_trks++;
		      
		    //-----------Muon end z cut------------------//
                    fmuonendzana[fmuendzana] = trkstopz;
	            fMichelcountmuonendzana[fmuendzana] = ftrueMichel;
	            fmuendzana++;
		  
		    if(trkstopz > _muendzmincut && trkstopz < _muendzmaxcut)
		    {
	              fmuendz_trks++;
	              _fmuendz_trks++;
	              if(ftrueMichel==1)_ftruemuendz_trks++;
*/		      
                    cout<<"ftrueMichel "<<ftrueMichel<<" daugx "<<daugx<<" daugy "<<daugy<<" daugz "<<daugz<<"\n";  


	          //*****************Storing the candidate muon track collection plane hits****************************//
    
	          float  dist=99999;
	          float  enddist=99999;
	          float  endhitkey=-999;
	          float  endpeaktime=-999;
	          float  endwireno=-999;
	          float  endtpcno=-999;
	          float  endchno=-999;
	          float  endhitchrg=-999;
	          float  endhitx=-999;
	          float  endhity=-999;
	          float  endhitz=-999;
	          float  endhitmult=-999;
	          float  endhitsigptime=-999;
	          float  endhitsigchrg=-999;
	          float  endhitsigpamp=-999;
	          float  endhitdof=-999;
	          float  endhitgof=-999;
	          float  endhitptminusRMS=-999;
	          float  endhitptplusRMS=-999;
 
                  double hits_key[kMaxCh]       ={-999};
                  double hits_charge[kMaxCh]    ={-999};
	          double hits_chno[kMaxCh]      ={-999};
                  double hits_wire[kMaxCh]      ={-999};
                  double hits_peakT[kMaxCh]     ={-999};
	          double hits_TPC[kMaxCh]       ={-999};
	          float hits_xpos[kMaxCh]       ={-999};
	          float hits_ypos[kMaxCh]       ={-999};
	          float hits_zpos[kMaxCh]       ={-999};
	          float hits_mult[kMaxCh]       ={-999};
	          float hits_sigptime[kMaxCh]   ={-999};
	          float hits_sigchrg[kMaxCh]    ={-999};
	          float hits_sigpamp[kMaxCh]     ={-999};
	          float hits_dof[kMaxCh]        ={-999};
	          float hits_gof[kMaxCh]        ={-999};
	          float hits_ptminusRMS[kMaxCh] ={-999};
	          float hits_ptplusRMS[kMaxCh]  ={-999};
	          float hits_cnnMichel[kMaxCh]  ={-999};
	          float hits_cnnEM[kMaxCh]      ={-999};
	          float hits_cnnTrack[kMaxCh]        ={-999};
 	          //*************Determine the position of the last hit and storing the collection plane hits for each track****************************//
		  
                  int trkcolhits=0;
		  std::vector <double> longtrk_hitkey;
		  
	          if(fmthm.isValid())
	          {
	            auto vhit=fmthm.at(i);
	            auto vmeta=fmthm.data(i);
		
      
	            //loop over all meta data hit
	            for (size_t ii = 0; ii<vhit.size(); ++ii)
	            {
		
	              bool fBadhit = false;
	              if (vmeta[ii]->Index() == std::numeric_limits<int>::max())
	              {
	                fBadhit = true;
	                //continue;
	              }
	              if (vmeta[ii]->Index()>=tracklist[i]->NumberTrajectoryPoints())
	              {
	                fBadhit = true;
		        //throw cet::exception("Calorimetry_module.cc") << "Requested track trajectory index "<<vmeta[ii]->Index()<<" exceeds the total number of trajectory points "<<tracklist[i]->NumberTrajectoryPoints()<<" for track index "<<i<<". Something is wrong with the track reconstruction. Please contact tjyang@fnal.gov!!";
	              }
	              if (!tracklist[i]->HasValidPoint(vmeta[ii]->Index()))
	              {
	                fBadhit = true;
	                //continue;
	              }

	    	      TVector3 loc;	  

	              if (fBadhit) continue;
		      else loc = tracklist[i]->LocationAtPoint<TVector3>(vmeta[ii]->Index());
		  
 	    	      if(vhit[ii]->WireID().Plane==2)
	              {
			
		        longtrk_hitkey.push_back(vhit[ii].key());

			if(loc.Z()!=0 && (loc.Z()<APAnegbound1 || loc.Z()>APAposbound1) && (loc.Z()<APAnegbound2 || loc.Z()>APAposbound2))
			{
	                  hits_key[trkcolhits]       =   vhit[ii].key();
	                  hits_charge[trkcolhits]    =   vhit[ii]->Integral();///(TMath::Exp(-(hits_peakT[trkcolhits1]-800)*500/tau)); //multiplied by 500nsec to convert time ticks to actual generation time
                          hits_wire[trkcolhits]      =   vhit[ii]->WireID().Wire;
	                  hits_peakT[trkcolhits]     =   vhit[ii]->PeakTime();
	                  hits_TPC[trkcolhits]       =   vhit[ii]->WireID().TPC;
		          hits_chno[trkcolhits]      =   vhit[ii]->Channel();
	                  hits_xpos[trkcolhits]      =   loc.X();
	                  hits_ypos[trkcolhits]      =   loc.Y();
	                  hits_zpos[trkcolhits]      =   loc.Z();
	                  hits_mult[trkcolhits]      =   vhit[ii]->Multiplicity();
	                  hits_sigptime[trkcolhits]  =   vhit[ii]->SigmaPeakTime();
	                  hits_sigchrg[trkcolhits]   =   vhit[ii]->SigmaIntegral();
	                  hits_sigpamp[trkcolhits]   =   vhit[ii]->SigmaPeakAmplitude();
	                  hits_dof[trkcolhits]       =   vhit[ii]->DegreesOfFreedom();
	                  hits_gof[trkcolhits]       =   vhit[ii]->GoodnessOfFit();      
		          hits_ptminusRMS[trkcolhits]=   vhit[ii]->PeakTimeMinusRMS(5.0);
		          hits_ptplusRMS[trkcolhits] =   vhit[ii]->PeakTimePlusRMS(5.0);

                          std::array<float,4> cnn_out = hitResults.getOutput(vhit[ii].key());
//			  cout<<"ii "<<ii<<" key "<<vhit[ii].key()<<"\n";
//                          double p_trk_or_sh = cnn_out[ hitResults.getIndex("track") ]+ cnn_out[ hitResults.getIndex("em") ]+ cnn_out[ hitResults.getIndex("michel") ]; 
//                          double cnn_score = cnn_out[ hitResults.getIndex("michel") ]; 		           
			  hits_cnnMichel[trkcolhits] =   cnn_out[ hitResults.getIndex("michel") ]; 
			  hits_cnnEM[trkcolhits]     =   cnn_out[ hitResults.getIndex("em") ]; 
			  hits_cnnTrack[trkcolhits]  =   cnn_out[ hitResults.getIndex("track") ]; 
     		    
	                  dist=sqrt(pow(hits_xpos[trkcolhits]-trkstopx,2)+pow(hits_ypos[trkcolhits]-trkstopy,2)+pow(hits_zpos[trkcolhits]-trkstopz,2));
 
	                  //cout<<"dist "<<dist<<" Hit num "<<trkcolhits1<<" Hit key "<<hits_key[trkcolhits1]<<" hit meta: channel ID "<<hits_channel[trkcolhits1]<<" hits_wire "<<hits_wire[trkcolhits]<<" hits_TPC "<<hits_TPC[trkcolhits]<<"  hits_peakT[trkcolhits] "<<hits_peakT[trkcolhits]<<"\n";
	                  //std::cout<<"hits_channel[trkcolhits] "<<hits_channel[trkcolhits]<<" hits_wire[trkcolhits] "<<hits_wire[trkcolhits]<<" hits_TPC[trkcolhits] "<<hits_TPC[trkcolhits]<<"\n";

	                  trkcolhits++;
					    
		          //---------storing the information for the end point of the track----------//

	                  if(abs(dist)<enddist)
	                  {
		            enddist=dist;
		            endhitkey=vhit[ii].key();
		            endwireno=vhit[ii]->WireID().Wire;
		            endpeaktime=vhit[ii]->PeakTime();
		            endtpcno=vhit[ii]->WireID().TPC;
		            endchno=vhit[ii]->Channel();
		            endhitchrg=vhit[ii]->Integral();
		            endhitx=loc.X();
		            endhity=loc.Y();
		            endhitz=loc.Z();
	                    endhitmult=vhit[ii]->Multiplicity();
	                    endhitsigptime=vhit[ii]->SigmaPeakTime();
	                    endhitsigchrg=vhit[ii]->SigmaIntegral();
	                    endhitsigpamp=vhit[ii]->SigmaPeakAmplitude();
	                    endhitdof=vhit[ii]->DegreesOfFreedom();
	                    endhitgof=vhit[ii]->GoodnessOfFit();
	                    endhitptminusRMS=vhit[ii]->PeakTimeMinusRMS(5.0);
	                    endhitptplusRMS=vhit[ii]->PeakTimePlusRMS(5.0);
			  }  
	                }	
	              }//if(vhit[ii]->WireID().Plane==2)
	            }//loop over vhit
	          }//fmthm valid
                   
		   
		   //-------Calculate the 
                  //**************Rearrange the wires such that the track end point is at the end of the wire order**********//
	      	      
                  //------find the first (and last) non-zero distances from the end wire location------//
              
	          double dist0= 999, dist1 = 999;
	      
	          // first non-zero hit distance
	          for(int jj = 0; jj<trkcolhits; jj++)
	          {
		    dist0 = sqrt(pow(hits_xpos[jj]-trkstopx,2)+pow(hits_ypos[jj]-trkstopy,2)+pow(hits_zpos[jj]-trkstopz,2));
		    break;  
	          }
	          // last non-zero hit distance
	          for(int jj = trkcolhits-1; jj>=0; jj--)
	          {
		    dist1 = sqrt(pow(hits_xpos[jj]-trkstopx,2)+pow(hits_ypos[jj]-trkstopy,2)+pow(hits_zpos[jj]-trkstopz,2));
		    break;
	          }
	      
                  //------Rearranging the wires-----------//
	          int a, c, e, l;
	          double b, d, f, g ,h,r, m, n, o, p, q, s,t,u,v,w;
	          //if(abs(dist1-enddist) > abs(dist0-enddist) )	
	          if(abs(dist1) > abs(dist0) )	
	          {
	            for(int jj=0; jj<trkcolhits-1; jj++)
		    {
                      for (int kk = jj + 1; kk < trkcolhits; kk++)
                      { 
		        a = hits_key[jj];
		        hits_key[jj] = hits_key[kk];
		        hits_key[kk] = a;
		    
		        b = hits_charge[jj];
		        hits_charge[jj] = hits_charge[kk];
		        hits_charge[kk] = b;

		        c = hits_wire[jj];
		        hits_wire[jj] = hits_wire[kk];
		        hits_wire[kk] = c;

		        d = hits_peakT[jj];
		        hits_peakT[jj] = hits_peakT[kk];
		        hits_peakT[kk] = d;

		        e = hits_TPC[jj];
		        hits_TPC[jj] = hits_TPC[kk];
		        hits_TPC[kk] = e;
		    
		        f = hits_chno[jj];
		        hits_chno[jj] = hits_chno[kk];
		        hits_chno[kk] = f;

		        g = hits_xpos[jj];
		        hits_xpos[jj] = hits_xpos[kk];
		        hits_xpos[kk] = g;
		    
		        h = hits_ypos[jj];
		        hits_ypos[jj] = hits_ypos[kk];
		        hits_ypos[kk] = h;
		    
		        r = hits_zpos[jj];
		        hits_zpos[jj] = hits_zpos[kk];
		        hits_zpos[kk] = r;

	                l = hits_mult[jj];
		        hits_mult[jj] = hits_mult[kk];
		        hits_mult[kk] = l;
			
	                m = hits_sigptime[jj];
		        hits_sigptime[jj] = hits_sigptime[kk];
		        hits_sigptime[kk] = m;

	                n = hits_sigchrg[jj];
		        hits_sigchrg[jj] = hits_sigchrg[kk];
		        hits_sigchrg[kk] = n;

	                o = hits_sigpamp[jj];
		        hits_sigpamp[jj] = hits_sigpamp[kk];
		        hits_sigpamp[kk] = o;

	                p = hits_dof[jj];
		        hits_dof[jj] = hits_dof[kk];
		        hits_dof[kk] = p;

	                q = hits_gof[jj];
		        hits_gof[jj] = hits_gof[kk];
		        hits_gof[kk] = q;

	                s = hits_ptminusRMS[jj];
		        hits_ptminusRMS[jj] = hits_ptminusRMS[kk];
		        hits_ptminusRMS[kk] = s;

	                t = hits_ptplusRMS[jj];
		        hits_ptplusRMS[jj] = hits_ptplusRMS[kk];
		        hits_ptplusRMS[kk] = t;

	                u = hits_cnnMichel[jj];
		        hits_cnnMichel[jj] = hits_cnnMichel[kk];
		        hits_cnnMichel[kk] = u;

	                v = hits_cnnEM[jj];
		        hits_cnnEM[jj] = hits_cnnEM[kk];
		        hits_cnnEM[kk] = v;

	                w = hits_cnnTrack[jj];
		        hits_cnnTrack[jj] = hits_cnnTrack[kk];
		        hits_cnnTrack[kk] = w;
		      }
		    }    		    
 	          }	
 	      
/*                  //---------------Excluding hits that may belong to a michel------------------------//
         	  double mean_charge = 0;
	          double trunc_charge[10000] = {-999};
	          float max_truncchrg = 0;
		  int trkcolhits = 0;
		  
                  for (int kk=1; kk < trkcolhits1; kk++)
                  { 
       	            if(kk<trkcolhits1-1)
	            mean_charge = (hits_charge[kk-1]+hits_charge[kk]+hits_charge[kk+1])/ 3;
	            if(kk==trkcolhits1-1)
	            mean_charge = (hits_charge[kk-1]+hits_charge[kk])/ 2;
	     
	            if(hits_charge[kk]>0.2*mean_charge && hits_charge[kk]<2*mean_charge)
	            trunc_charge[kk] = mean_charge;
	            else
	            trunc_charge[kk] = hits_charge[kk];

                    if(kk>=trkcolhits1-10 && kk<trkcolhits1)
	            {
                      if(trunc_charge[kk]>max_truncchrg)
	              {
	                max_truncchrg = trunc_charge[kk];
		        trkcolhits = kk;
		      }  
		    }
                  }
		  
		  //-------------storing long track hit keys------------------//
		  std::vector <double> longtrk_hitkey;
		  
                  for (int kk=0; kk < trkcolhits; kk++)
                  { 
		    longtrk_hitkey.push_back(hits_key[kk]);
		  }  
 		  
*/		  
/*                  //--------------------------------PH Test----------------------------------------//		
                    int PH =-999;			
                    double uphitcharge=0, downhitcharge=0;
                    //int upwire[20]={-999}, downwire[20]={-999};
                    double upchrg[40]={-999}, downchrg[40]={-999};
                    double avgupchrg=-999, avgdownchrg=-999;
                    int kk=0, ll=0, prevk=-999, prevl=-999;
                    double prevavgupchrg=-999, prevavgdownchrg=-999;
    	 
                    //upstream hits			 								
                    for(int i1=0; i1 < trkcolhits;i1++)
                    {
 	              if(i1>=trkcolhits-100 && i1<trkcolhits-60)    //Add up corrected charge on 40 hits
	              {		    
	                uphitcharge   += hits_charge[i1];     //adding charge on vertex side 20 hits of the track	
		        upchrg[kk]    = hits_charge[i1];
	                //std::cout<<"hit_wire: "<<hits_wire[i1]<<" "<<"hit_charge: "<<hits_charge[i1]<<"\n"; 
	                //upwire[kk]  = hits_wire[i1];
	                //std::cout<<kk<<" hit_charge: "<<upchrg[kk]<<"\n";
	                kk++; 
	              }
	            }
                    //-----------------get Truncated mean charge--------------------//
                    //upstram hits
                    avgupchrg  =  uphitcharge/kk;
                    //std::cout<<"avgupchrg: "<<avgupchrg<<"\n";
 
                    do
                    {   
                      prevk= kk;
	              std::cout<<"avg.up charge: "<<avgupchrg<<"prev k: "<<prevk<<"\n";
	              kk=0; 
	              uphitcharge = 0;
                      for(int i1=0; i1<prevk; i1++)
	              {
                        if(upchrg[i1]>0.2*avgupchrg && upchrg[i1]<5*avgupchrg)
	                {
                          uphitcharge += upchrg[i1];
                          upchrg[kk]    = upchrg[i1];
                          kk++;
       	                }
                      }
	              if(kk!=0)
	              {
	                prevavgupchrg =  avgupchrg;
	                avgupchrg  =  uphitcharge/kk;
	                std::cout<<"prev avg up charge: "<<prevavgupchrg<<" avg up charge: "<<avgupchrg<<"\n";	
                      }
	              else break;
	            }			
                    while(avgupchrg!=prevavgupchrg);


                    //downstream hits	
                    for(int i1=0; i1 < trkcolhits;i1++)
                    {
  	              if(i1<trkcolhits && i1>=trkcolhits-40)
	              {
	                downhitcharge += hits_charge[i1];   //adding charge on last 20 hits of the track	
			downchrg[ll]    = hits_charge[i1];
	                //std::cout<<i<<" hit_wire: "<<hits_wire[i1]<<" "<<"hit_charge: "<<hits_charge[i1]<<"\n";
	                //downwire[ll]    = hits_wire[i1];
	                //std::cout<<ll<<" hit_charge: "<<downchrg[ll]<<"\n";	 
	                ll++;                       
	              }    
                    }
                    //get Truncated mean
                    avgdownchrg  =  downhitcharge/ll;
                    do
                    {   
                      prevl= ll;
	              std::cout<<"avg down charge: "<<avgdownchrg<<"prev l: "<<prevl<<"\n";
	              ll=0; 
	              downhitcharge =0;
                      for(int i1=0; i1<prevl; i1++)
	              {
                        if(downchrg[i1]>0.2*avgdownchrg && downchrg[i1]<5*avgdownchrg)
	                {
                          downhitcharge += downchrg[i1];
                          downchrg[ll]    = downchrg[i1];
                          ll++;
       	                }
                      }
	              if(ll!=0)
	              {
	                prevavgdownchrg =  avgdownchrg;
	                avgdownchrg  =  downhitcharge/ll;
	                std::cout<<"prev avg down charge: "<<prevavgdownchrg<<" avg down charge: "<<avgdownchrg<<"\n";		
                      }
	              else break;
	            }			
                    while(avgdownchrg!=prevavgdownchrg);
   			  

                    if((avgdownchrg > avgupchrg) && avgdownchrg!=-999 && avgupchrg!=-999)	
                    {
                      PH=1;
                    }	
                    else if ((avgdownchrg <= avgupchrg) && avgdownchrg!=-999 && avgupchrg!=-999)
                    {
                      PH=-1;
                    }
	            else
	            {
	              PH=0;
	            }   
                    std::cout<<"PH "<<PH<<"\n";   
	    
                    fPHtestana[fPHana] = PH;
	            fMichelcountPHtestana[fPHana] = ftrueMichel;
	            fPHana++;
		    //cout<<"1st: ftrueMichel "<<ftrueMichel<<"\n";

	            if(PH==_PHpass)
	            {
	              _fPHtest_trks++;
	              if(ftrueMichel==1)_ftruePHtest_trks++;
		      
		      fPHtest_trks++;
	      
*/	      	        	

/*            //-----------Determine the number of distinguished coll plane hits-----------//
	            int resl = 1;
	            for(int kk=1; kk<trkcolhits; kk++)
	            { 
                      // Count distinguish elements one by one 
         
	              int ll = 0;
                      for (ll = 0; ll < kk; ll++) 
                      if (hits_wire[kk] == hits_wire[ll]) 
                      break; 
  
                      // If not printed earlier, then print it 
                      if (kk == ll) 
                      resl++; 
                    } 
                    cout<<"distinguish numbers of coll plane hits "<<resl<<"\n"; 	 
                    cout<<"ncolhits "<<ncolhits<<" trkcolhits "<<trkcolhits<<"\n";
		    
                    //--------------Distinguished collection plane hits cut-------------//
		    
	            ftrkdistcollhitsana[ftrkdistcolana] = float(resl)/float(trkcolhits);
	            fMichelcountdistcollana[ftrkdistcolana] = ftrueMichel;
                    ftrkdistcolana++;
		    cout<<"float(resl/trkcolhits) "<<float(resl)/float(trkcolhits)<<"\n";

	            if(float(resl)/float(trkcolhits) >= 0.5) // the selected track has to have at least 50% distinguish hits
	            {
                      fdistmorehits_trks++;
	              _fdistmorehits_trks++;
	              if(ftrueMichel==1)_ftruedistmorehits_trks++;
		      
*/

	               //--------Finding the hits close to the stopping point of the tracks-------//
	               //bool foundMichel = false;
                       int nearhitct=0;
	               int nearhits_key[kMaxHits]       = {-999};
		       float nearhits_peakT[kMaxHits]   = {-999};
		       float nearhits_charge[kMaxHits]  = {-999};
		       float nearhits_wire[kMaxHits]    = {-999};
		       float nearhits_chno[kMaxHits]    = {-999};
		       float nearhits_TPC[kMaxHits]     = {-999};
		       float nearhits_plane[kMaxHits]   = {-999};
		       float nearhits_xpos[kMaxHits]    = {-999};
		       float nearhits_ypos[kMaxHits]    = {-999};
		       float nearhits_zpos[kMaxHits]    = {-999};
		       float nearhits_mult[kMaxHits]    = {-999};
		       float nearhits_sigptime[kMaxHits]= {-999};
		       float nearhits_sigchrg[kMaxHits] = {-999};
		       float nearhits_sigpamp[kMaxHits]  = {-999};
		       float nearhits_dof[kMaxHits]     = {-999};
		       float nearhits_gof[kMaxHits]     = {-999};
		       float nearhits_ptminusRMS[kMaxHits]    = {-999};
		       float nearhits_ptplusRMS[kMaxHits]     = {-999};
		       float nearhits_cnnMichel[kMaxHits]     = {-999};
		       float nearhits_cnnEM[kMaxHits]         = {-999};
		       float nearhits_cnnTrack[kMaxHits]      = {-999};


		       //cout<<"trkcolhits "<<trkcolhits<<" longtrk_hitkey.size() "<<longtrk_hitkey.size()<<"\n";

                       auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
	               art::ServiceHandle<geo::Geometry const> geometry;

                       for(size_t ll=0; ll<hitlist.size();++ll) //loop over all hits
                       {
                         //auto tracks = thass.at(hitlist[ll].key());
			 
		 
	                 if(hitlist[ll]->WireID().Plane==2)
	                 {
	                   //----Make sure that this hit does not belong to the candidate muon track----//
			   
//	       if(hitlist[ll].key()==38025)cout<<"hitlist[ll].ptime "<<hitlist[ll]->PeakTime()<<" charge "<<hitlist[ll]->Integral()<<"\n";
	                   int same_trk = 0;
	                   for(size_t mm=0;mm<longtrk_hitkey.size();mm++)
	                   {
			     if(longtrk_hitkey.at(mm)==hitlist[ll].key())
		             { 
                               same_trk = 1;
		             }
                           }
		
                           //---Make sure that these hits don't belong to another long track-----//
				  
        	           int long_trk = 0;
		           for(size_t o1=0;o1<trkHitsKey.size();o1++)
		           {
		             if(trkHitsKey.at(o1)==hitlist[ll].key())
		             {
		               long_trk = 1;
		             }
		           }
/*
        	           int long_trk = 0;
		           for(Int_t o1=0;o1<longtrk_hits;o1++)
		           {
		             if(trkHitsKey[o1]==hitlist[ll].key())
		             {
		               long_trk = 1;
		             }
		           }
*/
///                          std::vector<art::Ptr<recob::SpacePoint>> sppt=fmsp.at(ll);
			  
			  
                           double diffpeaktt0     = hitlist[ll]->PeakTime() - (T00/1000)*2; //calculate the peak time in ticks
//                           double diffpeaktt0     = hitlist[ll]->PeakTime();
		           double allhitX = detprop->ConvertTicksToX(diffpeaktt0,hitlist[ll]->WireID().Plane,hitlist[ll]->WireID().TPC,hitlist[ll]->WireID().Cryostat); //convert ticks to usec (1tick = 0.5usec), and subtract the T0
	                   double Wirestart[3], Wireend[3];
                           geometry->WireEndPoints(hitlist[ll]->WireID().Cryostat, hitlist[ll]->WireID().TPC, hitlist[ll]->WireID().Plane, hitlist[ll]->WireID().Wire, Wirestart, Wireend);
			   double allhitZ = Wirestart[2];
			   double allhitY = 0;
                           if(fmsp.at(ll).size()) 
			   allhitY = fmsp.at(ll)[0]->XYZ()[1];
			   if(!fmsp.at(ll).size()) 
			   allhitY = 303.5;			   			    
 			   
	                   if(same_trk==0 && long_trk==0 && hitlist[ll]->WireID().TPC==endtpcno &&(allhitZ<APAnegbound1 || allhitZ>APAposbound1) && (allhitZ<APAnegbound2 || allhitZ>APAposbound2))
		           {
//		             if(abs(allhitX-endhitx)<_absxdiff )
//		             {	
//		               if( abs(allhitZ-endhitz)<=_abszdiff )

                               double hitdist = sqrt(pow(abs(allhitX-trkstopx),2)+pow(abs(allhitZ-trkstopz),2)); 
                               if(hitdist <= _hitdist)
		               {
	          		  nearhits_key[nearhitct]     = hitlist[ll].key();
		                  nearhits_peakT[nearhitct]   = hitlist[ll]->PeakTime();
		                  nearhits_charge[nearhitct]  = hitlist[ll]->Integral();
		                  nearhits_wire[nearhitct]    = hitlist[ll]->WireID().Wire;
		                  nearhits_chno[nearhitct]    = hitlist[ll]->Channel();
		                  nearhits_TPC[nearhitct]     = hitlist[ll]->WireID().TPC;
		                  nearhits_plane[nearhitct]   = hitlist[ll]->WireID().Plane;
				  nearhits_xpos[nearhitct]    = allhitX;
				  nearhits_ypos[nearhitct]    = allhitY;
				  nearhits_zpos[nearhitct]    = allhitZ;
	                          nearhits_mult[nearhitct]      = hitlist[ll]->Multiplicity();
	                          nearhits_sigptime[nearhitct]  = hitlist[ll]->SigmaPeakTime();
	                          nearhits_sigchrg[nearhitct]   = hitlist[ll]->SigmaIntegral();
	                          nearhits_sigpamp[nearhitct]	= hitlist[ll]->SigmaPeakAmplitude();
	                          nearhits_dof[nearhitct]	= hitlist[ll]->DegreesOfFreedom();
	                          nearhits_gof[nearhitct]       = hitlist[ll]->GoodnessOfFit();
		                  nearhits_ptminusRMS[nearhitct]= hitlist[ll]->PeakTimeMinusRMS(5.0);
		                  nearhits_ptplusRMS[nearhitct] = hitlist[ll]->PeakTimePlusRMS(5.0);

                                  std::array<float,4> cnn_out = hitResults.getOutput( hitlist[ll].key() );
//                                  double p_trk_or_sh = cnn_out[ hitResults.getIndex("track") ]+ cnn_out[ hitResults.getIndex("em") ]+ cnn_out[ hitResults.getIndex("michel") ]; 
//                                  double cnn_score = cnn_out[ hitResults.getIndex("michel") ]; 
		                  nearhits_cnnMichel[nearhitct] = cnn_out[ hitResults.getIndex("michel") ]; 
		                  nearhits_cnnEM[nearhitct]     = cnn_out[ hitResults.getIndex("em") ]; 
		                  nearhits_cnnTrack[nearhitct]  = cnn_out[ hitResults.getIndex("track") ]; 

 		                 
				  cout<<"hitlist[ll].key() "<<hitlist[ll].key()<<" peakT "<<hitlist[ll]->PeakTime()<<" wireID "<<hitlist[ll]->WireID().Wire<<" x, y, z "<<nearhits_xpos[nearhitct]<<" "<<nearhits_ypos[nearhitct]<<" "<<nearhits_zpos[nearhitct]<<"\n";		 
		                  cout<<"Last hit: track ID "<<track.ID()<<" wireno: "<<endwireno<<" peaktime:  "<<endpeaktime<<" tpcno: "<<endtpcno<<" hitchrg: "<<endhitchrg<<" This hit: wire: "<<hitlist[ll]->WireID().Wire<<" peaktime "<<hitlist[ll]->PeakTime()<<" TPC: "<<hitlist[ll]->WireID().TPC<<"\n";
	                         
				  nearhitct++;
	                       }
//		             }
		           }
	                 }//if(hitlist[ll]->WireID().Plane==2)
	               }//for(size_t ll=0; ll<hitlist.size();++ll)
	    		
			
	               cout<<" candidate muon coll. plane hits "<<trkcolhits<<"\n";
 	     
 		         //----------Looping over showers-------------//
		 
                         double shwr_dist        = 99999, shwr_dist0 = 99999; 		 
                         int    shwr_key         = -999;
                         int    shwr_ID          = -999;
                         double shwr_length      = -999;
                         double shwr_startx      = -999;
                         double shwr_starty      = -999;
                         double shwr_startz      = -999;
                         int    shwr_bestplane   = -999;
                         double shwr_startdcosx  = -999;
                         double shwr_startdcosy  = -999;
                         double shwr_startdcosz  = -999;
                         double shwr_openangle   = -999;
//                         double shwr_dEdx        = -999;
//                         double shwr_energy      = -999;
//                         double shwr_mipenergy   = -999;
 
                         cout<<"showerlist.size() "<<showerlist.size()<<"\n";
                         for(size_t jjj=0; jjj<showerlist.size();++jjj)
                         {
                           art::Ptr<recob::Shower> pshower(showerListHandle, jjj);
                           const recob::Shower& shower = *pshower;

                       std::vector<art::Ptr<recob::PFParticle>> pfpsh=pfp_shwr_assn.at(jjj);
                       if(pfpsh.size()) 
                       { 
        	
                        //---Only consider T0 tagged shower-----//

                         std::vector<art::Ptr<anab::T0>> t0sh=shwr_t0_assn_v.at(pfpsh[0].key());
       	                 if(t0sh.size() )
	                 {  
                           TVector3 const& shwr_start = shower.ShowerStart();
                           float shwr_startx1     = shwr_start.X();
                           float shwr_starty1     = shwr_start.Y();
                           float shwr_startz1     = shwr_start.Z();
  
                           shwr_dist = sqrt(pow(shwr_startx1 - trkstopx,2) + pow(shwr_startz1 - trkstopz,2));
//                           shwr_dist = sqrt(pow(shwr_startx1 - endhitx,2)+ pow(shwr_starty1 - endhity,2) + pow(shwr_startz1-endhitz,2));
//                           cout<<"shwr_dist "<<shwr_dist<<" shwr_dist0 "<<shwr_dist0<<"\n";

                           if( shwr_dist<shwr_dist0)
                           {			   
			     //cout<<"jjj "<<jjj<<" showerlist[jjj].key() "<<showerlist[jjj].key()<<" shwr_ID "<<shower.ID()<<"\n";
                             shwr_dist0       = shwr_dist;
                             shwr_key         = showerlist[jjj].key();
                             shwr_ID          = shower.ID();
                             shwr_length      = shower.Length();
                             shwr_startx      = shwr_startx1;
                             shwr_starty      = shwr_starty1;
                             shwr_startz      = shwr_startz1;
                             shwr_bestplane   = shower.best_plane();
                             TVector3 const& shwrdir_start = shower.Direction();
                             shwr_startdcosx  = shwrdir_start.X();
                             shwr_startdcosy  = shwrdir_start.Y();
                             shwr_startdcosz  = shwrdir_start.Z();  
                             shwr_openangle   = shower.OpenAngle();
//                             shwr_dEdx        = shower.dEdx().at(shower.best_plane());
//                             shwr_energy      = shower.Energy().at(shower.best_plane());
//                             shwr_mipenergy   = shower.MIPEnergy().at(shower.best_plane());  
/*
                             if (shower.Energy().size() == 3)
                             std::copy_n
                             (shower.Energy().begin(),    3, &shwr_energy);
                             if (shower.dEdx().size() == 3)
                             std::copy_n
                             (shower.dEdx().begin(),      3, &shwr_dEdx);
                             if (shower.MIPEnergy().size() == 3)
                             std::copy_n
                             (shower.MIPEnergy().begin(), 3, &shwr_mipenergy);
*/
                           }
			   }
			   }
                        }
			
 	               fnearhitcountana[fhitctana] = nearhitct;   
		       fnshwrdistana[fhitctana] = shwr_dist0;
	               fMichelcountnearhitana[fhitctana] = ftrueMichel;
		       ftrueEnearhitana[fhitctana] = mcd_energy;
                       fhitctana++;

	               ///*****************Nearby hit count cut******************//   
		       
		       cout<<"nearhitct "<<nearhitct<<"\n";
	               if(nearhitct>= _minhitcountmichel && nearhitct< _maxhitcountmichel)
	               {
		         fnearhits_trks++;
			 
	                 _fnearhits_trks++;
	                 if(ftrueMichel==1)_ftruenearhits_trks++;
	                 cout<<"Nearby hit count "<<nearhitct<<"\n";		

                        
                          std::vector<double> HitPeakTimesshwrall;
                          std::vector<double> HitPeakTimesshwrcol;
 
			  
                          std::vector<art::Ptr<recob::Hit >> shwrhits = fshwrHit.at(shwr_key);	

 //                          std::vector<const recob::Hit* > shwrhits = fshwrhit.at(shwr_key);	
   
                          //cout<<"shwrhits.size() "<<shwrhits.size()<<"\n";	 
     
//                         for(std::vector<const recob::Hit* >::iterator itr = shwrhits.begin(); itr < shwrhits.end(); itr++)
//                         {  
                          for(size_t itr = 0; itr<shwrhits.size(); ++itr)
		          {
//                           cout<<"inside shower hits plane: "<<shwrhits[itr]->WireID().Plane<<"\n";
                           HitPeakTimesshwrall.push_back(shwrhits[itr]->PeakTime());
			   
                           // looping over the collection plane  hits of the selected shower
                           if(shwrhits[itr]->WireID().Plane == 2)
                           {
                             HitPeakTimesshwrcol.push_back(shwrhits[itr]->PeakTime());
                           }
                         }
			 float MinHitPeakTimeshwrall = 0;
			 float MinHitPeakTimeshwrcol = 0;
			 float MaxHitPeakTimeshwrall = 0;
			 float MaxHitPeakTimeshwrcol = 0;
			 if(HitPeakTimesshwrall.size() && HitPeakTimesshwrcol.size())
			 {
                           MinHitPeakTimeshwrall = *(std::min_element(HitPeakTimesshwrall.begin(), HitPeakTimesshwrall.end()));
                           MinHitPeakTimeshwrcol = *(std::min_element(HitPeakTimesshwrcol.begin(), HitPeakTimesshwrcol.end()));
                           MaxHitPeakTimeshwrall = *(std::max_element(HitPeakTimesshwrall.begin(), HitPeakTimesshwrall.end()));
                           MaxHitPeakTimeshwrcol = *(std::max_element(HitPeakTimesshwrcol.begin(), HitPeakTimesshwrcol.end()));
			 }
			   cout<<"MinHitPeakTimeshwrall "<<MinHitPeakTimeshwrall<<" MinHitPeakTimeshwrcol "<<MinHitPeakTimeshwrcol<<"\n";
			   cout<<"MaxHitPeakTimeshwrall "<<MaxHitPeakTimeshwrall<<" MaxHitPeakTimeshwrcol "<<MaxHitPeakTimeshwrcol<<"\n";
			 
                        HitPeakTimesshwrall.clear();
			HitPeakTimesshwrcol.clear();
			
                        fshwrdistana[fshwrdisana] = shwr_dist0;
	                fMichelcountshwrdistana[fshwrdisana] = ftrueMichel;
			ftrueEshwrdistana[fshwrdisana] = mcd_energy;
	                fshwrdisana++;

			//---------The final selection cut---------------------//
 			 if(abs(shwr_dist0)< _shwr_dist && MinHitPeakTimeshwrcol > _minhitpeakcut  && MaxHitPeakTimeshwrcol < _maxhitpeakcut)
                         {
	                  _fshwr_trks++;
	                  if(ftrueMichel==1)_ftrueshwr_trks++;
			  
	                  cout<<"\n"<<"\n"<<"Selected the michel event!!!! run, subrun, event numbers: "<<frun<<" "<<fsubrun<<" "<<fevent<<" TrkackID "<<track.ID()<<" start x, y, z "<<trkbegx<<" "<<trkbegy<<"  "<<trkbegz<<" end x, y, z "<<trkstopx<<" "<<trkstopy<<"  "<<trkstopz<<"  Track length "<<track.Length()<<"\n";
 
                        cout<<"shwr_dist0 "<<shwr_dist0<<" fshwr_key "<<shwr_key<<" fshwr_startx "<<shwr_startx<<" fshwr_starty "<<shwr_starty<<" fshwr_startz "<<shwr_startz<<"\n";


 			 //----------Order near hits such that the distance from the muon end point increases--------------------//
			 cout<<" end xpos "<<endhitx<<" end zpos "<<endhitz<<"\n";
			 float hits_dis1 = 0, hits_dis2 = 0;
			 float aaa=0, bbb=0, ccc=0, ddd=0, eee=0, fff=0, ggg=0, hhh=0, iii=0, mmm=0, jjj=0, kkk=0, lll=0, nnn=0, ooo=0, ppp=0, qqq=0, rrr=0, sss=0, ttt=0, uuu=0;

                         for(int qq =0; qq < nearhitct; qq++)
		         {
                           for(int rr =qq+1; rr < nearhitct; rr++)
		           {
			     hits_dis1 = sqrt(pow((nearhits_xpos[qq]-endhitx),2)+pow((nearhits_zpos[qq]-endhitz),2));
			     hits_dis2 = sqrt(pow((nearhits_xpos[rr]-endhitx),2)+pow((nearhits_zpos[rr]-endhitz),2));
			     if(hits_dis1 > hits_dis2)
			     {
			       aaa = nearhits_key[qq];
	          	       nearhits_key[qq] = nearhits_key[rr];
			       nearhits_key[rr] = aaa;
			       
			       bbb = nearhits_charge[qq];
			       nearhits_charge[qq] = nearhits_charge[rr];
			       nearhits_charge[rr] = bbb;
			       
			       ccc = nearhits_wire[qq];
			       nearhits_wire[qq] = nearhits_wire[rr];
			       nearhits_wire[rr] = ccc;
			       
			       ddd = nearhits_chno[qq];
			       nearhits_chno[qq] = nearhits_chno[rr];
			       nearhits_chno[rr] = ddd;
			       
			       eee = nearhits_TPC[qq];
			       nearhits_TPC[qq] = nearhits_TPC[rr];
			       nearhits_TPC[rr] = eee;
			       
			       fff = nearhits_plane[qq];
			       nearhits_plane[qq] = nearhits_plane[rr];
			       nearhits_plane[rr] = fff;
			       
			       ggg = nearhits_xpos[qq];
			       nearhits_xpos[qq] = nearhits_xpos[rr];
			       nearhits_xpos[rr] = ggg;
			       
			       mmm = nearhits_ypos[qq];
			       nearhits_ypos[qq] = nearhits_ypos[rr];
			       nearhits_ypos[rr] = mmm;

			       hhh = nearhits_zpos[qq];
			       nearhits_zpos[qq] = nearhits_zpos[rr];
			       nearhits_zpos[rr] = hhh;
			       
			       iii = nearhits_peakT[qq];
			       nearhits_peakT[qq] = nearhits_peakT[rr];
			       nearhits_peakT[rr] = iii;

			       jjj = nearhits_mult[qq];
			       nearhits_mult[qq] = nearhits_mult[rr];
			       nearhits_mult[rr] = jjj;

			       kkk = nearhits_sigptime[qq];
			       nearhits_sigptime[qq] = nearhits_sigptime[rr];
			       nearhits_sigptime[rr] = kkk;

			       lll = nearhits_sigchrg[qq];
			       nearhits_sigchrg[qq] = nearhits_sigchrg[rr];
			       nearhits_sigchrg[rr] = lll;

			       nnn = nearhits_sigpamp[qq];
			       nearhits_sigpamp[qq] = nearhits_sigpamp[rr];
			       nearhits_sigpamp[rr] = nnn;

			       ooo = nearhits_dof[qq];
			       nearhits_dof[qq] = nearhits_dof[rr];
			       nearhits_dof[rr] = ooo;

			       ppp = nearhits_gof[qq];
			       nearhits_gof[qq] = nearhits_gof[rr];
			       nearhits_gof[rr] = ppp;

			       qqq = nearhits_ptminusRMS[qq];
			       nearhits_ptminusRMS[qq] = nearhits_ptminusRMS[rr];
			       nearhits_ptminusRMS[rr] = qqq;

			       rrr = nearhits_ptplusRMS[qq];
			       nearhits_ptplusRMS[qq] = nearhits_ptplusRMS[rr];
			       nearhits_ptplusRMS[rr] = rrr;

			       sss = nearhits_cnnMichel[qq];
			       nearhits_cnnMichel[qq] = nearhits_cnnMichel[rr];
			       nearhits_cnnMichel[rr] = sss;

			       ttt = nearhits_cnnEM[qq];
			       nearhits_cnnEM[qq] = nearhits_cnnEM[rr];
			       nearhits_cnnEM[rr] = ttt;

			       uuu = nearhits_cnnTrack[qq];
			       nearhits_cnnTrack[qq] = nearhits_cnnTrack[rr];
			       nearhits_cnnTrack[rr] = uuu;
			       }			       
			     }			 
			   }  
                         for(int qq =0; qq < nearhitct; qq++)
		         {
			 cout<<"nearby hits_dis "<<sqrt(pow((nearhits_xpos[qq]-endhitx),2)+pow((nearhits_zpos[qq]-endhitz),2))<<" TPC "<<nearhits_TPC[qq]<<" wire "<<nearhits_wire[qq]<<" peakt "<<nearhits_peakT[qq]<<" x, y, z "<<nearhits_xpos[qq]<<" "<<nearhits_ypos[qq]<<" "<<nearhits_zpos[qq]<<"\n";
			 }
                              //-------------checking the correlation of the michel hits----------------//	       	 
                
		         float sumwire=0, sumptime=0, meanwire=0, meanptime=0, diffWirePtime=0, diffWire=0, diffPtime=0, stDWire=0, stDPtime=0, CorrWirePtime=0;
                         for(int qq =0; qq < nearhitct; qq++)
		         {
		           sumwire += nearhits_wire[qq];
		           sumptime += nearhits_peakT[qq];
		         }
		
		         meanwire = sumwire/nearhitct;
		         meanptime = sumptime/nearhitct;
				
		         for(int qq =0; qq < nearhitct; qq++)
		         {
		           diffWirePtime += (nearhits_wire[qq]-meanwire)*(nearhits_peakT[qq]-meanptime);
		           diffWire      += (nearhits_wire[qq]-meanwire)*(nearhits_wire[qq]-meanwire);
		           diffPtime     += (nearhits_peakT[qq]-meanptime)*(nearhits_peakT[qq]-meanptime);
		         }  
                
		         stDWire         = sqrt(diffWire/nearhitct);
		         stDPtime        = sqrt(diffPtime/nearhitct);
		         CorrWirePtime   = diffWirePtime/((nearhitct)*stDWire*stDPtime);
                         //cout<<"CorrWirePtime "<<CorrWirePtime<<"\n";
			 			   

                         // -------- Storing selected candidate nuon quantities -------------//
			 
                         fsel_run          = evt.run();
                         fsel_subrun       = evt.subRun();
                         fsel_event        = evt.id().event();
                         fsel_evttime      = tts.AsDouble();
		         fsel_endhitkey    = endhitkey;	
		         fsel_endwire      = endwireno;
		         fsel_endchno      = endchno;
		         fsel_endtpcno     = endtpcno;
		         fsel_endhitchrg   = endhitchrg;
		         fsel_endptime     = endpeaktime;
		         fsel_endhitx      = endhitx;
		         fsel_endhity      = endhity;
		         fsel_endhitz      = endhitz;
		         fsel_endhitmult      = endhitmult;
		         fsel_endhitsigptime  = endhitsigptime;
		         fsel_endhitsigchrg   = endhitsigchrg;
		         fsel_endhitsigpamp    = endhitsigpamp;
		         fsel_endhitdof       = endhitdof;
		         fsel_endhitgof       = endhitgof;
		         fsel_endhitptplusRMS = endhitptplusRMS;
		         fsel_endhitptminusRMS = endhitptminusRMS;
			 fsel_ccrosser     = ccrosser;
			 fsel_minhitptime  = MinHitPeakTime;
			 fsel_maxhitptime  = MaxHitPeakTime;
			 fsel_ncolhits     = trkcolhits;
		         fsel_nearhitcount = nearhitct;
		         fsel_CorrWirePtime= CorrWirePtime;
                         fsel_trkrecotime  = double(T00/1000); //track reco time T0 (convert from nsec to usec)
		 

	                 fsel_dist_hit_end = enddist;
	                 fsel_dist_times   = dt_min;
                         fsel_trackthetaxz = std::atan2(dir_start.X(), dir_start.Z());
                         fsel_trackthetayz = std::atan2(dir_start.Y(), dir_start.Z());
                         fsel_trkstartx    = trkbegx;  
                         fsel_trkstarty    = trkbegy;  
                         fsel_trkstartz    = trkbegz;  
                         fsel_trkendx      = trkstopx;
                         fsel_trkendy      = trkstopy;
                         fsel_trkendz      = trkstopz;
                         fsel_trkstartcosx = dir_start.X();
                         fsel_trkstartcosy = dir_start.Y();
                         fsel_trkstartcosz = dir_start.Z();
                         fsel_trkendcosx   = dir_end.X();
                         fsel_trkendcosy   = dir_end.Y();
                         fsel_trkendcosz   = dir_end.Z();
                         fsel_trklen       = track.Length();
		         fsel_trktheta     = dir_start.Theta();
		         fsel_trkphi       = dir_start.Phi();
                         fsel_trkID        = track.ID();
//                         fsel_PHratio      = avgdownchrg/avgupchrg;
//                         fsel_PH           = PH;
		
		         fMichelcountselana= ftrueMichel;
			 
			 //---------storing selected shower reco quantities -------------//
                         fselshwr_key        =  shwr_key;
                         fselshwr_ID         =  shwr_ID;
                         fselshwr_length     =  shwr_length;
                         fselshwr_startx     =  shwr_startx;
                         fselshwr_starty     =  shwr_starty;
                         fselshwr_startz     =  shwr_startz;
                         fselshwr_bestplane  =  shwr_bestplane;
                         fselshwr_startdcosx =  shwr_startdcosx;
                         fselshwr_startdcosx =  shwr_startdcosy;
                         fselshwr_startdcosx =  shwr_startdcosz;  
                         fselshwr_openangle  =  shwr_openangle;
			 fselshwr_dist       =  shwr_dist0;
//                         fselshwr_dEdx       =  shwr_dEdx;
//                         fselshwr_energy     =  shwr_energy;
//                         fselshwr_mipenergy  =  shwr_mipenergy;  


                         //-------------Saving selected candidate muon hits info------------------------//
			 
			 ftrkcolhits = 0;
			 for(int mm=0; mm<trkcolhits; mm++)
			 {
	                   fhits_key[ftrkcolhits]      =  hits_key[mm];
	                   fhits_charge[ftrkcolhits]   =  hits_charge[mm];
                           fhits_wire[ftrkcolhits]     =  hits_wire[mm];
	                   fhits_peakT[ftrkcolhits]    =  hits_peakT[mm];
	                   fhits_TPC[ftrkcolhits]      =  hits_TPC[mm];
		           fhits_chno[ftrkcolhits]     =  hits_chno[mm];
	                   fhits_xpos[ftrkcolhits]     =  hits_xpos[mm];
	                   fhits_ypos[ftrkcolhits]     =  hits_ypos[mm];
	                   fhits_zpos[ftrkcolhits]     =  hits_zpos[mm];
	                   fhits_mult[ftrkcolhits]     =  hits_mult[mm];
	                   fhits_sigptime[ftrkcolhits] =  hits_sigptime[mm];
	                   fhits_sigchrg[ftrkcolhits]  =  hits_sigchrg[mm];
	                   fhits_sigpamp[ftrkcolhits]   =  hits_sigpamp[mm];
	                   fhits_dof[ftrkcolhits]      =  hits_dof[mm];
	                   fhits_gof[ftrkcolhits]      =  hits_gof[mm];
		           fhits_ptminusRMS[ftrkcolhits]= hits_ptminusRMS[mm];
		           fhits_ptplusRMS[ftrkcolhits]	 = hits_ptplusRMS[mm];
		           fhits_cnnMichel[ftrkcolhits]	 = hits_cnnMichel[mm];
		           fhits_cnnEM[ftrkcolhits]	 = hits_cnnEM[mm];
		           fhits_cnnTrack[ftrkcolhits]	 = hits_cnnTrack[mm];
			   
			   //fhits_chrg.at(fsel_mu).push_back(hits_charge[mm]);
			   //cout<<"fsel_mu "<<fsel_mu<<" ftrkcolhits "<<ftrkcolhits<<" fhits_charge "<<fhits_charge[ftrkcolhits]<<"\n";
			   
			   ftrkcolhits++;
			 }
			 
	  
                         //-----------Filling selected shower hits info---------------------//

                          fshwrcolhits = 0;
			  fshwrallhits = 0;
 
 //                          std::vector<const recob::Hit* > shwrhits = fshwrhit.at(shwr_key);	
 
//                       for(std::vector<const recob::Hit* >::iterator itr = shwrhits.begin(); itr < shwrhits.end(); itr++)
                         for(size_t itr = 0; itr<shwrhits.size(); ++itr)
                         {  
			 			   
                           //fshwrallhits_key[fshwrallhits]     =   shwrhits[itr]->ID();
	                   fshwrallhits_chno[fshwrallhits]      =   shwrhits[itr]->Channel();
	                   fshwrallhits_peakT[fshwrallhits]     =   shwrhits[itr]->PeakTime();
	                   fshwrallhits_charge[fshwrallhits]    =   shwrhits[itr]->Integral();///(TMath::Exp(-(hits_peakT[longcolhits]-800)*500/tau)); //multiplied by 500nsec to convert time ticks to actual generation time
                           fshwrallhits_wire[fshwrallhits]      =   shwrhits[itr]->WireID().Wire;
			   fshwrallhits_plane[fshwrallhits]     =   shwrhits[itr]->WireID().Plane;
	                   fshwrallhits_TPC[fshwrallhits]       =   shwrhits[itr]->WireID().TPC;
                           double diffpeaktt0allshwr            =   shwrhits[itr]->PeakTime() - (T00/1000)*2; //calculate the peak time in ticks
	                   fshwrallhits_xpos[fshwrallhits]      =   detprop->ConvertTicksToX(diffpeaktt0allshwr,shwrhits[itr]->WireID().Plane,shwrhits[itr]->WireID().TPC,shwrhits[itr]->WireID().Cryostat); //convert ticks to usec (1tick = 0.5usec), and subtract the T0
                           double Wirestartallshwr[3], Wireendallshwr[3];
                           geometry->WireEndPoints(shwrhits[itr]->WireID().Cryostat, shwrhits[itr]->WireID().TPC, shwrhits[itr]->WireID().Plane, shwrhits[itr]->WireID().Wire, Wirestartallshwr, Wireendallshwr);
//			   fshwrallhits_ypos[fshwrallhits]      =   Wirestartallshwr[1];
			   fshwrallhits_zpos[fshwrallhits]      =   Wirestartallshwr[2];
                           if(fmsp.at(shwrhits[itr].key()).size())
			   fshwrallhits_ypos[fshwrallhits]      =   fmsp.at(shwrhits[itr].key())[0]->XYZ()[1];
			   if(!fmsp.at(shwrhits[itr].key()).size())
			   fshwrallhits_ypos[fshwrallhits]      =   303.5;
			   fshwrallhits_mult[fshwrallhits]      =   shwrhits[itr]->Multiplicity();
			   fshwrallhits_sigptime[fshwrallhits]  =   shwrhits[itr]->SigmaPeakTime();
			   fshwrallhits_sigchrg[fshwrallhits]   =   shwrhits[itr]->SigmaIntegral();
			   fshwrallhits_sigpamp[fshwrallhits]	=   shwrhits[itr]->SigmaPeakAmplitude();
			   fshwrallhits_dof[fshwrallhits]       =   shwrhits[itr]->DegreesOfFreedom();
			   fshwrallhits_gof[fshwrallhits]       =   shwrhits[itr]->GoodnessOfFit();
		           fshwrallhits_ptminusRMS[fshwrallhits]= shwrhits[itr]->PeakTimeMinusRMS(5.0);
		           fshwrallhits_ptplusRMS[fshwrallhits]	 = shwrhits[itr]->PeakTimePlusRMS(5.0);

//			   cout<<" shwrhits[itr].key() "<<shwrhits[itr].key()<<" fshwrallhits_ypos "<<fshwrallhits_ypos[fshwrallhits]<<"\n";
			   
			   fshwrallhits++;
			   
                           // looping over the collection plane  hits of the selected shower
                           if(shwrhits[itr]->WireID().Plane == 2)
                           {


	                     fshwrhits_chno[fshwrcolhits]      =   shwrhits[itr]->Channel();
	                     fshwrhits_peakT[fshwrcolhits]     =   shwrhits[itr]->PeakTime();
//	                     chrghit->fshwrhits_charge[2][fshwrcolhits]    =   ((*itr)->Integral());///(TMath::Exp(-(hits_peakT[longcolhits]-800)*500/tau)); //multiplied by 500nsec to convert time ticks to actual generation time

	                     fshwrhits_charge[fshwrcolhits]    =   shwrhits[itr]->Integral();///(TMath::Exp(-(hits_peakT[longcolhits]-800)*500/tau)); //multiplied by 500nsec to convert time ticks to actual generation time
                             fshwrhits_wire[fshwrcolhits]      =   shwrhits[itr]->WireID().Wire;
			     fshwrhits_plane[fshwrcolhits]     =   shwrhits[itr]->WireID().Plane;
	                     fshwrhits_TPC[fshwrcolhits]       =   shwrhits[itr]->WireID().TPC;
                             double diffpeaktt0shwr           =    shwrhits[itr]->PeakTime() - (T00/1000)*2; //calculate the peak time in ticks
	                     fshwrhits_xpos[fshwrcolhits]      =   detprop->ConvertTicksToX(diffpeaktt0shwr,shwrhits[itr]->WireID().Plane,shwrhits[itr]->WireID().TPC,shwrhits[itr]->WireID().Cryostat); //convert ticks to usec (1tick = 0.5usec), and subtract the T0
                             double Wirestartshwr[3], Wireendshwr[3];
                             geometry->WireEndPoints(shwrhits[itr]->WireID().Cryostat, shwrhits[itr]->WireID().TPC, shwrhits[itr]->WireID().Plane, shwrhits[itr]->WireID().Wire, Wirestartshwr, Wireendshwr);
			     fshwrhits_zpos[fshwrcolhits]      =   Wirestartshwr[2];
                             if(fmsp.at(shwrhits[itr].key()).size())
			     fshwrhits_ypos[fshwrcolhits]      =   fmsp.at(shwrhits[itr].key())[0]->XYZ()[1];
			     if(!fmsp.at(shwrhits[itr].key()).size())
			     fshwrhits_ypos[fshwrcolhits]      =   303.5;
			     fshwrhits_mult[fshwrcolhits]      =   shwrhits[itr]->Multiplicity();
			     fshwrhits_sigptime[fshwrcolhits]  =   shwrhits[itr]->SigmaPeakTime();
			     fshwrhits_sigchrg[fshwrcolhits]   =   shwrhits[itr]->SigmaIntegral();
			     fshwrhits_sigpamp[fshwrcolhits]    =   shwrhits[itr]->SigmaPeakAmplitude();
			     fshwrhits_dof[fshwrcolhits]       =   shwrhits[itr]->DegreesOfFreedom();
			     fshwrhits_gof[fshwrcolhits]       =   shwrhits[itr]->GoodnessOfFit();
		             fshwrhits_ptminusRMS[fshwrcolhits]	  = shwrhits[itr]->PeakTimeMinusRMS(5.0);
		             fshwrhits_ptplusRMS[fshwrcolhits]	 = shwrhits[itr]->PeakTimePlusRMS(5.0);
	 

			     //cout<<"fshwrcolhits "<<fshwrcolhits[fsel_mu]<<" fshwrhits_charge "<<fshwrhits_charge[fsel_mu][fshwrcolhits[fsel_mu]]<<"\n";
	                     //cout<<"shwrhits_peakT "<<shwrhits_peakT[shwrcolhits]<<" shwrhits_charge "<<shwrhits_charge[shwrcolhits]<<" shwrhits_wire "<<shwrhits_wire[shwrcolhits]<<" shwrhits_channel "<<shwrhits_channel[shwrcolhits]<<" shwrhits_TPC "<<shwrhits_TPC[shwrcolhits]<<" shwrhits_xpos "<<shwrhits_xpos[shwrcolhits]<<"\n";
			     fshwrcolhits++;
                           }
                         }
			 	
                         //-------------Filling number of nearby hits--------------------//
                         fnearhitct = nearhitct;
			 
                         for(int lp=0; lp<fnearhitct; lp++)
		         {
	          	   fnearhits_key[lp]     = nearhits_key[lp];
		           fnearhits_peakT[lp]   = nearhits_peakT[lp];
		           fnearhits_charge[lp]  = nearhits_charge[lp];
		           fnearhits_wire[lp]    = nearhits_wire[lp];
		           fnearhits_chno[lp]    = nearhits_chno[lp];
		           fnearhits_TPC[lp]     = nearhits_TPC[lp];
		           fnearhits_plane[lp]   = nearhits_plane[lp];
			   fnearhits_xpos[lp]    = nearhits_xpos[lp];
		           fnearhits_ypos[lp]    = nearhits_ypos[lp];			   
		           fnearhits_zpos[lp]    = nearhits_zpos[lp];
 		           fnearhits_mult[lp]    = nearhits_mult[lp];
		           fnearhits_sigptime[lp]= nearhits_sigptime[lp];
		           fnearhits_sigchrg[lp] = nearhits_sigchrg[lp];
		           fnearhits_sigpamp[lp]  = nearhits_sigpamp[lp];
		           fnearhits_dof[lp]     = nearhits_dof[lp];
		           fnearhits_gof[lp]     = nearhits_gof[lp];
		           fnearhits_ptminusRMS[lp]= nearhits_ptminusRMS[lp];
		           fnearhits_ptplusRMS[lp] = nearhits_ptplusRMS[lp];
		           fnearhits_cnnMichel[lp] = nearhits_cnnMichel[lp];
		           fnearhits_cnnEM[lp]     = nearhits_cnnEM[lp];
		           fnearhits_cnnTrack[lp]  = nearhits_cnnTrack[lp];

			 }  
				  
			 		 

		          //------------Doing some calorimetry-----------------------------------------//
		  		  
                          std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(i);
			  
			  cout<<"calos.size() "<<calos.size()<<"\n";
			  
			    fhitsU = 0;
			    fhitsV = 0;
			    fhitsY = 0;

                          for(size_t ical = 0; ical<calos.size(); ++ical)
		          {
	                    if(!calos[ical]) continue;
	                    if(!calos[ical]->PlaneID().isValid) continue;
	                    int planenum = calos[ical]->PlaneID().Plane;
	                    if(planenum<0 || planenum>2) continue;
	                    const size_t NHits = calos[ical] -> dEdx().size();
	                    fntrkhits=int(NHits);
			    //cout<<"fntrkhits "<<fntrkhits<<" planenum "<<planenum<<"\n";
	                    for(size_t iHit = 0; iHit < NHits; ++iHit)
		            {
	                      const auto& TrkPos = (calos[ical] -> XYZ())[iHit];
			      
			      if(planenum == 0)
			      {

	                      ftrkdqdxU[fhitsU]=(calos[ical] -> dQdx())[iHit];
	                      ftrkdedxU[fhitsU]=(calos[ical] -> dEdx())[iHit];
	                      ftrkresrangeU[fhitsU]=(calos[ical]->ResidualRange())[iHit];
	                      ftrkhitxU[fhitsU]=TrkPos.X();
	                      ftrkhityU[fhitsU]=TrkPos.Y();
	                      ftrkhitzU[fhitsU]=TrkPos.Z();
	                      ftrkpitchU[fhitsU]=(calos[ical]->TrkPitchVec())[iHit];

			      fhitsU++;
			      }
			      if(planenum == 1)
			      {
	                      ftrkdqdxV[fhitsV]=(calos[ical] -> dQdx())[iHit];
	                      ftrkdedxV[fhitsV]=(calos[ical] -> dEdx())[iHit];
	                      ftrkresrangeV[fhitsV]=(calos[ical]->ResidualRange())[iHit];
	                      ftrkhitxV[fhitsV]=TrkPos.X();
	                      ftrkhityV[fhitsV]=TrkPos.Y();
	                      ftrkhitzV[fhitsV]=TrkPos.Z();
	                      ftrkpitchV[fhitsV]=(calos[ical]->TrkPitchVec())[iHit];
			      fhitsV++;
			      }
			      if(planenum == 2)
			      {
			      
	                      ftrkdqdxY[fhitsY]=(calos[ical] -> dQdx())[iHit];
	                      ftrkdedxY[fhitsY]=(calos[ical] -> dEdx())[iHit];
			      //cout<<"ftrkdqdxY["<<fhitsY<<"] "<<ftrkdqdxY[fhitsY]<<"\n";
	                      ftrkresrangeY[fhitsY]=(calos[ical]->ResidualRange())[iHit];
	                      ftrkhitxY[fhitsY]=TrkPos.X();
	                      ftrkhityY[fhitsY]=TrkPos.Y();
	                      ftrkhitzY[fhitsY]=TrkPos.Z();
	                      ftrkpitchY[fhitsY]=(calos[ical]->TrkPitchVec())[iHit];
			      fhitsY++;
			      }
	                    } // loop over iHit..
                          } // loop over ical 2nd time...
			  
			  //cout<<"fhitsU "<<fhitsU<<" fhitsV"<<fhitsV <<" fhitsY "<<fhitsY<<"\n";
  		   

/*                        //------------Making a cone to see which other hits might belong to a michel--------------//
			 
			 // first fit the already selected nearby michel hits to a straight line
                         int counter =0;
	                 int wire=0, wiresq=0;
                         double wirePtime=0, time=0;
                         double slope=-999, yintercept=-999, theta=-999;
	       
                         for(int i3=0; i3 <fnearhitct; i3++)
                         {
                           counter++;
                           wirePtime += fnearhits_zpos[i3]*fnearhits_xpos[i3];
                           time      += fnearhits_xpos[i3];
                           wire      += fnearhits_zpos[i3];
                           wiresq    += fnearhits_zpos[i3]*fnearhits_zpos[i3];  
                         }   

                         slope        = (fnearhits_xpos[fnearhitct-1] - fnearhits_xpos[0])/(fnearhits_zpos[fnearhitct-1] - fnearhits_zpos[0]);
                         yintercept   = (wiresq*time - wirePtime*wire)/(wiresq*counter - wire*wire);
//			 theta = std::atan(slope)*180/3.14;

//		         if(endhitx > 0)
//			 {
			   if(fnearhits_zpos[fnearhitct-1] > endhitz)
			   {
			     theta = std::atan(slope)*180/3.14;
			   }
			   if(fnearhits_xpos[fnearhitct-1] >= endhitx && fnearhits_zpos[fnearhitct-1] < endhitz)
			   {
			     theta = (std::atan(slope)*180/3.14)+180;
			   }
//			   if(fnearhits_xpos[fnearhitct-1] < endhitx && fnearhits_zpos[fnearhitct-1] > endhitz)
//			   {
// 			     theta = 360+std::atan(slope)*180/3.14;
//			   }
//			   if(fnearhits_xpos[fnearhitct-1] < endhitx && fnearhits_zpos[fnearhitct-1] < endhitz)
			   {
			     theta = (std::atan(slope)*180/3.14)-180;
			   }
//			 }  

			 cout<<" nume "<<(fnearhits_xpos[fnearhitct-1] - fnearhits_xpos[0])<<" denom "<<(fnearhits_zpos[fnearhitct-1] - endhitz)<<" slope "<<slope<<" theta "<<theta<<" yintercept "<<yintercept<<"\n";
                         //slope        = (wirePtime*counter - time*wire)/(wiresq*counter - wire*wire);
		 
		         //lets calculate the x-coordinate at distance d = 80 cm
		
		         //double xx = fnearhits_xpos[0] + sqrt((80*80)/(1+(slope*slope)));
		
		         // now find the corresponding z coordinate
			 //double zz = slope*(xx-fnearhits_xpos[0])+ fnearhits_zpos[0];
		         //double zz = slope*xx + yintercept;
		         //double magxz = sqrt(xx*xx+zz*zz);
			 //double dirline = std::atan(slope);
			 	
*/
   //---------------------------Taking the midpoint of the linear nearby hits----------//
    int counter =0;
    int zpos=0, zpossq=0;
    double xposzpos=0, xpos=0;
    double slope=-999, yintercept=-999;
    for(int i=0; i <fnearhitct; i++)
    {
      counter++;
      xposzpos += fnearhits_xpos[i]*fnearhits_zpos[i];
      xpos += fnearhits_xpos[i];
      zpos += fnearhits_zpos[i];
      zpossq += fnearhits_zpos[i]*fnearhits_zpos[i];  
    }   
    slope = (xposzpos*counter - xpos*zpos)/(zpossq*counter - zpos*zpos);
    yintercept = (zpossq*xpos - xposzpos*zpos)/(zpossq*counter - zpos*zpos);
    
    float xposend = slope*fnearhits_zpos[int(fnearhitct/2)]+yintercept;

                         /*********************For candidate muon truth information*******************************/

                         std::vector <double> trueHitsKey;
			 
	                 if(isMC)
                         {
                           // Get true MCParticle associated with recob::Track
                           const simb::MCParticle *particleP = truthUtil.GetMCParticleFromRecoTrack(track,evt,fTrackModuleLabel);
                           if(!particleP) continue;
                           const art::Ptr<simb::MCTruth> mcP=pi_serv->TrackIdToMCTruth_P(particleP->TrackId());
                           if(!mcP) continue;

/*                           std::vector<art::Ptr<recob::Hit>> allmuHits=fmtht.at(i); //storing hits for ith track
                           std::map<int,double> trkide;
	                   int trackid=-1;

                           for(size_t h=0; h<allmuHits.size();h++)
                           {
	                     art::Ptr<recob::Hit> hit=allmuHits[h];
	                     std::vector<sim::TrackIDE> eveIDs = bt_serv->HitToTrackIDEs(hit);
	                     for(size_t e=0;e<eveIDs.size(); ++e)
	                     {
	                       trkide[eveIDs[e].trackID] += eveIDs[e].energy;
	                     }
                           }
                           double  maxe = -1;
                           double tote = 0;
                           for(std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii)
                           {
	                     //	cout<<" trkid "<<ii->first<<"  energy deposited = "<<ii->second<<endl;
	                     tote += ii->second;
	                     if((ii->second)>maxe)
	                     {
	                       maxe = ii->second;
	                       trackid = ii->first;
                             }
                           }
*//*                         float total_energy=0.0;
                           for(size_t h=0; h<allmuHits.size();h++)
	                   {
	                     art::Ptr<recob::Hit> hit=allmuHits[h];
	                     if (hit->WireID().Plane!=2) continue;
	                     std::vector<sim::TrackIDE> eveIDs = bt_serv->HitToTrackIDEs(hit);
	                     for(size_t e=0;e<eveIDs.size(); ++e)
		             {
	                      if (eveIDs[e].trackID == trackid) total_energy +=  eveIDs[e].energy;
	                    }
                          }
 
                          const simb::MCParticle *particleP = pi_serv->TrackIdToParticle_P(trackid);
                          if(!particleP) continue;
                          const art::Ptr<simb::MCTruth> mcP=pi_serv->TrackIdToMCTruth_P(trackid);
                          if(!mcP) continue;
	                  //const simb::Origin_t parto(mctruthproto->GetParticle(iPartp));
*/
	                  TLorentzVector mcstart, mcend, mcstartdrifted, mcenddrifted;
	                  unsigned int pstarti, pendi, pstartdriftedi, penddriftedi; //mcparticle indices for starts and ends in tpc or drifted volumes
	                  double plen = length(*particleP, mcstart, mcend, pstarti, pendi);
	                  double plendrifted = driftedLength(*particleP, mcstartdrifted, mcenddrifted, pstartdriftedi, penddriftedi);

	                  // bool isActive = plen != 0;
	                  bool isDrifted = plendrifted!= 0;

                          unsigned int firstPointInAV = TrueParticleFirstPointInAV(fiducialBounds,*particleP);
                          unsigned int lastPointInAV = TrueParticleLastPointInAV(fiducialBounds,*particleP);
                          cout<<"track id "<<i<<" bt sim trackid "<<particleP->TrackId()<<"\n";
 
                          fmcsel_trkid      = particleP->TrackId();
                          fmcsel_vx         = particleP->Vx(firstPointInAV);
                          fmcsel_vy         = particleP->Vy(firstPointInAV);
                          fmcsel_vz         = particleP->Vz(firstPointInAV);
                          fmcsel_t          = particleP->T(firstPointInAV);
                          fmcsel_endx       = particleP->Vx(lastPointInAV);
                          fmcsel_endy       = particleP->Vy(lastPointInAV);
                          fmcsel_endz       = particleP->Vz(lastPointInAV);
                          fmcsel_endt       = particleP->T(lastPointInAV);
                          fmcsel_px         = particleP->Px();
                          fmcsel_py         = particleP->Py();
                          fmcsel_pz         = particleP->Pz();
                          fmcsel_momentum   = particleP->P();
                          fmcsel_energy     = particleP->E();
                          fmcsel_endpx      = particleP->EndPx();
                          fmcsel_endpy      = particleP->EndPy();
                          fmcsel_endpz      = particleP->EndPz();
                          fmcsel_endenergy  = particleP->EndE();
	                  fmcsel_pathlen    = plen;
			  fmcsel_length     = truthUtil.GetMCParticleLengthInTPCActiveVolume(*particleP,fiducialBounds[0],fiducialBounds[1],fiducialBounds[2],fiducialBounds[3],fiducialBounds[4],fiducialBounds[5]);

	                  if (isDrifted)
	                  {
                            fmcsel_vxdrifted         = particleP->Vx(firstPointInAV);
                            fmcsel_vydrifted         = particleP->Vy(firstPointInAV);
                            fmcsel_vzdrifted         = particleP->Vz(firstPointInAV);
                            fmcsel_tdrifted          = particleP->T(firstPointInAV);
                            fmcsel_endxdrifted       = particleP->Vx(lastPointInAV);
                            fmcsel_endydrifted       = particleP->Vy(lastPointInAV);
                            fmcsel_endzdrifted       = particleP->Vz(lastPointInAV);
                            fmcsel_endtdrifted       = particleP->T(lastPointInAV);
                            fmcsel_pxdrifted         = particleP->Px();
                            fmcsel_pydrifted         = particleP->Py();
                            fmcsel_pzdrifted         = particleP->Pz();
                            fmcsel_momentumdrifted   = particleP->P();
                            fmcsel_energydrifted     = particleP->E();
                            fmcsel_endpxdrifted      = particleP->EndPx();
                            fmcsel_endpydrifted      = particleP->EndPy();
                            fmcsel_endpzdrifted      = particleP->EndPz();
                            fmcsel_endenergydrifted  = particleP->EndE();
	                    fmcsel_pathlendrifted    = plendrifted;
	                  }
	                  fmcsel_endprocess = int(particleP->Trajectory().ProcessToKey(particleP->EndProcess()));
	                  fmcsel_theta      = particleP->Momentum().Theta();
	                  fmcsel_phi        = particleP->Momentum().Phi();
                          fmcsel_pdg        = particleP->PdgCode();
                          fmcsel_status_code= particleP->StatusCode();
                          fmcsel_mass       = particleP->Mass();
                          fmcsel_ND         = particleP->NumberDaughters();
                          fmcsel_mother     = particleP->Mother();
	                  fmcsel_origin     = mcP->Origin();
	                  fmcsel_process    = int(particleP->Trajectory().ProcessToKey(particleP->Process()));
                          fmcsel_rescatter  = particleP->Rescatter();
	
                          std::cout<<"true sel mu: Pdg code "<<fmcsel_pdg<<" track ID "<<fmcsel_trkid<<" no of daughters "<<fmcsel_ND<<" and origin_type "<<fmcsel_origin<<" fmcsel_endx "<<fmcsel_endx<<" fmcsel_endy "<<fmcsel_endy<<" fmcsel_endz "<<fmcsel_endz<<std::endl;
 
 		              fmcd_trkid      = mcd_trkid;
                              fmcd_vx         = mcd_vx;
                              fmcd_vy         = mcd_vy;
                              fmcd_vz         = mcd_vz;
                              fmcd_t          = mcd_t;
                              fmcd_endx       = mcd_endx;
                              fmcd_endy       = mcd_endy;
                              fmcd_endz       = mcd_endz;
                              fmcd_endt	      = mcd_endt;
                              fmcd_px         = mcd_px;
                              fmcd_py         = mcd_py;
                              fmcd_pz         = mcd_pz;
                              fmcd_momentum   = mcd_momentum;
                              fmcd_energy     = mcd_energy;
                              fmcd_endpx      = mcd_endpx;
                              fmcd_endpy      = mcd_endpy;
                              fmcd_endpz      = mcd_endpz;
                              fmcd_endenergy  = mcd_endenergy;
	                      fmcd_pathlen    = mcd_pathlen;
                                fmcd_vxdrifted         = mcd_vxdrifted;
                                fmcd_vydrifted         = mcd_vydrifted;
                                fmcd_vzdrifted         = mcd_vzdrifted;
                                fmcd_tdrifted          = mcd_tdrifted;
                                fmcd_endxdrifted       = mcd_endxdrifted;
                                fmcd_endydrifted       = mcd_endydrifted;
                                fmcd_endzdrifted       = mcd_endzdrifted;
                                fmcd_endtdrifted       = mcd_endtdrifted;
                                fmcd_pxdrifted         = mcd_pxdrifted;
                                fmcd_pydrifted	       = mcd_pydrifted;
                                fmcd_pzdrifted         = mcd_pzdrifted;
                                fmcd_momentumdrifted   = mcd_momentumdrifted;
                                fmcd_energydrifted     = mcd_energydrifted;
                                fmcd_endpxdrifted      = mcd_endpxdrifted;
                                fmcd_endpydrifted      = mcd_endpydrifted;
                                fmcd_endpzdrifted      = mcd_endpzdrifted;
                                fmcd_endenergydrifted  = mcd_endenergydrifted;
	                        fmcd_pathlendrifted    = mcd_pathlendrifted;
	                      fmcd_endprocess = mcd_endprocess;
	                      fmcd_theta      = mcd_theta;
	                      fmcd_phi        = mcd_phi;
                              fmcd_pdg        = mcd_pdg;
                              fmcd_status_code= mcd_status_code;
                              fmcd_mass       = mcd_mass;
                              fmcd_ND         = mcd_ND;
                              fmcd_mother     = mcd_mother;
	                      fmcd_origin     = mcd_origin;
	                      fmcd_process    = mcd_process;
                              fmcd_rescatter  = mcd_rescatter;
			      
	                  //-----------Getting the average y position of the Michel track-----------------//
                          double yposMi = 0; int countMi = 0;
			  int trueparhitscol_key[500] = {0};

 	                  std::vector<art::Ptr<recob::Hit>> hitsfromMi = bt_serv->TrackIdToHits_Ps(fmcd_trkid, hitlist);
	                  for(size_t e=0;e<hitsfromMi.size(); ++e)
	                  {			    
			    if(hitsfromMi[e]->WireID().Plane==2)
			    {
                              if(!fmsp.at(hitsfromMi[e].key()).size()) continue;
			      
			      trueparhitscol_key[countMi]   = hitsfromMi[e].key();
  			      yposMi      +=   fmsp.at(hitsfromMi[e].key())[0]->XYZ()[1];
			      countMi ++;
			    }
			  }  

 	                  std::vector<art::Ptr<recob::Hit>> hitsfromMi1 = bt_serv->TrackIdToHits_Ps(-(fmcd_trkid), hitlist);
	                  for(size_t e=0;e<hitsfromMi1.size(); ++e)
	                  {			    
			    int found_dupMi = 0;
			    if(hitsfromMi1[e]->WireID().Plane==2)
			    {
			      for(int kk=0; kk<countMi; kk++)
			      {
			        trueparhitscol_key[countMi] = hitsfromMi1[e].key();
			        if(trueparhitscol_key[kk]==trueparhitscol_key[countMi]) found_dupMi = 1;
			      }
			       
			      if(found_dupMi == 0)
			      {

                              if(!fmsp.at(hitsfromMi1[e].key()).size()) continue;
			      
  			      yposMi      +=   fmsp.at(hitsfromMi1[e].key())[0]->XYZ()[1];
			      countMi ++;
			      }
			    }
			  }  
			  
			  float avgYpos = float(yposMi)/int(countMi) ;   
			  std::cout<<"avgMichelYpos "<<avgYpos<<"\n";

	                  //-----------Tyring to solve overlapping hits issues -----------------//
			    
			    int hitmultptimeminus[20000] = {0};
			    int hitmultptimeplus[20000] = {0};
			    int prevptimeminus = 6000;
			    int prevptimeplus = 0;
			    int prevchanno = 0;
			    int countii = 0;
			    std::vector<int> mult_key; 
			    
			    int MitrackID_key[500] = {0}; int ctMi = 0;
			    
 	                    std::vector<art::Ptr<recob::Hit>> hitsfromPar = bt_serv->TrackIdToHits_Ps(fmcd_trkid, hitlist);

	                    for(size_t e=0;e<hitsfromPar.size(); ++e)
	                    {			    
			       if(hitsfromPar[e]->WireID().Plane==2 && hitsfromPar[e]->Multiplicity()>1)
			       {
			         MitrackID_key[ctMi]   = hitsfromPar[e].key();
				 ctMi++;
				 
			         int currchanno = hitsfromPar[e]->Channel();
			         if(countii==0 || currchanno !=prevchanno)
				 {
//				   cout<<"inside loop"<<"\n";
				   prevptimeminus = hitsfromPar[e]->PeakTimeMinusRMS(3);
				   prevptimeplus = hitsfromPar[e]->PeakTimePlusRMS(3);
				   prevchanno = hitsfromPar[e]->Channel();
				   hitmultptimeminus[hitsfromPar[e]->Channel()] = int(hitsfromPar[e]->PeakTimeMinusRMS(3));
				   hitmultptimeplus[hitsfromPar[e]->Channel()] = int(hitsfromPar[e]->PeakTimePlusRMS(3));
				  
				 }  
				   
				 if(hitsfromPar[e]->PeakTimeMinusRMS(3) < prevptimeminus)				           				         hitmultptimeminus[hitsfromPar[e]->Channel()] = int(hitsfromPar[e]->PeakTimeMinusRMS(3));
				 
				 if(hitsfromPar[e]->PeakTimePlusRMS(3) > prevptimeplus)
				 hitmultptimeplus[hitsfromPar[e]->Channel()] = int(hitsfromPar[e]->PeakTimePlusRMS(3));
				 
				 prevptimeminus = hitsfromPar[e]->PeakTimeMinusRMS(3);
				 prevptimeplus = hitsfromPar[e]->PeakTimePlusRMS(3);
				 prevchanno = hitsfromPar[e]->Channel();
				 
				 mult_key.push_back(hitsfromPar[e].key());
				 
//				 cout<<"hitkey "<<hitsfromPar[e].key()<<" prevptimeminus "<<prevptimeminus<<" prevptimeplus "<<prevptimeplus<<" prevchanno "<<prevchanno<<" hitmultptimeminus["<<hitsfromPar[e]->Channel()<<"] "<<hitmultptimeminus[hitsfromPar[e]->Channel()]<<" hitmultptimeplus["<<hitsfromPar[e]->Channel()<<"] "<<hitmultptimeplus[hitsfromPar[e]->Channel()]<<"\n";
				 
				 countii++;
			      }
                            }
			    
			    //Now from negative particle ID 
			    
 	                    std::vector<art::Ptr<recob::Hit>> hitsfromPar1 = bt_serv->TrackIdToHits_Ps(-(fmcd_trkid), hitlist);
			    
	                    for(size_t e=0;e<hitsfromPar1.size(); ++e)
	                    {
			       int found_dupMi = 0;		    
			       if(hitsfromPar1[e]->WireID().Plane==2 && hitsfromPar1[e]->Multiplicity()>1)
			       {
			         for(int kk=0; kk<ctMi; kk++)
			         {
				   MitrackID_key[ctMi] = hitsfromPar1[e].key();
			           if(MitrackID_key[kk]==MitrackID_key[ctMi]) found_dupMi = 1;
			         }
			       
			         if(found_dupMi == 0)
			         {
			           int currchanno = hitsfromPar1[e]->Channel();
			           if(countii==0 || currchanno !=prevchanno)
				   {
//				   cout<<"inside neg track ID loop"<<"\n";
				   prevptimeminus = hitsfromPar1[e]->PeakTimeMinusRMS(3);
				   prevptimeplus = hitsfromPar1[e]->PeakTimePlusRMS(3);
				   prevchanno = hitsfromPar1[e]->Channel();
				   hitmultptimeminus[hitsfromPar1[e]->Channel()] = int(hitsfromPar1[e]->PeakTimeMinusRMS(3));
				   hitmultptimeplus[hitsfromPar1[e]->Channel()] = int(hitsfromPar1[e]->PeakTimePlusRMS(3));
				  
//				   ctMi++;
				   }  //if(countii==0 || currchanno !=prevchanno)
				   
				 if(hitsfromPar1[e]->PeakTimeMinusRMS(3) < prevptimeminus)				           				 hitmultptimeminus[hitsfromPar1[e]->Channel()] = int(hitsfromPar1[e]->PeakTimeMinusRMS(3));
				 
				 if(hitsfromPar1[e]->PeakTimePlusRMS(3) > prevptimeplus)
				 hitmultptimeplus[hitsfromPar1[e]->Channel()] = int(hitsfromPar1[e]->PeakTimePlusRMS(3));
				 
				 prevptimeminus = hitsfromPar1[e]->PeakTimeMinusRMS(3);
				 prevptimeplus = hitsfromPar1[e]->PeakTimePlusRMS(3);
				 prevchanno = hitsfromPar1[e]->Channel();
				 
				 mult_key.push_back(hitsfromPar1[e].key());
				 
//				 cout<<"hitkey "<<hitsfromPar1[e].key()<<" prevptimeminus "<<prevptimeminus<<" prevptimeplus "<<prevptimeplus<<" prevchanno "<<prevchanno<<" hitmultptimeminus["<<hitsfromPar1[e]->Channel()<<"] "<<hitmultptimeminus[hitsfromPar1[e]->Channel()]<<" hitmultptimeplus["<<hitsfromPar1[e]->Channel()<<"] "<<hitmultptimeplus[hitsfromPar1[e]->Channel()]<<"\n";
				 
				 countii++;
			      }//if(found_dupMi == 0)
                            }//if(hitsfromPar1[e]->WireID().Plane==2...)
			  }//for(size_t e=0;e<hitsfromPar1.size(); ++e)  
			  
			  std::vector<int> mult_chan;  
			  std::vector<int> mult_ptimeminus;
			  std::vector<int> mult_ptimeplus;
			  std::vector<float> mult_trueE;
			  
			  //----Now counting energy of hits having higher multiplicity-------//
                          for(int ii = 0; ii< 20000; ii++)
                          {
			    if(hitmultptimeminus[ii]!=0)
			    {
	                      cout<<"multhit: channal "<<ii<<" ptimeminus "<<hitmultptimeminus[ii]<<" ptimeplus "<<hitmultptimeplus[ii]<<endl;
			      mult_chan.push_back(ii);
			      mult_ptimeminus.push_back(hitmultptimeminus[ii]);
			      mult_ptimeplus.push_back(hitmultptimeplus[ii]);

	                      fReadOutWindowSize = detprop->ReadOutWindowSize();
                              fNumberTimeSamples = detprop->NumberTimeSamples(); 	       
                              art::Handle< std::vector<sim::SimChannel> > simchannelHandle;
			      float totEmult = 0;	   

                              if(evt.getByLabel("largeant", simchannelHandle))
	                      {
		 
                                //Loop over simchannels 
                                for(auto const& simchannel : (*simchannelHandle))
		                {
                                  if(fGeometry->View(simchannel.Channel()) != 2) continue;	
                                  auto const& alltimeslices = simchannel.TDCIDEMap();
		   
                                  // Loop over ticks 6000
                                  for(auto const& tslice : alltimeslices)
		                  {
		                    int spill = 1;
                                    if (fReadOutWindowSize != fNumberTimeSamples) 
		                    {
                                      if (tslice.first < spill * fReadOutWindowSize || tslice.first > (spill + 1) * fReadOutWindowSize) continue;
                                    }
                                    else if (tslice.first < 0 || tslice.first > fReadOutWindowSize) continue;		 

                                    auto const& simide = tslice.second;
		                    int sim_channel = int(simchannel.Channel());
		     
                                    // Loop over energy deposits
                                    for(auto const& eDep : simide)
		                    {
			              if(sim_channel == ii && (eDep.trackID == fmcd_trkid || eDep.trackID == -fmcd_trkid))
				      if(tslice.first >= hitmultptimeminus[ii] && tslice.first <= hitmultptimeplus[ii])
			              {	
			                totEmult += eDep.energy;
                                        std::cout<<"totEmult "<<totEmult<<" tslice.first "<<tslice.first<<" MeV"<<"\n";	
			              }
				    }  //for(auto const& eDep : simide)
			          } //for(auto const& tslice : alltimeslices)  
		                }//for(auto const& simchannel : (*simchannelHandle))
                              }//if(evt.getByLabel("largeant", simchannelHandle))
			      
			      mult_trueE.push_back(totEmult);

			    }//if(hitmultptimeminus[ii]!=0)
			  }//for(int ii = 0; ii< 20000; ii++)    
			    


	                  //-----------getting a subset of the reco hits that are matched to MC particles listed in trkIDs-----------------//
			    std::vector<float> hit_key_fromHits;
			    std::vector<float> chan_charge_fromHits;
			    std::vector<int> chan_no_fromHits;
			    std::vector<float> energy_fromHits;
			    std::vector<int> ptime_minus_rms_fromHits;
			    std::vector<int> ptime_plus_rms_fromHits;
			    std::vector<float> ptime_startTick_fromHits;
			    std::vector<float> ptime_endTick_fromHits;

// 	                  std::vector<art::Ptr<recob::Hit>> hitsfromPar = bt_serv->TrackIdToHits_Ps(fmcd_trkid, hitlist);
                          std::map<int,double> mtrkide1;
	                  int mtrackid1=-1;

		            ftrueparhitallcount = 0, ftrueparhitcolcount = 0;
	                    for(size_t e=0;e<hitsfromPar.size(); ++e)
	                    { 
			     double truepardiffpeaktt0 = hitsfromPar[e]->PeakTime() - (T00/1000)*2;
		             double trueparhitXpos = detprop->ConvertTicksToX(truepardiffpeaktt0,hitsfromPar[e]->WireID().Plane,hitsfromPar[e]->WireID().TPC,hitsfromPar[e]->WireID().Cryostat); //convert ticks to usec (1tick = 0.5usec), and subtract the T0
                             double trueparWirestart[3], trueparWireend[3];
                             geometry->WireEndPoints(hitsfromPar[e]->WireID().Cryostat, hitsfromPar[e]->WireID().TPC, hitsfromPar[e]->WireID().Plane, hitsfromPar[e]->WireID().Wire, trueparWirestart, trueparWireend);
                             double trueparhitZpos = trueparWirestart[2];
//                             if(!fmsp.at(hitsfromPar[e].key()).size()) continue;
//			     double trueparhitYpos      =   fmsp.at(hitsfromPar[e].key())[0]->XYZ()[1];
                             double trueparhitYpos = 0;
                             if(fmsp.at(hitsfromPar[e].key()).size())
			     trueparhitYpos      =   fmsp.at(hitsfromPar[e].key())[0]->XYZ()[1];
			     if(!fmsp.at(hitsfromPar[e].key()).size()) trueparhitYpos = avgYpos;


	 double magnearhitvec = sqrt(pow(xposend-fsel_endhitx,2)+pow(fnearhits_zpos[int(fnearhitct/2)]-fsel_endhitz,2));
	 double maghitvec = sqrt(pow(trueparhitXpos-fsel_endhitx,2)+pow(trueparhitZpos-fsel_endhitz,2));
	 double angle_theta = acos(((xposend-fsel_endhitx)*(trueparhitXpos-fsel_endhitx) + (fnearhits_zpos[int(fnearhitct/2)]-fsel_endhitz)*(trueparhitZpos-fsel_endhitz))/(magnearhitvec*maghitvec));
	 double angle_theta_deg = (angle_theta*180)/3.14;
	 double maghitveccostheta = maghitvec* cos(angle_theta);
	 double distance = sqrt(pow(trueparhitXpos-fsel_trkendx,2)+pow(trueparhitZpos-fsel_trkendz,2));
			

		               ftrueparhitsall_key[ftrueparhitallcount]   = hitsfromPar[e].key();
		               ftrueparhitsall_peakT[ftrueparhitallcount] = hitsfromPar[e]->PeakTime();
		               ftrueparhitsall_charge[ftrueparhitallcount]= hitsfromPar[e]->Integral();
		               ftrueparhitsall_wire[ftrueparhitallcount]  = hitsfromPar[e]->WireID().Wire;
		               ftrueparhitsall_chno[ftrueparhitallcount]  = hitsfromPar[e]->Channel();
		               ftrueparhitsall_TPC[ftrueparhitallcount]   = hitsfromPar[e]->WireID().TPC;
		               ftrueparhitsall_plane[ftrueparhitallcount] = hitsfromPar[e]->WireID().Plane;
			       ftrueparhitsall_xpos[ftrueparhitallcount]  = trueparhitXpos;
			       ftrueparhitsall_ypos[ftrueparhitallcount]  = trueparhitYpos;
			       ftrueparhitsall_zpos[ftrueparhitallcount]  = trueparhitZpos;
			       ftrueparhitsall_mult[ftrueparhitallcount]      = hitsfromPar[e]->Multiplicity();
			       ftrueparhitsall_sigptime[ftrueparhitallcount]  = hitsfromPar[e]->SigmaPeakTime();
			       ftrueparhitsall_sigchrg[ftrueparhitallcount]   = hitsfromPar[e]->SigmaIntegral();
			       ftrueparhitsall_sigpamp[ftrueparhitallcount]    = hitsfromPar[e]->SigmaPeakAmplitude();
			       ftrueparhitsall_dof[ftrueparhitallcount]       = hitsfromPar[e]->DegreesOfFreedom();
			       ftrueparhitsall_gof[ftrueparhitallcount]       = hitsfromPar[e]->GoodnessOfFit();
		               ftrueparhitsall_ptminusRMS[ftrueparhitallcount]	  = hitsfromPar[e]->PeakTimeMinusRMS(5.0);
		               ftrueparhitsall_ptplusRMS[ftrueparhitallcount]	 = hitsfromPar[e]->PeakTimePlusRMS(5.0);

			       ftrueparhitallcount++;
			       
			       if(hitsfromPar[e]->WireID().Plane==2)// && hitsfromPar[e]->Multiplicity()<2)
			       {
			    
			    //--------look at true info of this hit--------//
			    
			    art::Ptr<recob::Hit> mhit1=hitlist[hitsfromPar[e].key()];
	                    std::vector<sim::TrackIDE> meveIDs1 = bt_serv->HitToTrackIDEs(mhit1);
                            
			    cout<<"\n";
//			    cout<<"key "<<hitsfromPar[e].key()<<" angle_theta_deg "<<angle_theta_deg<<" maghitveccostheta "<<maghitveccostheta<<" mult "<<hitsfromPar[e]->Multiplicity()<<" charge "<<hitsfromPar[e]->Integral()<<"\n";
			    float maxen1 = -1; int maxenid1  =0;  float maxenfrac1=-1;
			    std::map<int,int> mtrkid11;
			    std::map<int,double> mtrkide11;
			    std::map<int,double> mtrkidefrac11;			    

	                    for(size_t e=0;e<meveIDs1.size(); ++e)
	                    { 
//	                        mtrkide1[abs(meveIDs1[e].trackID)] += meveIDs1[e].energy;
				mtrkid11[abs(meveIDs1[e].trackID)] = abs(meveIDs1[e].trackID);
	                        mtrkide11[abs(meveIDs1[e].trackID)] += meveIDs1[e].energy;
				mtrkidefrac11[abs(meveIDs1[e].trackID)] += meveIDs1[e].energyFrac;			
		                cout<<"trackID "<<meveIDs1[e].trackID<<" energy "<<meveIDs1[e].energy<<" energyFrac "<<meveIDs1[e].energyFrac<<"\n";
				if(mtrkide11[abs(meveIDs1[e].trackID)]>maxen1)
				{
				  maxen1 = mtrkide11[abs(meveIDs1[e].trackID)];
				  maxenid1 = abs(meveIDs1[e].trackID);
				  maxenfrac1 = mtrkidefrac11[abs(meveIDs1[e].trackID)];
				  cout<<"maxen1 "<<maxen1<<" maxenid1 "<<maxenid1<<" RealID "<<meveIDs1[e].trackID<<" maxenfrac1 "<<maxenfrac1<<"\n";
				}
	                    }
			    
                          //see what other particles take away the energy fraction;
                          for(std::map<int,double>::iterator ii = mtrkidefrac11.begin(); ii!=mtrkidefrac11.end(); ++ii)
                          {              
	                    const simb::MCParticle *mparticleP1 = pi_serv->TrackIdToParticle_P(abs(ii->first));
                            // if(!mparticleP) continue;

                            // const art::Ptr<simb::MCTruth> mmcP=pi_serv->TrackIdToMCTruth_P(mtrackid);
                            // if(!mmcP) continue;

                            cout<<"TrackID "<<ii->first<<" energyFrac "<<ii->second<<" pdg "<<mparticleP1->PdgCode()<<" E "<<mparticleP1->E()<<"\n";
                          }
			  
			  ftrueMiEFrac[ftrueparhitcolcount] = mtrkidefrac11[abs(fmcd_trkid)];
			  
 
			    if(hitsfromPar[e]->Multiplicity()<2)
			    {
			    //for counting true MC energy		    
	                    for(size_t e=0;e<meveIDs1.size(); ++e)
	                    { 
	                        mtrkide1[abs(meveIDs1[e].trackID)] += meveIDs1[e].energy;
//		                cout<<"abs(meveIDs1[e].trackID) "<<abs(meveIDs1[e].trackID)<<" mtrkide1[abs(meveIDs1[e].trackID)] "<<mtrkide1[abs(meveIDs1[e].trackID)]<<"\n";
			    }	
			    ////////////////////////////////
                               hit_key_fromHits.push_back(hitsfromPar[e].key());
                               chan_charge_fromHits.push_back(hitsfromPar[e]->Integral());
                               chan_no_fromHits.push_back(hitsfromPar[e]->Channel());
			       energy_fromHits.push_back(maxen1);
			       ptime_minus_rms_fromHits.push_back(hitsfromPar[e]->PeakTimeMinusRMS(3));
			       ptime_plus_rms_fromHits.push_back(hitsfromPar[e]->PeakTimePlusRMS(3));
			       ptime_startTick_fromHits.push_back(hitsfromPar[e]->StartTick());
			       ptime_endTick_fromHits.push_back(hitsfromPar[e]->EndTick());
			    }//if(hitsfromPar[e]->Multiplicity()<2)   

                            if(abs(maxenid1) != abs(fmcd_trkid)) continue; // to make sure that the maximum energy deposition is coming from the true Michel 
			    
		               ftrueparhitscol_key[ftrueparhitcolcount]   = hitsfromPar[e].key();
		               ftrueparhitscol_peakT[ftrueparhitcolcount] = hitsfromPar[e]->PeakTime();
		               ftrueparhitscol_charge[ftrueparhitcolcount]= hitsfromPar[e]->Integral();
		               ftrueparhitscol_wire[ftrueparhitcolcount]  = hitsfromPar[e]->WireID().Wire;
		               ftrueparhitscol_chno[ftrueparhitcolcount]  = hitsfromPar[e]->Channel();
		               ftrueparhitscol_TPC[ftrueparhitcolcount]   = hitsfromPar[e]->WireID().TPC;
		               ftrueparhitscol_plane[ftrueparhitcolcount] = hitsfromPar[e]->WireID().Plane;
			       ftrueparhitscol_xpos[ftrueparhitcolcount]  = trueparhitXpos;
			       ftrueparhitscol_ypos[ftrueparhitcolcount]  = trueparhitYpos;
			       ftrueparhitscol_zpos[ftrueparhitcolcount]  = trueparhitZpos;
	                       ftrueparhitscol_angledeg[ftrueparhitcolcount] = angle_theta_deg;
	                       ftrueparhitscol_maghitveccostheta[ftrueparhitcolcount] = maghitveccostheta;
	                       ftrueparhitscol_distance[ftrueparhitcolcount] = distance;
			       ftrueparhitscol_mult[ftrueparhitcolcount]      = hitsfromPar[e]->Multiplicity();
			       ftrueparhitscol_sigptime[ftrueparhitcolcount]  = hitsfromPar[e]->SigmaPeakTime();
			       ftrueparhitscol_sigchrg[ftrueparhitcolcount]   = hitsfromPar[e]->SigmaIntegral();
			       ftrueparhitscol_sigpamp[ftrueparhitcolcount]    = hitsfromPar[e]->SigmaPeakAmplitude();
			       ftrueparhitscol_dof[ftrueparhitcolcount]       = hitsfromPar[e]->DegreesOfFreedom();
			       ftrueparhitscol_gof[ftrueparhitcolcount]       = hitsfromPar[e]->GoodnessOfFit();
		               ftrueparhitscol_ptminusRMS[ftrueparhitcolcount]     = hitsfromPar[e]->PeakTimeMinusRMS(5.0);
		               ftrueparhitscol_ptplusRMS[ftrueparhitcolcount]     = hitsfromPar[e]->PeakTimePlusRMS(5.0);
                       	       
			       trueHitsKey.push_back(hitsfromPar[e].key());
     
			       ftrueparhitcolcount++;
			       }
//                               std::cout<<"hit key from true par "<<hitsfromPar[e].key()<<"\n";
//			     } //if(same_trk == 0) 
                           }

                            cout<<"\n Now from -(fmcd_trkid)"<<"\n";
// 	                    std::vector<art::Ptr<recob::Hit>> hitsfromPar1 = bt_serv->TrackIdToHits_Ps(-(fmcd_trkid), hitlist);
	                    for(size_t e=0;e<hitsfromPar1.size(); ++e)
	                    { 

			     double truepardiffpeaktt0 = hitsfromPar1[e]->PeakTime() - (T00/1000)*2;
		             double trueparhitXpos = detprop->ConvertTicksToX(truepardiffpeaktt0,hitsfromPar1[e]->WireID().Plane,hitsfromPar1[e]->WireID().TPC,hitsfromPar1[e]->WireID().Cryostat); //convert ticks to usec (1tick = 0.5usec), and subtract the T0
                             double trueparWirestart[3], trueparWireend[3];
                             geometry->WireEndPoints(hitsfromPar1[e]->WireID().Cryostat, hitsfromPar1[e]->WireID().TPC, hitsfromPar1[e]->WireID().Plane, hitsfromPar1[e]->WireID().Wire, trueparWirestart, trueparWireend);
                             double trueparhitZpos = trueparWirestart[2];
//                             if(!fmsp.at(hitsfromPar1[e].key()).size()) continue;
//			     double trueparhitYpos      =   fmsp.at(hitsfromPar1[e].key())[0]->XYZ()[1];
                             double trueparhitYpos = 0;
                             if(fmsp.at(hitsfromPar1[e].key()).size())
			     trueparhitYpos      =   fmsp.at(hitsfromPar1[e].key())[0]->XYZ()[1];
			     if(!fmsp.at(hitsfromPar1[e].key()).size()) trueparhitYpos = avgYpos;


	 double magnearhitvec = sqrt(pow(xposend-fsel_endhitx,2)+pow(fnearhits_zpos[int(fnearhitct/2)]-fsel_endhitz,2));
	 double maghitvec = sqrt(pow(trueparhitXpos-fsel_endhitx,2)+pow(trueparhitZpos-fsel_endhitz,2));
	 double angle_theta = acos(((xposend-fsel_endhitx)*(trueparhitXpos-fsel_endhitx) + (fnearhits_zpos[int(fnearhitct/2)]-fsel_endhitz)*(trueparhitZpos-fsel_endhitz))/(magnearhitvec*maghitvec));
	 double angle_theta_deg = (angle_theta*180)/3.14;
	 double maghitveccostheta = maghitvec* cos(angle_theta);
	 double distance = sqrt(pow(trueparhitXpos-fsel_trkendx,2)+pow(trueparhitZpos-fsel_trkendz,2));
			

		               ftrueparhitsall_key[ftrueparhitallcount]   = hitsfromPar1[e].key();

			       int found_dup = 0;
			       for(int kk=0; kk<ftrueparhitallcount; kk++)
			       {
			         if(ftrueparhitsall_key[ftrueparhitallcount]==ftrueparhitsall_key[kk]) found_dup = 1;
			       }
			       
			       if(found_dup == 0)
			       {
		               ftrueparhitsall_peakT[ftrueparhitallcount] = hitsfromPar1[e]->PeakTime();
		               ftrueparhitsall_charge[ftrueparhitallcount]= hitsfromPar1[e]->Integral();
		               ftrueparhitsall_wire[ftrueparhitallcount]  = hitsfromPar1[e]->WireID().Wire;
		               ftrueparhitsall_chno[ftrueparhitallcount]  = hitsfromPar1[e]->Channel();
		               ftrueparhitsall_TPC[ftrueparhitallcount]   = hitsfromPar1[e]->WireID().TPC;
		               ftrueparhitsall_plane[ftrueparhitallcount] = hitsfromPar1[e]->WireID().Plane;
			       ftrueparhitsall_xpos[ftrueparhitallcount]  = trueparhitXpos;
			       ftrueparhitsall_ypos[ftrueparhitallcount]  = trueparhitYpos;
			       ftrueparhitsall_zpos[ftrueparhitallcount]  = trueparhitZpos;
			       ftrueparhitsall_mult[ftrueparhitallcount]      = hitsfromPar1[e]->Multiplicity();
			       ftrueparhitsall_sigptime[ftrueparhitallcount]  = hitsfromPar1[e]->SigmaPeakTime();
			       ftrueparhitsall_sigchrg[ftrueparhitallcount]   = hitsfromPar1[e]->SigmaIntegral();
			       ftrueparhitsall_sigpamp[ftrueparhitallcount]    = hitsfromPar1[e]->SigmaPeakAmplitude();
			       ftrueparhitsall_dof[ftrueparhitallcount]       = hitsfromPar1[e]->DegreesOfFreedom();
			       ftrueparhitsall_gof[ftrueparhitallcount]       = hitsfromPar1[e]->GoodnessOfFit();
		               ftrueparhitsall_ptminusRMS[ftrueparhitallcount]	  = hitsfromPar1[e]->PeakTimeMinusRMS(5.0);
		               ftrueparhitsall_ptplusRMS[ftrueparhitallcount]	 = hitsfromPar1[e]->PeakTimePlusRMS(5.0);
			       

			       ftrueparhitallcount++;
			       }


			       if(hitsfromPar1[e]->WireID().Plane==2)// && hitsfromPar1[e]->Multiplicity()<2)
			       {
		               ftrueparhitscol_key[ftrueparhitcolcount]   = hitsfromPar1[e].key();
			       int found_dup1 = 0;
			       for(int kk=0; kk<ftrueparhitcolcount; kk++)
			       {
			         if(ftrueparhitscol_key[ftrueparhitcolcount]==ftrueparhitscol_key[kk]) found_dup1 = 1;
			       }
			       
			       if(found_dup1 == 0)
			       {

			    
			    //--------look at true info of this hit--------//
			    
			    art::Ptr<recob::Hit> mhit2=hitlist[hitsfromPar1[e].key()];
	                    std::vector<sim::TrackIDE> meveIDs2 = bt_serv->HitToTrackIDEs(mhit2);

                            cout<<"\n";
//			    cout<<"key "<<hitsfromPar1[e].key()<<" angle_theta_deg "<<angle_theta_deg<<" maghitveccostheta "<<maghitveccostheta<<" mult "<<hitsfromPar1[e]->Multiplicity()<<" charge "<<hitsfromPar1[e]->Integral()<<"\n";
			    			    
			    float maxen2 = -1; int maxenid2  =0; float maxenfrac2=-1;
			    std::map<int,double> mtrkide12;
			    std::map<int,double> mtrkidefrac12;

	                    for(size_t e=0;e<meveIDs2.size(); ++e)
	                    { 
//	                        mtrkide1[abs(meveIDs2[e].trackID)] += meveIDs2[e].energy;
	                        mtrkide12[abs(meveIDs2[e].trackID)] += meveIDs2[e].energy;
				mtrkidefrac12[abs(meveIDs2[e].trackID)] += meveIDs2[e].energyFrac;			
		                cout<<"trackID "<<meveIDs2[e].trackID<<" energy "<<meveIDs2[e].energy<<" energyFrac "<<meveIDs2[e].energyFrac<<"\n";
//		                cout<<"abs(meveIDs2[e].trackID) "<<abs(meveIDs2[e].trackID)<<" mtrkide12[abs(meveIDs2[e].trackID)] "<<mtrkide12[abs(meveIDs2[e].trackID)]<<"\n";
				if(mtrkide12[abs(meveIDs2[e].trackID)]>maxen2)
				{
				  maxen2 = mtrkide12[abs(meveIDs2[e].trackID)];
				  maxenid2 = abs(meveIDs2[e].trackID);
				  maxenfrac2 = mtrkidefrac12[abs(meveIDs2[e].trackID)];
				  cout<<"maxen2 "<<maxen2<<" maxenid2 "<<maxenid2<<" RealID "<<meveIDs2[e].trackID<<" maxenfrac2 "<<maxenfrac2<<"\n";
				}  
	                    }

                          //see what other particles take away the energy fraction;
                          for(std::map<int,double>::iterator ii = mtrkidefrac12.begin(); ii!=mtrkidefrac12.end(); ++ii)
                          {              
	                    const simb::MCParticle *mparticleP2 = pi_serv->TrackIdToParticle_P(abs(ii->first));
                            // if(!mparticleP) continue;

              
	                    // const art::Ptr<simb::MCTruth> mmcP=pi_serv->TrackIdToMCTruth_P(mtrackid);
                            // if(!mmcP) continue;

                            cout<<"TrackID "<<ii->first<<" energyFrac "<<ii->second<<" pdg "<<mparticleP2->PdgCode()<<" E "<<mparticleP2->E()<<"\n";
                          }

			     ftrueMiEFrac[ftrueparhitcolcount] = mtrkidefrac12[abs(fmcd_trkid)];

			    if(hitsfromPar1[e]->Multiplicity()<2)
			    {
                            //for counting true MC energy
	                    for(size_t e=0;e<meveIDs2.size(); ++e)
	                    { 
	                        mtrkide1[abs(meveIDs2[e].trackID)] += meveIDs2[e].energy;
//		                cout<<"abs(meveIDs2[e].trackID) "<<abs(meveIDs2[e].trackID)<<" mtrkide1[abs(meveIDs2[e].trackID)] "<<mtrkide12[abs(meveIDs2[e].trackID)]<<"\n";
			    }
			    ////////////////////////////////
			       hit_key_fromHits.push_back(hitsfromPar1[e].key());
                               chan_charge_fromHits.push_back(hitsfromPar1[e]->Integral());
                               chan_no_fromHits.push_back(hitsfromPar1[e]->Channel());
			       energy_fromHits.push_back(maxen2);
			       ptime_minus_rms_fromHits.push_back(hitsfromPar1[e]->PeakTimeMinusRMS(3));
			       ptime_plus_rms_fromHits.push_back(hitsfromPar1[e]->PeakTimePlusRMS(3));
			       ptime_startTick_fromHits.push_back(hitsfromPar1[e]->StartTick());
			       ptime_endTick_fromHits.push_back(hitsfromPar1[e]->EndTick());
			    } // if(hitsfromPar1[e]->Multiplicity()<2)  

                            if(abs(maxenid2) != abs(fmcd_trkid)) continue; // to make sure that the maximum energy deposition is coming from the true Michel 
		               ftrueparhitscol_peakT[ftrueparhitcolcount] = hitsfromPar1[e]->PeakTime();
		               ftrueparhitscol_charge[ftrueparhitcolcount]= hitsfromPar1[e]->Integral();
		               ftrueparhitscol_wire[ftrueparhitcolcount]  = hitsfromPar1[e]->WireID().Wire;
		               ftrueparhitscol_chno[ftrueparhitcolcount]  = hitsfromPar1[e]->Channel();
		               ftrueparhitscol_TPC[ftrueparhitcolcount]   = hitsfromPar1[e]->WireID().TPC;
		               ftrueparhitscol_plane[ftrueparhitcolcount] = hitsfromPar1[e]->WireID().Plane;
			       ftrueparhitscol_xpos[ftrueparhitcolcount]  = trueparhitXpos;
			       ftrueparhitscol_ypos[ftrueparhitcolcount]  = trueparhitYpos;
			       ftrueparhitscol_zpos[ftrueparhitcolcount]  = trueparhitZpos;
	                       ftrueparhitscol_angledeg[ftrueparhitcolcount] = angle_theta_deg;
	                       ftrueparhitscol_maghitveccostheta[ftrueparhitcolcount] = maghitveccostheta;
	                       ftrueparhitscol_distance[ftrueparhitcolcount] = distance;
			       ftrueparhitscol_mult[ftrueparhitcolcount]      = hitsfromPar1[e]->Multiplicity();
			       ftrueparhitscol_sigptime[ftrueparhitcolcount]  = hitsfromPar1[e]->SigmaPeakTime();
			       ftrueparhitscol_sigchrg[ftrueparhitcolcount]   = hitsfromPar1[e]->SigmaIntegral();
			       ftrueparhitscol_sigpamp[ftrueparhitcolcount]    = hitsfromPar1[e]->SigmaPeakAmplitude();
			       ftrueparhitscol_dof[ftrueparhitcolcount]       = hitsfromPar1[e]->DegreesOfFreedom();
			       ftrueparhitscol_gof[ftrueparhitcolcount]       = hitsfromPar1[e]->GoodnessOfFit();
		               ftrueparhitscol_ptminusRMS[ftrueparhitcolcount]     = hitsfromPar1[e]->PeakTimeMinusRMS(5.0);
		               ftrueparhitscol_ptplusRMS[ftrueparhitcolcount]     = hitsfromPar1[e]->PeakTimePlusRMS(5.0);
                       	       
			       trueHitsKey.push_back(hitsfromPar1[e].key());
     
			       ftrueparhitcolcount++;
			       }
			       }//if(hitsfromPar1[e]->WireID().Plane==2)
//                               std::cout<<"hit key from true par "<<hitsfromPar[e].key()<<"\n";
//			     } //if(same_trk == 0) 
                           }//for(size_t e=0;e<hitsfromPar1.size(); ++e)
			   cout<<"ftrueparhitcolcount "<<ftrueparhitcolcount<<"\n";

                          double  mmaxe1 = -1;
                          //double mtote1 = 0;
                          for(std::map<int,double>::iterator ii = mtrkide1.begin(); ii!=mtrkide1.end(); ++ii)
                          {
			    if(abs(ii->first) == fmcd_trkid)
			    {
	                    	cout<<"\ntrkid "<<ii->first<<" energy deposited from hitMult=1 "<<ii->second<<endl;
//	                    if((ii->second)>mmaxe1)
//	                    {
	                      mmaxe1 = ii->second;
	                      mtrackid1 = abs(ii->first);
//	                    }
                            }
                          }
			  
	                 //--------Total Michel energy deposited on hits including hit multiplicity >= 2--------------//
                         double tothitmultE = 0;
	                 for(size_t nn=0; nn<mult_chan.size(); nn++)
	                 {
			   tothitmultE += mult_trueE.at(nn);
                           chan_no_fromHits.push_back(mult_chan.at(nn));
			   ptime_minus_rms_fromHits.push_back(mult_ptimeminus.at(nn)); 
			   ptime_plus_rms_fromHits.push_back(mult_ptimeplus.at(nn)); 
			   energy_fromHits.push_back(mult_trueE.at(nn)); 
			   cout<<"tothitmultE "<<tothitmultE<<"\n";
	                 }
			 
	                 cout<<" Total true Michel energy deposited on hits "<<(tothitmultE+mmaxe1)<<"\n";

                         const simb::MCParticle *mparticleP1 = pi_serv->TrackIdToParticle_P(mtrackid1);
                        // if(!mparticleP) continue;

                        // const art::Ptr<simb::MCTruth> mmcP=pi_serv->TrackIdToMCTruth_P(mtrackid);
                        // if(!mmcP) continue;

 cout<<"mtrackid1 "<<mtrackid1<<" energy "<<(tothitmultE+mmaxe1)<<" pdg "<<mparticleP1->PdgCode()<<" E "<<mparticleP1->E()<<"\n";
 
                              fmcd_trueselhitsE     = (tothitmultE+mmaxe1);
			      ftrueEdepo = truthUtil.GetDepEnergyMC(evt, fGeometry, mtrackid1, 2 );
			      cout<<"pdg "<<mparticleP1->PdgCode()<<" trueEdepo "<<ftrueEdepo<<"\n";
			   
	                  //-----------getting the MC true information of the selected nearby michel hits-----------------//
	       
                          std::map<int,double> mtrkide;
	                  int mtrackid=-1;
                          for(int qq =0; qq < fnearhitct; qq++)
		          {
		            //cout<<"hitKey[qq] "<<hitKey[qq]<<"\n";
		            art::Ptr<recob::Hit> mhit=hitlist[fnearhits_key[qq]];    
		            cout<<"michel hit key "<<mhit.key()<<" peakT "<<mhit->PeakTime()<<" WireID "<<mhit->WireID().Wire<<"\n";
	                    std::vector<sim::TrackIDE> meveIDs = bt_serv->HitToTrackIDEs(mhit);

	                    for(size_t e=0;e<meveIDs.size(); ++e)
	                    { 
		              if(abs(meveIDs[e].trackID)!=fmcsel_trkid) //the hits should not belong to the same mother track (muon)
		              {
	                        mtrkide[meveIDs[e].trackID] += meveIDs[e].energy;
		                cout<<"meveIDs[e].trackID "<<meveIDs[e].trackID<<" mtrkide[meveIDs[e].trackID] "<<mtrkide[meveIDs[e].trackID]<<"\n";
		              }
	                    }
                          }
                          double  mmaxe = -1;
                          double mtote = 0;
                          for(std::map<int,double>::iterator ii = mtrkide.begin(); ii!=mtrkide.end(); ++ii)
                          {
	                    //	cout<<" trkid "<<ii->first<<"  energy deposited = "<<ii->second<<endl;
	                    mtote += ii->second;
	                    if((ii->second)>mmaxe)
	                    {
	                      mmaxe = ii->second;
	                      mtrackid = abs(ii->first);
	                    }
                          }
/*                        float mtotal_energy=0.0;
                          for(int qq =0; qq < hitcount; qq++)
		          {
	                    art::Ptr<recob::Hit> hit=hitlist[hitKey.at(qq)];   
	                    if (hit->WireID().Plane!=2) continue;
	                    std::vector<sim::TrackIDE> eveIDs = bt_serv->HitToTrackIDEs(hit);
	                    for(size_t e=0;e<eveIDs.size(); ++e)
		           {
	                     if (eveIDs[e].trackID == mtrackid) mtotal_energy +=  eveIDs[e].energy;
	                   }
                         }
*/ 
                         const simb::MCParticle *mparticleP = pi_serv->TrackIdToParticle_P(mtrackid);
                         if(!mparticleP) continue;

                         const art::Ptr<simb::MCTruth> mmcP=pi_serv->TrackIdToMCTruth_P(mtrackid);
                         if(!mmcP) continue;

	                 TLorentzVector mmcstart, mmcend, mmcstartdrifted, mmcenddrifted;
	                 unsigned int mpstarti, mpendi, mpstartdriftedi, mpenddriftedi; //mcparticle indices for starts and ends in tpc or drifted volumes
	                 double mplen = length(*mparticleP, mmcstart, mmcend, mpstarti, mpendi);
	                 double mplendrifted = driftedLength(*mparticleP, mmcstartdrifted, mmcenddrifted, mpstartdriftedi, mpenddriftedi);
	                 // bool isActive = plen != 0;
	                 bool misDrifted = mplendrifted!= 0;

                         cout<<"Michel hits: bt:trackid "<<mtrackid<<"\n";
                         fmchits_trkid      = mparticleP->TrackId();
                         fmchits_vx         = mparticleP->Vx();
                         fmchits_vy         = mparticleP->Vy();
                         fmchits_vz         = mparticleP->Vz();
                         fmchits_t          = mparticleP->T();
                         fmchits_endx       = mparticleP->EndX();
                         fmchits_endy       = mparticleP->EndY();
                         fmchits_endz       = mparticleP->EndZ();
                         fmchits_endt       = mparticleP->EndT();
                         fmchits_px         = mparticleP->Px();
                         fmchits_py         = mparticleP->Py();
                         fmchits_pz         = mparticleP->Pz();
                         fmchits_momentum   = mparticleP->P();
                         fmchits_energy     = mparticleP->E();
                         fmchits_endpx      = mparticleP->EndPx();
                         fmchits_endpy      = mparticleP->EndPy();
                         fmchits_endpz      = mparticleP->EndPz();
                         fmchits_endenergy  = mparticleP->EndE();
	                 fmchits_pathlen    = mplen;
	                 if (misDrifted)
	                 {
                           fmchits_vxdrifted         = mparticleP->Vx();
                           fmchits_vydrifted         = mparticleP->Vy();
                           fmchits_vzdrifted         = mparticleP->Vz();
                           fmchits_tdrifted          = mparticleP->T();
                           fmchits_endxdrifted       = mparticleP->EndX();
                           fmchits_endydrifted       = mparticleP->EndY();
                           fmchits_endzdrifted       = mparticleP->EndZ();
                           fmchits_endtdrifted       = mparticleP->EndT();
                           fmchits_pxdrifted         = mparticleP->Px();
                           fmchits_pydrifted         = mparticleP->Py();
                           fmchits_pzdrifted         = mparticleP->Pz();
                           fmchits_momentumdrifted   = mparticleP->P();
                           fmchits_energydrifted     = mparticleP->E();
                           fmchits_endpxdrifted      = mparticleP->EndPx();
                           fmchits_endpydrifted      = mparticleP->EndPy();
                           fmchits_endpzdrifted      = mparticleP->EndPz();
                           fmchits_endenergydrifted  = mparticleP->EndE();
	                   fmchits_pathlendrifted    = mplendrifted;
	                 }
	                 fmchits_endprocess = int(mparticleP->Trajectory().ProcessToKey(mparticleP->EndProcess()));
	                 fmchits_theta      = mparticleP->Momentum().Theta();
	                 fmchits_phi        = mparticleP->Momentum().Phi();
                         fmchits_pdg        = mparticleP->PdgCode();
                         fmchits_status_code= mparticleP->StatusCode();
                         fmchits_mass       = mparticleP->Mass();
                         fmchits_ND         = mparticleP->NumberDaughters();
                         fmchits_mother     = mparticleP->Mother();
	                 fmchits_origin     = mmcP->Origin();
	                 fmchits_process    = int(mparticleP->Trajectory().ProcessToKey(mparticleP->Process()));
                         fmchits_rescatter  = mparticleP->Rescatter();
	
                         std::cout<<"MC true Michel hits: Pdg code "<<fmchits_pdg<<" track ID "<<fmchits_trkid<<" no of daughters "<<fmchits_ND<<" and origin_type "<<fmchits_origin<<std::endl;


	   //--------Calculate the average missing energy-----------//

	       fReadOutWindowSize = detprop->ReadOutWindowSize();
               fNumberTimeSamples = detprop->NumberTimeSamples();
	       std::cout<<"fReadOutWindowSize "<<fReadOutWindowSize<<" fNumberTimeSamples "<<fNumberTimeSamples<<"\n";

               double totEMiDepo =0;
	       double totNelectrons = 0;
	       float avg_missing_energy = 0, avg_missing_numelec=0;
	       favg_missing_energy = 0, favg_missing_numelec = 0;
	       
               art::Handle< std::vector<sim::SimChannel> > simchannelHandle;
               if(evt.getByLabel("largeant", simchannelHandle))
	       {
		 
                 //Loop over simchannels 
                 for(auto const& simchannel : (*simchannelHandle))
		 {
		 
                   if(fGeometry->View(simchannel.Channel()) != 2) continue;
		   
                   auto const& alltimeslices = simchannel.TDCIDEMap();
		   
//                   bool savethis = false;
		   
//	           std::cout<<"Simchan_charge["<<simchannel.Channel()<<"] "<<chan_charge[simchannel.Channel()]<<"\n";
		   
//                   if( chan_charge[simchannel.Channel()] ==  0) savethis = true;
		   
		   double MichelEperCh = 0;
		   double MissingEperCh = 0;
		   
                   // Loop over ticks 60000
                   for(auto const& tslice : alltimeslices)
		   {
		     int spill = 1;
                     if (fReadOutWindowSize != fNumberTimeSamples) 
		     {
                       if (tslice.first < spill * fReadOutWindowSize || tslice.first > (spill + 1) * fReadOutWindowSize) continue;
                     }
                     else if (tslice.first < 0 || tslice.first > fReadOutWindowSize) continue;		 

                     auto const& simide = tslice.second;
		     float simETick = 0;
		     

		     //---For calculation of missing energy, only look for TDCs where we dont see a signal------//
		     bool missed_signal=true;
                     int sim_channel = int(simchannel.Channel());
//		     std::cout<<"chan_charge_fromHits.size() "<<chan_charge_fromHits.size()<<"\n";

		     //---------- recob:hit signal--------------//
		     for(size_t nn=0; nn<chan_no_fromHits.size(); nn++)
		     {
		       if(chan_no_fromHits.at(nn) == sim_channel && ptime_minus_rms_fromHits.at(nn)<=tslice.first && ptime_plus_rms_fromHits.at(nn)>=tslice.first) missed_signal=false;
		       
		     }
//		     cout<<"TDC "<<tslice.first<<" charge_on_tick "<<simchannel.Charge(tslice.first)<<" energy_on_tick "<<simchannel.Energy(tslice.first)<<" MeV, missed_signal "<<missed_signal<<"\n";
		      
//		     if(missed_signal)
//		     {
//		       fmisssimETick = (simchannel.Energy(tslice.first));
//		       fTickTree->Fill();
//		     }

                     // Loop over energy deposits
                     for(auto const& eDep : simide)
		     {

		       if(eDep.trackID == fmcd_trkid || eDep.trackID == -fmcd_trkid)
		       {
		         MichelEperCh += eDep.energy;
			 totEMiDepo += eDep.energy;
			 totNelectrons += eDep.numElectrons;
//			 std::cout<<"eDep.energy "<<eDep.energy<<"\n";
//                         if( savethis ) 
                         if( missed_signal ) 
			 {
			   simETick += eDep.energy;
			   avg_missing_energy += eDep.energy;
			   avg_missing_numelec += eDep.numElectrons;
			   MissingEperCh += eDep.energy;
			   			   
//			   std::cout<<"eDep.energy "<<eDep.energy<<" avg_missing_energy "<<avg_missing_energy<<" MeV\n";
			 }  
		       }
                     }//for(auto const& eDep : simide
		     
/*		     //-----Count the missing Charge-------//
                     auto const& wires = evt.getValidHandle<std::vector<recob::Wire> >(fWireModuleLabel);	       
//		     std::vector<int> chan_fromWires;
//		     std::vector<float> signal_fromWires;
                     for(auto & wire : * wires)
	             {
                       int wire_channel_no = wire.Channel();
                       int wire_plane = fGeometry->View(wire.Channel()); 
                       if( wire_plane != 2) continue;
		       
//		       chan_fromWires.push_back(wire_channel_no);

//			 float totwsig = 0;
                        if(wire_channel_no==sim_channel && i<=tslice.first && i>=tslice.first)
			{
			  missed_signal2 = false;
                         for(size_t i = 0; i < wire.Signal().size(); ++i)
	                 {
			   totwsig += wire.Signal()[i];	
                         } 
//			 signal_fromWires.push_back(totwsig); 
			} 
	              } // for(auto & wire : * wires)  
		     
*/		     
		     //cout<<"missing simE/tick "<<simETick<<" MeV\n";
		     //fmisssimETick->Fill(simETick);
                   }//for(auto const& tslice : alltimeslices)
//		   cout<<"MichelEperCh["<<simchannel.Channel()<<"] "<<MichelEperCh<<" MissingEperCh "<<MissingEperCh<<"\n";
//		   int simchan = simchannel.Channel();
//		   ftotEnergyperChan->Fill(simchan,MichelEperCh);
//		   fmissingEnergyperChan->Fill(simchan,MissingEperCh);
                 } //for(auto const& simchannel : (*simchannelHandle))  
               }//if(evt.getByLabel("largeant", simchannelHandle))
	       
		cout<<"chan_no_fromHits.size() "<<chan_no_fromHits.size()<<" chan_charge_fromHits.size() "<<chan_charge_fromHits.size()<<"\n";

	       favg_missing_energy = avg_missing_energy;
	       favg_missing_numelec = avg_missing_numelec;
	       cout<<"\ntotEMiDepo "<<totEMiDepo<<" avg_missing_E "<<favg_missing_energy<<" MeV\n";
	       cout<<"avg_missing_numelec "<<favg_missing_numelec<<"\n";

 cout<<"\nGEANTMiE "<<fmcd_energy*1000<<" totEMiDepo "<<totEMiDepo<<" fmcd_trueselhitsE "<<fmcd_trueselhitsE<<" avg_missing_E "<<favg_missing_energy<<" sel+missingE "<<(fmcd_trueselhitsE+favg_missing_energy)<<" MeV\n"<<"\n";
 
     	        } //if(isMC)  
 	     


                        //------------Making a cone to see which other hits might belong to a michel (2nd method)--------------//
			
			
//			double magnearhitvec = sqrt(pow(fnearhits_xpos[fnearhitct-1]-endhitx,2)+pow(fnearhits_zpos[fnearhitct-1]-endhitz,2));


 /*      //True michel reco hits 
       std::vector <double> trueHitsKey;
     
       cout<<"\nftrueparhitcolcount "<<ftrueparhitcolcount<<"\n";
       for(int kk=0; kk<ftrueparhitcolcount; kk++)
       {
	 trueHitsKey.push_back(ftrueparhitscol_key[kk]);
       } 
			
*/

                         //------------Looping over all hits and checking if they belong to the michel cone -------------//
		         fmhitcount = 0;
			 			                 		
                         for(size_t ll=0; ll<hitlist.size();++ll)
                         {
/*			   if(thass.isValid())
			   {
			     if(thass.at(ll).size()!=0)
			     {
			       int hit_trkid = thass.at(ll)[0]->ID();
			       double hit_trklen = thass.at(ll)[0]->Length();
			   //auto tracks = thass.at(hitlist[ll].key());
*/	 	 	   if(hitlist[ll]->WireID().Plane==2)
	                   {
			   
	                     //----Make sure that this hit does not belong to the  candidate muon track----//
                             int same_trk = 0;
	                     for(size_t mm=0;mm<longtrk_hitkey.size();mm++)
	                     {
			       if(longtrk_hitkey.at(mm)==hitlist[ll].key())
		               { 
                                 same_trk = 1;
		               }
                             }
		
                             //---To make sure that these hits dont belong to another long track-----//
				  
        	             int long_trk = 0;
		             for(size_t o1=0;o1<trkHitsKey.size();o1++)
		             {
		               if(trkHitsKey.at(o1)==hitlist[ll].key())
		               {
		                 long_trk = 1;
		               }
		             }
/*
        	           int long_trk = 0;
		           for(Int_t o1=0;o1<longtrk_hits;o1++)
		           {
		             if(trkHitsKey[o1]==hitlist[ll].key())
		             {
		               long_trk = 1;
		             }
		           }
*/
                             int corr_hit = 0;
	                     for(size_t o1=0;o1<trueHitsKey.size();o1++)
	                     {
	                       if(trueHitsKey.at(o1)==hitlist[ll].key())
	                       {
	                         corr_hit = 1;
	                         cout<<"yes corr hit "<<"\n";

	                       }
                             }

			     double mdiffpeaktt0 = hitlist[ll]->PeakTime() - (T00/1000)*2;
//			     if(hitlist[ll]->WireID().TPC==2 || hitlist[ll]->WireID().TPC==6 || hitlist[ll]->WireID().TPC==10)
//		               mdiffpeaktt0 = hitlist[ll]->PeakTime() - (T00/1000)*2; //calculate the peak time in ticks
//			     if(hitlist[ll]->WireID().TPC==1 || hitlist[ll]->WireID().TPC==5 || hitlist[ll]->WireID().TPC==9)
//		               mdiffpeaktt0 = hitlist[ll]->PeakTime() + (T00/1000)*2; //calculate the peak time in ticks
		             double mhitXposcone = detprop->ConvertTicksToX(mdiffpeaktt0,hitlist[ll]->WireID().Plane,hitlist[ll]->WireID().TPC,hitlist[ll]->WireID().Cryostat); //convert ticks to usec (1tick = 0.5usec), and subtract the T0
//		           double mhitXposcone = detprop->ConvertTicksToX(mdiffpeaktt0,hitlist[ll]->WireID()); //convert ticks to usec (1tick = 0.5usec), and subtract the T0
                             double mWirestart[3], mWireend[3];
                             geometry->WireEndPoints(hitlist[ll]->WireID().Cryostat, hitlist[ll]->WireID().TPC, hitlist[ll]->WireID().Plane, hitlist[ll]->WireID().Wire, mWirestart, mWireend);
                             double mhitZposcone = mWirestart[2];
			     double mhitYposcone = 0;
                             if(fmsp.at(ll).size())
			     mhitYposcone      =   fmsp.at(ll)[0]->XYZ()[1];
			     if(!fmsp.at(ll).size()) mhitYposcone = 303.5; //average y position in the detecor FV
		  
	

/////////////////////////////////////////////////////////////////////

/*
		             //double magxzhit  = sqrt(mhitXposcone*mhitXposcone+mhitZposcone*mhitZposcone);
		             double disthits  = sqrt(pow(mhitXposcone-fnearhits_xpos[0],2)+pow(mhitZposcone-fnearhits_zpos[0],2));
		             //double anglehits = std::acos((mhitXposcone*xx + mhitZposcone*zz)/(abs(magxz)*abs(magxzhit)));
			     

			     // Constructing another line
			   if( mhitZposcone > endhitz)
			   {
			     theta2 = std::atan(slope2)*180/3.14;
			   }
			   if(mhitXposcone >= endhitx && mhitZposcone < endhitz)
			   {
			     theta2 = (std::atan(slope2)*180/3.14)+180;
			   }
//			   if(mhitXposcone < endhitx && mhitZposcone > endhitz)
//			   {
// 			     theta2 = 360+std::atan(slope2)*180/3.14;
//			   }
//			   if(mhitXposcone < endhitx && mhitZposcone < endhitz)
			   {
			     theta2 = (std::atan(slope2)*180/3.14)-180;
			   }

			     double anglehits = abs(theta-theta2);

		             if(long_trk ==0 && same_trk==0 && disthits < 20 && (anglehits) < 30)
*/

/////////////////////////////////////////////////////////////////////

/*			double maghitvec = sqrt(pow(mhitXposcone-endhitx,2)+pow(mhitZposcone-endhitz,2));
			
			double angle_theta = acos(((fnearhits_xpos[fnearhitct-1]-endhitx)*(mhitXposcone-endhitx) + (fnearhits_zpos[fnearhitct-1]-endhitz)*(mhitZposcone-endhitz))/(magnearhitvec*maghitvec));
			double angle_theta_deg = (angle_theta*180)/3.14;
*/			
	 double magnearhitvec = sqrt(pow(xposend-fsel_endhitx,2)+pow(fnearhits_zpos[int(fnearhitct/2)]-fsel_endhitz,2));
	 double maghitvec = sqrt(pow(mhitXposcone-fsel_endhitx,2)+pow(mhitZposcone-fsel_endhitz,2));
	 double angle_theta = acos(((xposend-fsel_endhitx)*(mhitXposcone-fsel_endhitx) + (fnearhits_zpos[int(fnearhitct/2)]-fsel_endhitz)*(mhitZposcone-fsel_endhitz))/(magnearhitvec*maghitvec));
	 double angle_theta_deg = (angle_theta*180)/3.14;
	 double maghitveccostheta = maghitvec* cos(angle_theta);
	 double distance = sqrt(pow(mhitXposcone-fsel_trkendx,2)+pow(mhitZposcone-fsel_trkendz,2));
			
	 if(corr_hit==1)
	 {
	   cout<<"angle_theta_deg "<<angle_theta_deg<<" maghitveccostheta "<<maghitveccostheta<<" distance "<<distance<<"\n";
         }

			if(long_trk ==0 && same_trk==0 /*&& abs(maghitveccostheta) < 80 && angle_theta_deg<30 && maghitvec > 0 && maghitvec < 80*/)
		     
		        {
/* 	                  fmhitangledegana[fhitangleana] = angle_theta_deg;   
	                  fMichelcountangledegana[fhitangleana] = corr_hit;
                          fhitangleana++;

	               ///--------------cone angle degree cut-----------------//   
		       
	               if(angle_theta_deg< _maxangledeg)
	               {
		         fanglehits_trks++;
			 
	                 _fanglehits_trks++;
	                 if(corr_hit==1)_ftrueanglehits_trks++;
*/
		               fmhits_key[fmhitcount]    = hitlist[ll].key();
		               fmhits_peakT[fmhitcount] = hitlist[ll]->PeakTime();
		               fmhits_charge[fmhitcount]= hitlist[ll]->Integral();
		               fmhits_wire[fmhitcount]  = hitlist[ll]->WireID().Wire;
		               fmhits_chno[fmhitcount]  = hitlist[ll]->Channel();
		               fmhits_TPC[fmhitcount]   = hitlist[ll]->WireID().TPC;
		               fmhits_plane[fmhitcount] = hitlist[ll]->WireID().Plane;
			       fmhits_xpos[fmhitcount]  = mhitXposcone;
			       fmhits_ypos[fmhitcount]  = mhitYposcone;
			       fmhits_zpos[fmhitcount]  = mhitZposcone;
		               fmhits_mult[fmhitcount]    = hitlist[ll]->Multiplicity();
		               fmhits_sigptime[fmhitcount]= hitlist[ll]->SigmaPeakTime();
		               fmhits_sigchrg[fmhitcount] = hitlist[ll]->SigmaIntegral();
		               fmhits_sigpamp[fmhitcount]  = hitlist[ll]->SigmaPeakAmplitude();
		               fmhits_dof[fmhitcount]     = hitlist[ll]->DegreesOfFreedom();
		               fmhits_gof[fmhitcount]     = hitlist[ll]->GoodnessOfFit();
		               fmhits_ptminusRMS[fmhitcount]     = hitlist[ll]->PeakTimeMinusRMS(5.0);
		               fmhits_ptplusRMS[fmhitcount]     = hitlist[ll]->PeakTimePlusRMS(5.0);
	                       fmhits_angledeg[fmhitcount] = angle_theta_deg;
	                       fmhits_maghitveccostheta[fmhitcount] = maghitveccostheta;
	                       fmhits_distance[fmhitcount] = distance;
			       fmhits_longtrk[fmhitcount] = long_trk;
			       fmhits_sametrk[fmhitcount] = same_trk;			       
			       fmhits_corrhit[fmhitcount] = corr_hit;
			       
//			      cout<<"ll "<<ll<<" key "<<hitlist[ll].key()<<"\n";

			       std::array<float,4> cnn_out = hitResults.getOutput( hitlist[ll].key() );
//			       double p_trk_or_sh = cnn_out[ hitResults.getIndex("track") ]+ cnn_out[ hitResults.getIndex("em") ] + cnn_out[ hitResults.getIndex("michel")];
//			       double cnn_score = cnn_out[ hitResults.getIndex("michel") ]; 
			       fmhits_cnnMichel[fmhitcount]     = cnn_out[ hitResults.getIndex("michel") ];
			       fmhits_cnnEM[fmhitcount]         = cnn_out[ hitResults.getIndex("em") ];
			       fmhits_cnnTrack[fmhitcount]      = cnn_out[ hitResults.getIndex("track") ];

//cout<<"fmhits_longtrk[fmhitcount] "<<fmhits_longtrk[fmhitcount]<<"\n";
if(fmhits_corrhit[fmhitcount]==1)cout<<"hitID "<<hitlist[ll].key()<<" long_trk "<<fmhits_longtrk[fmhitcount]<<" same_trk "<<fmhits_sametrk[fmhitcount]<<" corr_hit "<<fmhits_corrhit[fmhitcount]<<" charge "<<fmhits_charge[fmhitcount]<<" michel "<<cnn_out[ hitResults.getIndex("michel") ]<<" trk "<<cnn_out[ hitResults.getIndex("track") ]<<" em "<<cnn_out[ hitResults.getIndex("em") ]<<"\n";
			       //cout<<"hitID "<<hit_trkid<<" tracklength "<<hit_trklen<<"\n";
//		               cout<<" sel hits: fmhitcount "<<fmhitcount<<" xpos "<<fmhits_xpos[fmhitcount]<<" zpos "<<fmhits_zpos[fmhitcount]<<" disthits "<<disthits<<" slope2 "<<slope2<<" theta2 "<<theta2<<" anglehits "<<anglehits<<"\n";
//			cout<<"found cone hit key "<<hitlist[ll].key()<<" TPC "<<hitlist[ll]->WireID().TPC<<" wireID "<<hitlist[ll]->WireID().Wire<<" peakT "<<hitlist[ll]->PeakTime()<<" x, y, z "<<mhitXposcone<<" "<<mhitYposcone<<" "<<mhitZposcone<<"\n";	

//	 if(corr_hit!=0 && corr_hit!=1)cout<<"corr_hit "<<corr_hit<<"\n";

			       fmhitcount++;
		 
                            }
	                   }//if(hitlist[ll]->WireID().Plane==2)
	                 }//for(size_t ll=0; ll<hitlist.size();++ll)

                trueHitsKey.clear();
                longtrk_hitkey.clear();
/* 			 
     for(int kk=0; kk<ftrueparhitcolcount; kk++)
       {
         int near_hit = 0;
	 double magnearhitvec1 = sqrt(pow(xposend-fsel_endhitx,2)+pow(fnearhits_zpos[int(fnearhitct/2)]-fsel_endhitz,2));
	 double maghitvec1 = sqrt(pow(ftrueparhitscol_xpos[kk]-fsel_endhitx,2)+pow(ftrueparhitscol_zpos[kk]-fsel_endhitz,2));
	 double angle_theta1 = acos(((xposend-fsel_endhitx)*(ftrueparhitscol_xpos[kk]-fsel_endhitx) + (fnearhits_zpos[int(fnearhitct/2)]-fsel_endhitz)*(ftrueparhitscol_zpos[kk]-fsel_endhitz))/(magnearhitvec1*maghitvec1));
	 double angle_theta_deg1 = (angle_theta1*180)/3.14;
	 double maghitveccostheta1 = maghitvec1* cos(angle_theta1);
	 double distance1 = sqrt(pow(ftrueparhitscol_xpos[kk]-fsel_trkendx,2)+pow(ftrueparhitscol_zpos[kk]-fsel_trkendz,2));

	 for(size_t o1=0;o1<nearHitsKey.size();o1++)
	 {
	   if(nearHitsKey.at(o1)==ftrueparhitscol_key[kk])
	   {
	     near_hit = 1;
	   } 
	 }
	 if(near_hit==0)
	 {
           
	   cout<<"true hits "<<"\n";
	   cout<<"angle_theta_deg "<<angle_theta_deg1<<" maghitveccostheta "<<maghitveccostheta1<<" distance "<<distance1<<"\n";
           anglethetamisstruehits->Fill(angle_theta_deg1,maghitveccostheta1);
         }
       }//for(int kk=0; kk<ftrueparhitcolcount; kk++)
*/


      //-----------recover some michel hits from the last 10 muon hits----------------------//
	  float trunc_charge[10000] = {-999};
	  double mean_charge = 0;
	  double max_truncchrg = 0;
	  int   nnewcolhits = 0; 
	  	  
         for(int kk=1; kk<ftrkcolhits; kk++)
          {
	     if(kk<ftrkcolhits-1)
	     mean_charge = (fhits_charge[kk-1]+fhits_charge[kk]+fhits_charge[kk+1])/ 3;
	     if(kk==ftrkcolhits-1)
	     mean_charge = (fhits_charge[kk-1]+fhits_charge[kk])/ 2;
	     
	     if(fhits_charge[kk]>0.2*mean_charge && fhits_charge[kk]<2*mean_charge)
	     trunc_charge[kk] = mean_charge;
	     else
	     trunc_charge[kk] = fhits_charge[kk];

             if(kk>=ftrkcolhits-10 && kk<ftrkcolhits)
	     {
	    
	       if(trunc_charge[kk]>max_truncchrg)
	       {
	         max_truncchrg = trunc_charge[kk];
		 nnewcolhits = kk;
	       }
	     }
	    
          }//for(int kk=1; kk<ftrkcolhits; kk++)
//	  cout<<"extraHits "<<ftrkcolhits-1-nnewcolhits<<"\n";
 
 
	  std::vector <double> cone_key;
	  std::vector <double> cone_charge;
	  std::vector <int>    cone_wire;
	  std::vector <int>    cone_tpc;
	  std::vector <double> cone_ptime;
	  std::vector <double> addmu_key;
	  std::vector <double> cone_xpos;
	  std::vector <double> cone_ypos;
	  std::vector <double> cone_zpos;
	  std::vector <double> cone_angledeg;
	  std::vector <double> cone_maghitveccostheta;
	  std::vector <double> cone_distance;
	  std::vector <int>    cone_mult;
	  std::vector <double> cone_sigptime;
	  std::vector <double> cone_sigchrg;
	  std::vector <double> cone_sigpamp;
	  std::vector <double> cone_dof;
	  std::vector <double> cone_gof;
	  std::vector <double> cone_ptminusRMS;
	  std::vector <double> cone_ptplusRMS;
	  std::vector <double> cone_status;
	  std::vector <double> cone_cnnMichel;
	  std::vector <double> cone_cnnEM;
	  std::vector <double> cone_cnnTrack;
	  

	  int add_hits = 0;

          for(int kk=0; kk<ftrkcolhits; kk++)
          {
	  
	 double magnearhitvec = sqrt(pow(xposend-fsel_endhitx,2)+pow(fnearhits_zpos[int(fnearhitct/2)]-fsel_endhitz,2));
	 double maghitvec = sqrt(pow(fhits_xpos[kk]-fsel_endhitx,2)+pow(fhits_zpos[kk]-fsel_endhitz,2));
	 double angle_theta = 0;
	 if(fhits_xpos[kk]!=fsel_endhitx)
	 angle_theta = acos(((xposend-fsel_endhitx)*(fhits_xpos[kk]-fsel_endhitx) + (fnearhits_zpos[int(fnearhitct/2)]-fsel_endhitz)*(fhits_zpos[kk]-fsel_endhitz))/(magnearhitvec*maghitvec));
	 double angle_theta_deg = (angle_theta*180)/3.14;
	 double maghitveccostheta = maghitvec* cos(angle_theta);
	 double distance = sqrt(pow(fhits_xpos[kk]-fsel_trkendx,2)+pow(fhits_ypos[kk]-fsel_trkendy,2)+pow(fhits_zpos[kk]-fsel_trkendz,2));

//	  cout<<" ftrkcolhits "<<ftrkcolhits<<" nnewcolhits "<<nnewcolhits<<"\n";

            if(kk>nnewcolhits && kk<ftrkcolhits)
	    if(distance<=1.5 || (angle_theta_deg>5 && angle_theta_deg<40))
	    {
	      add_hits++;
	      cone_key.push_back(fhits_key[kk]);
	      cone_charge.push_back(fhits_charge[kk]);
	      cone_wire.push_back(fhits_wire[kk]);
	      cone_tpc.push_back(fhits_TPC[kk]);
	      cone_ptime.push_back(fhits_peakT[kk]);
	      cone_xpos.push_back(fhits_xpos[kk]);
	      cone_ypos.push_back(fhits_ypos[kk]);
	      cone_zpos.push_back(fhits_zpos[kk]);
	      cone_angledeg.push_back(angle_theta_deg);
	      cone_maghitveccostheta.push_back(maghitveccostheta);
	      cone_distance.push_back(distance);
	      cone_mult.push_back(fhits_mult[kk]);
	      cone_sigptime.push_back(fhits_sigptime[kk]);
	      cone_sigchrg.push_back(fhits_sigchrg[kk]);
	      cone_sigpamp.push_back(fhits_sigpamp[kk]);
	      cone_dof.push_back(fhits_dof[kk]);
	      cone_gof.push_back(fhits_gof[kk]);
	      cone_ptminusRMS.push_back(fhits_ptminusRMS[kk]);
	      cone_ptplusRMS.push_back(fhits_ptplusRMS[kk]);
	      cone_status.push_back(2);
	      cone_cnnMichel.push_back(fhits_cnnMichel[kk]);
	      cone_cnnEM.push_back(fhits_cnnEM[kk]);
	      cone_cnnTrack.push_back(fhits_cnnTrack[kk]);
	    }  
          }
 

       //--------analyzing michel cone hits MC-----------------//
       int inicone_hits = 0;
       for(int kk=0; kk<fmhitcount; kk++)
       {
	 double magnearhitvec = sqrt(pow(xposend-fsel_endhitx,2)+pow(fnearhits_zpos[int(fnearhitct/2)]-fsel_endhitz,2));
	 double maghitvec = sqrt(pow(fmhits_xpos[kk]-fsel_endhitx,2)+pow(fmhits_zpos[kk]-fsel_endhitz,2));
	 double angle_theta = acos(((xposend-fsel_endhitx)*(fmhits_xpos[kk]-fsel_endhitx) + (fnearhits_zpos[int(fnearhitct/2)]-fsel_endhitz)*(fmhits_zpos[kk]-fsel_endhitz))/(magnearhitvec*maghitvec));
	 double angle_theta_deg = (angle_theta*180)/3.14;
	 double maghitveccostheta = maghitvec* cos(angle_theta);
	 double distance = sqrt(pow(fmhits_xpos[kk]-fsel_trkendx,2)+pow(fmhits_ypos[kk]-fsel_trkendy,2)+pow(fmhits_zpos[kk]-fsel_trkendz,2));

          float y0 = (-(_cone_length/_cone_angle1)*angle_theta_deg)+_cone_length;
	 
         if(maghitveccostheta <= y0 && angle_theta_deg<_cone_angle1)
	 {
	  inicone_hits++;
	  cone_key.push_back(fmhits_key[kk]); 
	  cone_charge.push_back(fmhits_charge[kk]);
	  cone_wire.push_back(fmhits_wire[kk]);
	  cone_tpc.push_back(fmhits_TPC[kk]);
	  cone_ptime.push_back(fmhits_peakT[kk]);
	  cone_xpos.push_back(fmhits_xpos[kk]);
	  cone_ypos.push_back(fmhits_ypos[kk]);
	  cone_zpos.push_back(fmhits_zpos[kk]);
	  cone_angledeg.push_back(angle_theta_deg);
	  cone_maghitveccostheta.push_back(maghitveccostheta);
	  cone_distance.push_back(distance);
	  cone_mult.push_back(fmhits_mult[kk]);
	  cone_sigptime.push_back(fmhits_sigptime[kk]);
	  cone_sigchrg.push_back(fmhits_sigchrg[kk]);
	  cone_sigpamp.push_back(fmhits_sigpamp[kk]);
	  cone_dof.push_back(fmhits_dof[kk]);
	  cone_gof.push_back(fmhits_gof[kk]);
	  cone_ptminusRMS.push_back(fmhits_ptminusRMS[kk]);
	  cone_ptplusRMS.push_back(fmhits_ptplusRMS[kk]);
	  cone_status.push_back(1);
	  cone_cnnMichel.push_back(fmhits_cnnMichel[kk]);
	  cone_cnnEM.push_back(fmhits_cnnEM[kk]);
	  cone_cnnTrack.push_back(fmhits_cnnTrack[kk]);
	 } 
        }
        
	//-----Sort hits based on increasing peak times--------//
			       
/*        for(size_t o1=0;o1<cone_key.size();o1++)
	{
	  cout<<"before ptime "<<cone_ptime.at(o1)<<"\n";
	}  
*/
        float o3=0,o4=0,o5=0,o6=0,o7=0,o8=0,o9=0,o10=0,o11=0,o12=0,o13=0,o14=0,o15=0,o16=0,o17=0,o18=0,o19=0,o20=0,o21=0,o22=0,o23=0,o24=0,o25=0;
	
        for(size_t o1=0;o1<cone_key.size();o1++)
	{
	  for(size_t o2=o1+1; o2 < cone_key.size(); o2++)
	  {
	    if(cone_ptime.at(o1) > cone_ptime.at(o2))
	    {
	      o3 = cone_key.at(o1);
	      cone_key.at(o1) = cone_key.at(o2);
	      cone_key.at(o2) = o3;

	      o4 = cone_charge.at(o1);
	      cone_charge.at(o1) = cone_charge.at(o2);
	      cone_charge.at(o2) = o4;

	      o5 = cone_wire.at(o1);
	      cone_wire.at(o1) = cone_wire.at(o2);
	      cone_wire.at(o2) = o5;

	      o6 = cone_tpc.at(o1);
	      cone_tpc.at(o1) = cone_tpc.at(o2);
	      cone_tpc.at(o2) = o6;
	      
	      o7 = cone_ptime.at(o1);
	      cone_ptime.at(o1) = cone_ptime.at(o2);
	      cone_ptime.at(o2) = o7;

	      o8 = cone_xpos.at(o1);
	      cone_xpos.at(o1) = cone_xpos.at(o2);
	      cone_xpos.at(o2) = o8;

	      o9 = cone_ypos.at(o1);
	      cone_ypos.at(o1) = cone_ypos.at(o2);
	      cone_ypos.at(o2) = o9;

	      o10 = cone_zpos.at(o1);
	      cone_zpos.at(o1) = cone_zpos.at(o2);
	      cone_zpos.at(o2) = o10;

	      o11 = cone_angledeg.at(o1);
	      cone_angledeg.at(o1) = cone_angledeg.at(o2);
	      cone_angledeg.at(o2) = o11;

	      o12 = cone_maghitveccostheta.at(o1);
	      cone_maghitveccostheta.at(o1) = cone_maghitveccostheta.at(o2);
	      cone_maghitveccostheta.at(o2) = o12;

	      o13 = cone_distance.at(o1);
	      cone_distance.at(o1) = cone_distance.at(o2);
	      cone_distance.at(o2) = o13;

	      o14 = cone_mult.at(o1);
	      cone_mult.at(o1) = cone_mult.at(o2);
	      cone_mult.at(o2) = o14;

	      o15 = cone_sigptime.at(o1);
	      cone_sigptime.at(o1) = cone_sigptime.at(o2);
	      cone_sigptime.at(o2) = o15;

	      o16 = cone_sigchrg.at(o1);
	      cone_sigchrg.at(o1) = cone_sigchrg.at(o2);
	      cone_sigchrg.at(o2) = o16;

	      o17 = cone_sigpamp.at(o1);
	      cone_sigpamp.at(o1) = cone_sigpamp.at(o2);
	      cone_sigpamp.at(o2) = o17;

	      o18 = cone_dof.at(o1);
	      cone_dof.at(o1) = cone_dof.at(o2);
	      cone_dof.at(o2) = o18;

	      o19 = cone_gof.at(o1);
	      cone_gof.at(o1) = cone_gof.at(o2);
	      cone_gof.at(o2) = o19;

	      o20 = cone_ptminusRMS.at(o1);
	      cone_ptminusRMS.at(o1) = cone_ptminusRMS.at(o2);
	      cone_ptminusRMS.at(o2) = o20;

	      o21 = cone_ptplusRMS.at(o1);
	      cone_ptplusRMS.at(o1) = cone_ptplusRMS.at(o2);
	      cone_ptplusRMS.at(o2) = o21;

	      o22 = cone_status.at(o1);
	      cone_status.at(o1) = cone_status.at(o2);
	      cone_status.at(o2) = o22;

	      o23 = cone_cnnMichel.at(o1);
	      cone_cnnMichel.at(o1) = cone_cnnMichel.at(o2);
	      cone_cnnMichel.at(o2) = o23;

	      o24 = cone_cnnEM.at(o1);
	      cone_cnnEM.at(o1) = cone_cnnEM.at(o2);
	      cone_cnnEM.at(o2) = o24;

	      o25 = cone_cnnTrack.at(o1);
	      cone_cnnTrack.at(o1) = cone_cnnTrack.at(o2);
	      cone_cnnTrack.at(o2) = o25;
	    }
	  }
	}      
	  
/*        for(size_t o1=0;o1<cone_key.size();o1++)
	{
	  cout<<"after ptime "<<cone_ptime.at(o1)<<"\n";
	}  
*/	  
	//-----Stack up all cone + additional mu hits

        fmichel_conesize = cone_key.size();
	
        std::map<int,std::vector<double>> hit_t1;
        std::map<int,std::vector<double>> hit_t2;
        std::map<int,std::vector<double>> hit_y;
        std::map<int,std::vector<double>> hit_x;
        std::map<int,std::vector<double>> hit_key;
        std::map<int,std::vector<double>> hit_wire;
        std::map<int,std::vector<double>> hit_tpc;
        std::map<int,std::vector<double>> hit_chargehit;
        std::map<int,std::vector<double>> hit_ptime;
        std::map<int,std::vector<double>> hit_angledeg;
        std::map<int,std::vector<double>> hit_maghitveccostheta;
        std::map<int,std::vector<double>> hit_distance;
        std::map<int,std::vector<double>> hit_mult;
        std::map<int,std::vector<double>> hit_sigptime;
        std::map<int,std::vector<double>> hit_sigchrg;
        std::map<int,std::vector<double>> hit_sigpamp;
        std::map<int,std::vector<double>> hit_dof;
        std::map<int,std::vector<double>> hit_gof;
        std::map<int,std::vector<double>> hit_ptminusRMS;
        std::map<int,std::vector<double>> hit_ptplusRMS;
        std::map<int,std::vector<double>> hit_status;
        std::map<int,std::vector<double>> hit_cnnMichel;
        std::map<int,std::vector<double>> hit_cnnEM;
        std::map<int,std::vector<double>> hit_cnnTrack;
	
	fmichel_wcount = 0;
	
	for(size_t o1=0;o1<cone_key.size();o1++)
	{
	  std::vector<art::Ptr<recob::Wire>> wirescone = wFromHits.at(cone_key.at(o1));
          hit_t1[wirescone[0]->Channel()].push_back(cone_ptminusRMS.at(o1));
          hit_t2[wirescone[0]->Channel()].push_back(cone_ptplusRMS.at(o1));
          hit_key[wirescone[0]->Channel()].push_back(cone_key.at(o1));
          hit_wire[wirescone[0]->Channel()].push_back(cone_wire.at(o1));
          hit_tpc[wirescone[0]->Channel()].push_back(cone_tpc.at(o1));
          hit_chargehit[wirescone[0]->Channel()].push_back(cone_charge.at(o1));
          hit_ptime[wirescone[0]->Channel()].push_back(cone_ptime.at(o1));
          hit_angledeg[wirescone[0]->Channel()].push_back(cone_angledeg.at(o1));
          hit_maghitveccostheta[wirescone[0]->Channel()].push_back(cone_maghitveccostheta.at(o1));
          hit_distance[wirescone[0]->Channel()].push_back(cone_distance.at(o1));
          hit_mult[wirescone[0]->Channel()].push_back(cone_mult.at(o1));
          hit_sigptime[wirescone[0]->Channel()].push_back(cone_sigptime.at(o1));
          hit_sigchrg[wirescone[0]->Channel()].push_back(cone_sigchrg.at(o1));
          hit_sigpamp[wirescone[0]->Channel()].push_back(cone_sigpamp.at(o1));
          hit_dof[wirescone[0]->Channel()].push_back(cone_dof.at(o1));
          hit_gof[wirescone[0]->Channel()].push_back(cone_gof.at(o1));
          hit_ptminusRMS[wirescone[0]->Channel()].push_back(cone_ptminusRMS.at(o1));
          hit_ptplusRMS[wirescone[0]->Channel()].push_back(cone_ptplusRMS.at(o1));
          hit_status[wirescone[0]->Channel()].push_back(cone_status.at(o1));
          hit_cnnMichel[wirescone[0]->Channel()].push_back(cone_cnnMichel.at(o1));
          hit_cnnEM[wirescone[0]->Channel()].push_back(cone_cnnEM.at(o1));
          hit_cnnTrack[wirescone[0]->Channel()].push_back(cone_cnnTrack.at(o1));
       //std::cout<<wires[0]->Channel()<<" "<<sh_hits[j]->PeakTime()<<" "<<sh_hits[j]->PeakTimePlusRMS(5.0)<<" "<<sh_hits[j]->PeakTimeMinusRMS(5.0)<<std::endl;
//          double diffconepeaktt0     = cone_ptime.at(o1) - (T00/1000)*2;
//         hit_x[wirescone[0]->Channel()].push_back(detprop->ConvertTicksToX(diffconepeaktt0,2,cone_tpc.at(o1),0));
          hit_x[wirescone[0]->Channel()].push_back(cone_xpos.at(o1));	 	  
//          std::vector<art::Ptr<recob::SpacePoint>> sp = fmsp.at(cone_key.at(o1)); 
//          if(!sp.empty())
//	  {
//            hit_y[wirescone[0]->Channel()].push_back(sp[0]->XYZ()[1]);
//          }
//          else 
	  hit_y[wirescone[0]->Channel()].push_back(cone_ypos.at(o1)); 
        }   
 
	
        //one wire can have various hits so lets remove duplicate wires
        auto const& wires = evt.getValidHandle<std::vector<recob::Wire> >(fWireModuleLabel);
        //auto const& wires = evt.getValidHandle<std::vector<recob::Wire> >("wclsdatasp:gauss"); //new reco 
        auto w1 = hit_t1.begin();
        auto w2 = hit_t2.begin();
        auto x  = hit_x.begin();
        auto y  = hit_y.begin();
	auto key = hit_key.begin(); 
	auto wi = hit_wire.begin();
	auto chargehit = hit_chargehit.begin();
	auto tpc = hit_tpc.begin();
	auto ptime = hit_ptime.begin();
	auto angledeg = hit_angledeg.begin();
	auto maghitveccostheta = hit_maghitveccostheta.begin();
	auto distance = hit_distance.begin();
	auto mult = hit_mult.begin();
	auto sigptime = hit_sigptime.begin();
	auto sigchrg = hit_sigchrg.begin();
	auto sigpamp = hit_sigpamp.begin();
	auto dof = hit_dof.begin();
	auto gof = hit_gof.begin();
	auto ptminusRMS = hit_ptminusRMS.begin();
	auto ptplusRMS = hit_ptplusRMS.begin();
	auto status = hit_status.begin();
	auto cnnMichel = hit_cnnMichel.begin();
	auto cnnEM = hit_cnnEM.begin();
	auto cnnTrack = hit_cnnTrack.begin();
	
        while( w1 != hit_t1.end())
	{
          int it_w = w1->first;
          int n_hits = w1->second.size();
           std::cout<<"it_w "<<it_w<<" n_hits "<<n_hits<<" w1->second[0] "<<w1->second[0]<<"  w2->second[0] "<<w2->second[0]<<" w1->second[n_hits-1] "<<w1->second[n_hits-1]<<" w2->second[n_hits-1] "<<w2->second[n_hits-1]<<std::endl;
          double t1 =  w1->second[0]; //first hit
          double t2 =  w2->second[n_hits-1]; //last hit
//          double t2 =  w2->second[n_hits-1]; //last hit
          double x_w, y_w, key_w, wi_w, chargehit_w, tpc_w, ptime_w, angledeg_w, maghitveccostheta_w, distance_w, mult_w, sigptime_w, 
	  sigchrg_w, sigpamp_w, dof_w, gof_w, ptminusRMS_w, ptplusRMS_w, status_w, cnnMichel_w, cnnEM_w, cnnTrack_w;
/*          if( x->second.size() < 1 )
	  {
            x_w = (x->second[0]-x->second[n_hits-1])/2.0;
            y_w = (y->second[0]-y->second[n_hits-1])/2.0;
          }
          else 
	  {
            x_w = x->second[0];
            y_w = y->second[0];
          }
*/
           x_w = x->second[0];
           y_w = y->second[0];
	  key_w = key->second[0];
          wi_w = wi->second[0];
	  tpc_w = tpc->second[0];
	  chargehit_w = chargehit->second[0];
	  ptime_w = ptime->second[0];
	  angledeg_w = angledeg->second[0];
	  maghitveccostheta_w = maghitveccostheta->second[0];
	  distance_w = distance->second[0];
	  mult_w = mult->second[0];
	  sigptime_w = sigptime->second[0];
	  sigchrg_w = sigchrg->second[0];
	  sigpamp_w = sigpamp->second[0];
	  dof_w = dof->second[0];
	  gof_w = gof->second[0];
	  ptminusRMS_w = ptminusRMS->second[0];
	  ptplusRMS_w = ptplusRMS->second[0];
	  status_w = status->second[0];
	  cnnMichel_w = cnnMichel->second[0];
	  cnnEM_w = cnnEM->second[0];
	  cnnTrack_w = cnnTrack->second[0];
	  
          for(auto & wire : * wires)
	  {
            int channel_no = wire.Channel();
            int plane = fGeometry->View(wire.Channel()); 
            if( plane != 2 ) continue;
            std::vector< geo::WireID > wireID= fGeometry->ChannelToWire(channel_no);
            const geo::WireGeo* pwire = fGeometry->WirePtr(wireID.at(0)); //for collection plane there is one wire per channel
            TVector3 xyzWire = pwire->GetCenter<TVector3>();
            if( it_w == channel_no )
	    {  
              double charge =0.0;
              for(size_t i = 0; i < wire.Signal().size(); ++i)
	      {
                if( i > t1 && i < t2 ) charge += wire.Signal()[i];
              }   
              fmichel_zpos[fmichel_wcount] = (xyzWire.Z()); 
              fmichel_ypos[fmichel_wcount] = (y_w); 
              fmichel_xpos[fmichel_wcount] = (x_w); 
              fmichel_chrg[fmichel_wcount] = (charge);
              fmichel_chno[fmichel_wcount] = (channel_no);
	      fmichel_key[fmichel_wcount] = (key_w);
	      fmichel_wire[fmichel_wcount] = (wi_w);
	      fmichel_chargehit[fmichel_wcount] = (chargehit_w);
	      fmichel_tpc[fmichel_wcount]  = (tpc_w);
	      fmichel_ptime[fmichel_wcount]= (ptime_w);
	      fmichel_angledeg[fmichel_wcount]= (angledeg_w);
	      fmichel_maghitveccostheta[fmichel_wcount]= (maghitveccostheta_w);
	      fmichel_distance[fmichel_wcount]= (distance_w);
	      fmichel_mult[fmichel_wcount]= (mult_w);
	      fmichel_sigptime[fmichel_wcount]= (sigptime_w);
	      fmichel_sigchrg[fmichel_wcount]= (sigchrg_w);
	      fmichel_sigpamp[fmichel_wcount]= (sigpamp_w);
	      fmichel_dof[fmichel_wcount]= (dof_w);
	      fmichel_gof[fmichel_wcount]= (gof_w);
	      fmichel_ptminusRMS[fmichel_wcount]= (ptminusRMS_w);
	      fmichel_ptplusRMS[fmichel_wcount]= (ptplusRMS_w);
	      fmichel_status[fmichel_wcount]= (status_w);
	      fmichel_cnnMichel[fmichel_wcount]= (cnnMichel_w);
	      fmichel_cnnEM[fmichel_wcount]= (cnnEM_w);
	      fmichel_cnnTrack[fmichel_wcount]= (cnnTrack_w);
              fmichel_wcount++;
	      
//cout<<"xyzWire.Z() "<<xyzWire.Z()<<" y "<<y_w<<" x "<<x_w<<" charge "<<charge<<" t1 "<<t1<<" t2 "<<t2<<"\n"<<"\n";

	      cout<<"wi_w "<<wi_w<<" tpc_w "<<tpc_w<<" ptime_w "<<ptime_w<<" mult_w "<<mult_w<<" michel_z "<<xyzWire.Z()<<" y "<<y_w<<" x "<<x_w<<" charge "<<charge<<" hit_chrge "<<chargehit_w<<" t1 "<<t1<<" t2 "<<t2<<" status "<<status_w<<"\n"<<"\n";
              break;
            }   
          }// all recob::wire
          w1 ++; w2 ++;
          x ++; y ++;
	  key++; wi++; chargehit++; tpc++; ptime++; angledeg++; maghitveccostheta++; distance++; mult++; sigptime++; sigchrg++; sigpamp++;
	  dof++; gof++; ptminusRMS++; ptplusRMS++, status++;
        }//wire from cone hits 




	  cone_key.clear(); 
	  cone_charge.clear();
	  cone_wire.clear();
	  cone_tpc.clear();
	  cone_ptime.clear();
	  cone_xpos.clear();
	  cone_ypos.clear();
	  cone_zpos.clear();
	  cone_angledeg.clear();
	  cone_maghitveccostheta.clear();
	  cone_distance.clear();
	  cone_mult.clear();
	  cone_sigptime.clear();
	  cone_sigchrg.clear();
	  cone_sigpamp.clear();
	  cone_dof.clear();
	  cone_gof.clear();
	  cone_ptminusRMS.clear();
	  cone_ptplusRMS.clear();
	  cone_status.clear();
	  cone_cnnMichel.clear();
	  cone_cnnEM.clear();
	  cone_cnnTrack.clear();

	                  //-----------getting the MC true information of the selected cone michel hits-----------------//
	       
 		       if(isMC) 
		       { 
                         std::map<int,double> mconetrkide;
	                  int mconetrackid=-1;
                          for(int qq =0; qq < fmhitcount; qq++)
		          {
		            //cout<<"hitKey[qq] "<<hitKey[qq]<<"\n";
		            art::Ptr<recob::Hit> mconehit=hitlist[fmhits_key[qq]];    
//			    cout<<"cone hit plane "<<mconehit->WireID().Plane<<"\n";
//		            cout<<"michel cone hit key "<<mconehit.key()<<" peakT "<<mconehit->PeakTime()<<" WireID "<<mconehit->WireID().Wire<<"\n";
	                    std::vector<sim::TrackIDE> meveIDs = bt_serv->HitToTrackIDEs(mconehit);

	                    for(size_t e=0;e<meveIDs.size(); ++e)
	                    { 
		              if(abs(meveIDs[e].trackID)!=fmcsel_trkid) //the hits should not belong to the same mother track (muon)
		              {
	                        mconetrkide[meveIDs[e].trackID] += meveIDs[e].energy;
//		                cout<<"meveIDs[e].trackID "<<meveIDs[e].trackID<<" mconetrkide[meveIDs[e].trackID] "<<mconetrkide[meveIDs[e].trackID]<<"\n";
		              }
	                    }
                          }
                          double  mconemaxe = -1;
                          double mconetote = 0;
                          for(std::map<int,double>::iterator ii = mconetrkide.begin(); ii!=mconetrkide.end(); ++ii)
                          {
	                    //	cout<<" trkid "<<ii->first<<"  energy deposited = "<<ii->second<<endl;
	                    mconetote += ii->second;
	                    if((ii->second)>mconemaxe)
	                    {
	                      mconemaxe = ii->second;
	                      mconetrackid = abs(ii->first);
	                    }
                          }
/*                        float mtotal_energy=0.0;
                          for(int qq =0; qq < hitcount; qq++)
		          {
	                    art::Ptr<recob::Hit> hit=hitlist[hitKey.at(qq)];   
	                    if (hit->WireID().Plane!=2) continue;
	                    std::vector<sim::TrackIDE> eveIDs = bt_serv->HitToTrackIDEs(hit);
	                    for(size_t e=0;e<eveIDs.size(); ++e)
		           {
	                     if (eveIDs[e].trackID == mtrackid) mtotal_energy +=  eveIDs[e].energy;
	                   }
                         }
*/ 
                         const simb::MCParticle *mconeparticleP = pi_serv->TrackIdToParticle_P(mconetrackid);
                         if(!mconeparticleP) continue;

                         const art::Ptr<simb::MCTruth> mconemcP=pi_serv->TrackIdToMCTruth_P(mconetrackid);
                         if(!mconemcP) continue;

	                 TLorentzVector mconemcstart, mconemcend, mconemcstartdrifted, mconemcenddrifted;
	                 unsigned int mconepstarti, mconependi, mconepstartdriftedi, mconependdriftedi; //mcparticle indices for starts and ends in tpc or drifted volumes
	                 double mconeplen = length(*mconeparticleP, mconemcstart, mconemcend, mconepstarti, mconependi);
	                 double mconeplendrifted = driftedLength(*mconeparticleP, mconemcstartdrifted, mconemcenddrifted, mconepstartdriftedi, mconependdriftedi);
	                 // bool isActive = plen != 0;
	                 bool mconeisDrifted = mconeplendrifted!= 0;

                         cout<<"Michel cone hits: bt:trackid "<<mconetrackid<<"\n";
                         fmcconehits_trkid      = mconeparticleP->TrackId();
                         fmcconehits_vx         = mconeparticleP->Vx();
                         fmcconehits_vy         = mconeparticleP->Vy();
                         fmcconehits_vz         = mconeparticleP->Vz();
                         fmcconehits_t          = mconeparticleP->T();
                         fmcconehits_endx       = mconeparticleP->EndX();
                         fmcconehits_endy       = mconeparticleP->EndY();
                         fmcconehits_endz       = mconeparticleP->EndZ();
                         fmcconehits_endt       = mconeparticleP->EndT();
                         fmcconehits_px         = mconeparticleP->Px();
                         fmcconehits_py         = mconeparticleP->Py();
                         fmcconehits_pz         = mconeparticleP->Pz();
                         fmcconehits_momentum   = mconeparticleP->P();
                         fmcconehits_energy     = mconeparticleP->E();
                         fmcconehits_endpx      = mconeparticleP->EndPx();
                         fmcconehits_endpy      = mconeparticleP->EndPy();
                         fmcconehits_endpz      = mconeparticleP->EndPz();
                         fmcconehits_endenergy  = mconeparticleP->EndE();
	                 fmcconehits_pathlen    = mconeplen;
	                 if (mconeisDrifted)
	                 {
                           fmcconehits_vxdrifted         = mconeparticleP->Vx();
                           fmcconehits_vydrifted         = mconeparticleP->Vy();
                           fmcconehits_vzdrifted         = mconeparticleP->Vz();
                           fmcconehits_tdrifted          = mconeparticleP->T();
                           fmcconehits_endxdrifted       = mconeparticleP->EndX();
                           fmcconehits_endydrifted       = mconeparticleP->EndY();
                           fmcconehits_endzdrifted       = mconeparticleP->EndZ();
                           fmcconehits_endtdrifted       = mconeparticleP->EndT();
                           fmcconehits_pxdrifted         = mconeparticleP->Px();
                           fmcconehits_pydrifted         = mconeparticleP->Py();
                           fmcconehits_pzdrifted         = mconeparticleP->Pz();
                           fmcconehits_momentumdrifted   = mconeparticleP->P();
                           fmcconehits_energydrifted     = mconeparticleP->E();
                           fmcconehits_endpxdrifted      = mconeparticleP->EndPx();
                           fmcconehits_endpydrifted      = mconeparticleP->EndPy();
                           fmcconehits_endpzdrifted      = mconeparticleP->EndPz();
                           fmcconehits_endenergydrifted  = mconeparticleP->EndE();
	                   fmcconehits_pathlendrifted    = mconeplendrifted;
	                 }
	                 fmcconehits_endprocess = int(mconeparticleP->Trajectory().ProcessToKey(mconeparticleP->EndProcess()));
	                 fmcconehits_theta      = mconeparticleP->Momentum().Theta();
	                 fmcconehits_phi        = mconeparticleP->Momentum().Phi();
                         fmcconehits_pdg        = mconeparticleP->PdgCode();
                         fmcconehits_status_code= mconeparticleP->StatusCode();
                         fmcconehits_mass       = mconeparticleP->Mass();
                         fmcconehits_ND         = mconeparticleP->NumberDaughters();
                         fmcconehits_mother     = mconeparticleP->Mother();
	                 fmcconehits_origin     = mconemcP->Origin();
	                 fmcconehits_process    = int(mconeparticleP->Trajectory().ProcessToKey(mconeparticleP->Process()));
                         fmcconehits_rescatter  = mconeparticleP->Rescatter();
	
                         std::cout<<"MC true Michel cone hits: Pdg code "<<fmcconehits_pdg<<" track ID "<<fmcconehits_trkid<<" no of daughters "<<fmcconehits_ND<<" and origin_type "<<fmcconehits_origin<<std::endl;
 	     
                         /*********************MC shower true information*******************************/
/*                            
                           // Get true MCParticle associated with recob::Track
                           const simb::MCParticle *shparticleP = truthUtil.GetMCParticleFromRecoShower(shower,evt,fShowerModuleLabel);
                           if(!shparticleP) continue;
                           const art::Ptr<simb::MCTruth> shmcP=pi_serv->TrackIdToMCTruth_P(shparticleP->TrackId());
                           if(!shmcP) continue;
*/
			 std::vector<art::Ptr<recob::Hit>> shwrHits= fshwrHit.at(fselshwr_key); //storing hits for ith track
                         std::map<int,double> shwride;
	                 int shwrid=-1;

                         for(size_t h=0; h<shwrHits.size();h++)
                         {
	                   art::Ptr<recob::Hit> shhit=shwrHits[h];
	                   std::vector<sim::TrackIDE> sheveIDs = bt_serv->HitToTrackIDEs(shhit);
	                   for(size_t e=0;e<sheveIDs.size(); ++e)
	                   {
		              if(abs(sheveIDs[e].trackID)!=fmcsel_trkid) //the hits should not belong to the same mother particle (muon)
		              {
	                        shwride[sheveIDs[e].trackID] += sheveIDs[e].energy;
			      }
	                   }
                         }
                         double  shmaxe = -1;
                         double shtote = 0;
                         for(std::map<int,double>::iterator ii = shwride.begin(); ii!=shwride.end(); ++ii)
                         {
	                   //	cout<<" trkid "<<ii->first<<"  energy deposited = "<<ii->second<<endl;
	                   shtote += ii->second;
	                   if((ii->second)>shmaxe)
	                   {
	                     shmaxe = ii->second;
	                     shwrid = ii->first;
	                   }
                         }
  
                         const simb::MCParticle *shparticleP = pi_serv->TrackIdToParticle_P(shwrid);
                         if(!shparticleP) continue;
                         const art::Ptr<simb::MCTruth> shmcP=pi_serv->TrackIdToMCTruth_P(shwrid);
                         if(!shmcP) continue;
	                 //const simb::Origin_t parto(mctruthproto->GetParticle(iPartp));

	                 TLorentzVector shmcstart, shmcend, shmcstartdrifted, shmcenddrifted;
	                 unsigned int shpstarti, shpendi, shpstartdriftedi, shpenddriftedi; //mcparticle indices for starts and ends in tpc or drifted volumes
	                 double shplen = length(*shparticleP, shmcstart, shmcend, shpstarti, shpendi);
	                 double shplendrifted = driftedLength(*shparticleP, shmcstartdrifted, shmcenddrifted, shpstartdriftedi, shpenddriftedi);
	                 bool isshDrifted = shplendrifted!= 0;

                         fmcshwr_trkid     = shparticleP->TrackId();
                         fmcshwr_vx        = shparticleP->Vx();
                         fmcshwr_vy        = shparticleP->Vy();
                         fmcshwr_vz        = shparticleP->Vz();
                         fmcshwr_t         = shparticleP->T();
                         fmcshwr_endx      = shparticleP->EndX();
                         fmcshwr_endy      = shparticleP->EndY();
                         fmcshwr_endz      = shparticleP->EndZ();
                         fmcshwr_endt      = shparticleP->EndT();
                         fmcshwr_px        = shparticleP->Px();
                         fmcshwr_py        = shparticleP->Py();
                         fmcshwr_pz        = shparticleP->Pz();
                         fmcshwr_momentum  = shparticleP->P();
                         fmcshwr_energy    = shparticleP->E();
                         fmcshwr_endpx     = shparticleP->EndPx();
                         fmcshwr_endpy     = shparticleP->EndPy();
                         fmcshwr_endpz     = shparticleP->EndPz();
                         fmcshwr_endenergy = shparticleP->EndE();
	                 fmcshwr_pathlen   = shplen;
	                 if (isshDrifted)
	                 {
                           fmcshwr_vxdrifted        = shparticleP->Vx();
                           fmcshwr_vydrifted        = shparticleP->Vy();
                           fmcshwr_vzdrifted        = shparticleP->Vz();
                           fmcshwr_tdrifted         = shparticleP->T();
                           fmcshwr_endxdrifted      = shparticleP->EndX();
                           fmcshwr_endydrifted      = shparticleP->EndY();
                           fmcshwr_endzdrifted      = shparticleP->EndZ();
                           fmcshwr_endtdrifted      = shparticleP->EndT();
                           fmcshwr_pxdrifted        = shparticleP->Px();
                           fmcshwr_pydrifted        = shparticleP->Py();
                           fmcshwr_pzdrifted        = shparticleP->Pz();
                           fmcshwr_momentumdrifted  = shparticleP->P();
                           fmcshwr_energydrifted    = shparticleP->E();
                           fmcshwr_endpxdrifted     = shparticleP->EndPx();
                           fmcshwr_endpydrifted     = shparticleP->EndPy();
                           fmcshwr_endpzdrifted     = shparticleP->EndPz();
                           fmcshwr_endenergydrifted = shparticleP->EndE();
	                   fmcshwr_pathlendrifted   = shplendrifted;
	                 }
	                 fmcshwr_endprocess = int(shparticleP->Trajectory().ProcessToKey(shparticleP->EndProcess()));
	                 fmcshwr_theta      = shparticleP->Momentum().Theta();
	                 fmcshwr_phi        = shparticleP->Momentum().Phi();
                         fmcshwr_pdg        = shparticleP->PdgCode();
                         fmcshwr_status_code= shparticleP->StatusCode();
                         fmcshwr_mass       = shparticleP->Mass();
                         fmcshwr_ND         = shparticleP->NumberDaughters();
                         fmcshwr_mother     = shparticleP->Mother();
	                 fmcshwr_origin     = shmcP->Origin();
	                 fmcshwr_process    = int(shparticleP->Trajectory().ProcessToKey(shparticleP->Process()));
                         fmcshwr_rescatter  = shparticleP->Rescatter();
	
                         std::cout<<"MC shower: Pdg code "<<fmcshwr_pdg<<" track ID "<<fmcshwr_trkid<<" no of daughters "<<fmcshwr_ND<<" and origin_type "<<fmcshwr_origin<<" energy "<<fmcshwr_energy<<std::endl;
		       }//if(isMC)  


                        //------------------ storing selected flash information	----------------------------//	
//                         for (size_t lm = 0; lm < ophitlistInt.size(); ++lm)
//                         {
//                           art::Ptr<recob::OpHit> pophitInt(ophitListHandleInt, lm);
//                           const recob::OpHit& ophit = *pophitInt;
//		           cout<<"ophit ch: "<<ophit.OpChannel()<<" peaktime "<<ophit.PeakTime()<<" frame "<<ophit.Frame()<<" width "<<ophit.Width()<<" area "<<ophit.Area()<<" amp "<<ophit.Amplitude()<<" pe "<<ophit.PE()<<"\n";
//			 } 

                 	double TPC_trigger_offset = 0.0;

		        // Get trigger to TPC Offset
		        auto const* detclock = lar::providerFrom<detinfo::DetectorClocksService>();
		        TPC_trigger_offset = detclock->TriggerOffsetTPC();
		        std::cout << "TPC time offset from trigger: " << TPC_trigger_offset << "\n";

			float flsh_trk_time = 99999;
			totflash = 0;
    
                        for (size_t nn = 0; nn < flashlistInt.size() && nn < kMaxFlashes ; ++nn)
                        {
                          art::Ptr<recob::OpFlash> pflashInt(flashListHandleInt, nn);
                          const recob::OpFlash& flashInt = *pflashInt;

	                 // Find the closest in time external flash
                         double minflashtime = 9999999999999, minflashtime1 = 999; 
			 fext_trigger_time = 999;
			 
                         for (size_t nm = 0; nm < flashlistExt.size() && nm < kMaxFlashes ; ++nm)
                         {
			   //loop over extrnal flashes
                           art::Ptr<recob::OpFlash> pmflashExt(flashListHandleExt,nm);
                           const recob::OpFlash& mflashExt = *pmflashExt;
			   
//			   cout<<"mflashExt.Time() "<<mflashExt.Time()<<"\n";
		           minflashtime1 = abs(flashInt.Time() - mflashExt.Time());
		           if(minflashtime1 < minflashtime )
		           {
		             minflashtime = minflashtime1;
		             fext_trigger_time = mflashExt.Time();
		           }
		         }
                                                  			 			   
/*
                         // Get OpFlash and associated ophits
                         std::vector< art::Ptr<recob::OpHit> > matchedopHits = fophitflsh.at(nn);

			 double npe = 0;
	                 for(size_t lm=0; lm<matchedopHits.size();lm++)
                         {
	                   art::Ptr<recob::OpHit> ophit = matchedopHits[lm];
			   npe += ophit->PE();
		           cout<<"ophit ch: "<<ophit->OpChannel()<<" peaktime "<<ophit->PeakTime()<<" frame "<<ophit->Frame()<<" width "<<ophit->Width()<<" area "<<ophit->Area()<<" amp "<<ophit->Amplitude()<<" pe "<<ophit->PE()<<"\n";
                        
			 }
			 cout<<"tot pe "<<npe<<"\n";
*/			 


//	                    float fall_flash_time_diffs = (flashInt.Time()-fext_trigger_time) - (fsel_trkrecotime/*+TPC_trigger_offset*/); // Time diff relative to trigger
float fall_flash_time_diffs = (flashInt.Time()-fext_trigger_time) - (fsel_trkrecotime); // Time diff relative to trigger			  
//			    flash_distana[totflash] = min_dist;
//			    flash_peana[totflash] = flashInt.TotalPE();
			  
//			    double dist_trk_flash = sqrt(pow((flashInt.YCenter()-fsel_endhity),2)+pow((flashInt.ZCenter()-fsel_endhitz),2));
//			    cout<<"Int-ext_time "<<(flashInt.Time()-fext_trigger_time)<<" fsel_trkrecotime "<<fsel_trkrecotime<<" dist_trk_flash "<<dist_trk_flash<<" sel flash PE "<<flashInt.TotalPE()<<" time_diff "<<abs(fall_flash_time_diffs)<<" y-center "<<flashInt.YCenter()<<" z-center "<<flashInt.ZCenter()<<" y-width "<<flashInt.YWidth()<<" z-width "<<flashInt.ZWidth()<<"\n";

			    if(abs(fall_flash_time_diffs)< flsh_trk_time && flashInt.TotalPE()>20)
			    {
 	  	            flsh_trk_time = abs(fall_flash_time_diffs);
			    
			    
			    fflashsel_track_dist_int  = flsh_trk_time;
                            fflashsel_time_int        = flashInt.Time(); // Time relative to trigger
                            fflashsel_pe_int          = flashInt.TotalPE();
                            fflashsel_ycenter_int     = flashInt.YCenter(); // Geometric center in y
                            fflashsel_zcenter_int     = flashInt.ZCenter(); // Geometric center in z
                            fflashsel_ywidth_int      = flashInt.YWidth(); // Geometric width in y
                            fflashsel_zwidth_int      = flashInt.ZWidth(); // Geometric width in z
                            fflashsel_timewidth_int   = flashInt.TimeWidth(); // Width of the flash in time
                            fflashsel_abstime_int     = flashInt.AbsTime(); // Time by PMT readout clock (?)
                            fflashsel_frame_int       = flashInt.Frame(); // Frame number 
                         }
			 totflash++;
		       }// loop over flashes
		       cout<<"\n"<<"sel flash track time diff "<<fflashsel_track_dist_int<<" sel flash PE "<<fflashsel_pe_int<<"\n"<<"\n";
		       


/*                       //-------------------Looking at the waveform information----------------------------//
		       
		       //Extenral waveforms
		       
                       double wfexttime[10] = {0};
	               int  mm=0;
	               for(auto const& wfext : *WfListHandleExt)
	               {
//                       cout<<"wfext.size() "<<wfext.size()<<"\n";
	  
                         wfexttime[mm] = ( wfext.TimeStamp());
//		         std::cout << "wfext.TimeStamp() external " << wfext.TimeStamp() <<"\n";	
		         mm++;		

	               }
		       
		       // Internal waveforms
		       
	               for(auto const& wf : *WfListHandleInt)
	               {
//                       cout<<"wf.size() "<<wf.size()<<"\n";
	   
//		         std::cout << "wf.TimeStamp() " << wf.TimeStamp() << " T0 track time " << t0s.at(0)->Time()<<" trigger offset "<<TPC_trigger_offset <<" ftrack_time "<<ftrack_time<<"\n\n";

//	      	         if ( (wf.TimeStamp() - ftrack_time) < 0  ) continue;
		
	      	         fwftimeint[ftotintwf]     = wf.TimeStamp();
	      	         fwfchan[ftotintwf]        = wf.ChannelNumber();
			 fwfsel_endhitx[ftotintwf] = fsel_endhitx;
			 fwfsel_endhity[ftotintwf] = fsel_endhity;
			 fwfsel_endhitz[ftotintwf] = fsel_endhitz;
		         fwfsel_endhitkey[ftotintwf] = fsel_endhitkey;	
		         fwfsel_endwire[ftotintwf]   = fsel_endwire;
		         fwfsel_endchno[ftotintwf]   = fsel_endchno;
		         fwfsel_endtpcno[ftotintwf]  = fsel_endtpcno;
		         fwfsel_endhitchrg[ftotintwf]= fsel_endhitchrg;
		         fwfsel_endptime[ftotintwf]  = fsel_endptime;
		
//		         cout<<"abs(fwftiming - (ftrack_time)) "<<abs(fwftiming - (ftrack_time))<<"\n";	
			
              
	                 // Find the closest in time external falsh
                         double minwftime = 9999999999999, minwftime1 = 999, wfminexttime =999;
                         for(int nm =0; nm<mm; nm++)
		         {
		           minwftime1 = abs(wf.TimeStamp() - wfexttime[nm]);
		           if(minwftime1 < minwftime )
		           {
		             minwftime = minwftime1;
		             wfminexttime = wfexttime[nm];
		           }
		         }
			 
			 fwftimeext[ftotintwf]      = wfminexttime;
			 ftrkrecotime[ftotintwf]    = fsel_trkrecotime;
			 fwftime[ftotintwf]         = fwftimeint[ftotintwf] - wfminexttime;
			 fwftracktimediff[ftotintwf]= abs(fwftime[ftotintwf]-ftrkrecotime[ftotintwf]);

		
//		         double wftimerel = wftimeextint ;
//		         double timetrack  = double(t0s.at(0)->Time()/1000)-TPC_trigger_offset;
		         cout<<"trackid "<<i<<" wf channel "<<fwfchan[ftotintwf]<<" wftime "<<fwftime[ftotintwf]<<" ftrack_time "<<ftrkrecotime[ftotintwf]<<" wf and track time difference "<<fwftracktimediff[ftotintwf]<<"\n";

			 ftotintwf++;

	               }
		       
*/




 
 
 /*                       // storing selected flash information
			float flsh_trk_dist = 99999;
			float flsh_pe = -99999;
			
			totflash = 0;
    
                        for (size_t nn = 0; nn < flashlistInt.size() && nn < kMaxFlashes ; ++nn)
                        {
                          art::Ptr<recob::OpFlash> pflashInt(flashListHandleInt, nn);
                          const recob::OpFlash& flashInt = *pflashInt;
			  float min_dist = sqrt(pow((flashInt.YCenter()-fsel_endhity),2)+pow((flashInt.ZCenter()-fsel_endhitz),2));
			  
			  flash_distana[totflash] = min_dist;
			  flash_peana[totflash] = flashInt.TotalPE();
			  
			  if(min_dist< flsh_trk_dist && flashInt.TotalPE()>flsh_pe)
			  {
 	  	            flsh_trk_dist = min_dist;
			    flsh_pe = flashInt.TotalPE();
			    
			    fflashsel_track_dist_int  = flsh_trk_dist;
                            fflashsel_time_int        = flashInt.Time(); // Time relative to trigger
                            fflashsel_pe_int          = flashInt.TotalPE();
                            fflashsel_ycenter_int     = flashInt.YCenter(); // Geometric center in y
                            fflashsel_zcenter_int     = flashInt.ZCenter(); // Geometric center in z
                            fflashsel_ywidth_int      = flashInt.YWidth(); // Geometric width in y
                            fflashsel_zwidth_int      = flashInt.ZWidth(); // Geometric width in z
                            fflashsel_timewidth_int   = flashInt.TimeWidth(); // Width of the flash in time
                            fflashsel_abstime_int     = flashInt.AbsTime(); // Time by PMT readout clock (?)
                            fflashsel_frame_int       = flashInt.Frame(); // Frame number 
		  
                         }
			 totflash++;
		       }// loop over flashes
		       cout<<"\n"<<"sel flash track dist "<<fflashsel_track_dist_int<<" sel flash PE "<<fflashsel_pe_int<<"\n"<<"\n";
*/	 
		       fSelTree->Fill();
		       fsel_mu++;
		      
	             }//if(abs(shwr_dist0)< _shwr_dist)
		   }//if(fnearhitct>= _hitcountmichel)
//                 }//if(PH == _PHpass)
//                longtrk_hitkey.clear();
//	       }//if(float(resl/trkcolhits) >= 0.5)
//                }//if(trkstopz > _muendzmincut && trkstopz < _muendzmaxcut)
//	       }//if(trkstopy > _muendycut)
              }//if(MaxHitPeakTime < 4800)
	     }//if(MinHitPeakTime > 200) 
            }//if(track.Lenght()> _longtrkLen) 	     
	   }//if(count==tracklist.size()-1) //unbroken tracks
	 }//remove electron diverter bounds
        }// if(ccrosser == 1)
       }//if(IsPointInFV(fiducialBounds,recoEndPoint)) // End Point within the fiducial volume.
      }//if(IsPointInBounds(activeBounds_eff,recoStartPoint)) // Start point within the Active volume bounds.

	 
       //loop over uncontained tracks
       if((pos.X()<Xnegbound || pos.X()>Xposbound || pos.Y()<Ynegbound || pos.Y()>Yposbound || pos.Z()<Znegbound || pos.Z()>Zposbound) && 
	    (end.X()<Xnegbound || end.X()<Xposbound || end.Y()<Ynegbound || end.Y()>Yposbound || end.Z()<Znegbound || end.Z()>Zposbound) && track.Length()>100)
       {
         funcont_trks++;
       }
	  
       //if track cross two z boundaries
       if((pos.Z()<Znegbound || pos.Z()>Zposbound) && (end.Z()<Znegbound || end.Z()>Zposbound) && track.Length()>100)
       {
         fcrossZ_trks++;
       }	 

     }//if(t0s.size())
   }//if(pfps.size())
 }//Loop over all tracks 
  
 trkHitsKey.clear();
   		 
 fEventTree->Fill();
  
}//end of analyze function   
  
  /////////////////// Defintion of reset function ///////////
  void MichelStudy::reset()
  {
     
    fdau_pdg                         = -99999;
    fdau_energy                      = -99999;
       
    frun                             = -99999;
    fsubrun                          = -99999;
    fevent                           = -99999;
    fevttime                         = -99999;
    fyear_month_date                 = -99999;
    fhour_min_sec                    = -99999;
    fext_trigger_time                = -99999;
    fDriftVelocity                   = -99999;

    fall_trks                        = 0;
    fPFP_trks                        = 0;
    fT0_trks                         = 0;
    fstartinbound_trks               = 0;
    fendinFV_trks                    = 0;
    fccrosser_trks                   = 0;
    fnoelecdiv_bound                 = 0;
    funbroken_trks                   = 0; //these are unbroken stopping tracks
    flongcm_trks                     = 0;
    fminhitpt_trks                   = 0;
    fmaxhitpt_trks                   = 0;    
//    fmuendy_trks                     = 0;    
//    fmuendz_trks                     = 0;    
//    fdistmorehits_trks                   = 0;
//    fPHtest_trks                     = 0;
    fnearhits_trks                   = 0;
    fsel_mu                          = 0;
//    fseldau_mu                       = -99999;
    funcont_trks                     = 0;
    fcrossZ_trks                     = 0;
    fbackward_trks                   = 0;
    fEndinTPC_trks                   = 0;
//    ftrkcolhits                      = 0;
    ftrueMichel                      = 0;
    fno_flashes_ext                  = 0;
    fno_flashes_int                  = 0;

    fpfpana                          = 0;
    ft0ana                           = 0;
    fstartinboundana                 = 0;
    fendinFVana                      = 0;
    fccrossana                       = 0;
    fstopzana                        = 0;
    fdistana                         = 0;
    fbrokcoana                       = 0;
    ftrklenana                       = 0;
    fminhitptana                     = 0;
    fmaxhitptana                     = 0;
//    fmuendyana                       = 0;
//    fmuendzana                       = 0;
    ftrkdistcolana                       = 0;
//    fPHana                           = 0;
    fhitctana                        = 0;
    fshwrdisana                      = 0;
     
    fnearhitct                     = 0;
    fmhitcount                     = 0;
    ftrueparhitallcount            = 0;
    ftrueparhitcolcount            = 0;
    fshwrallhits                   = 0;
    
    fhitsU                         = 0;
    fhitsV                         = 0;
    fhitsY                         = 0;

    ftrkcolhits                    = 0;
    fshwrcolhits                   = 0;
    fntrkhits                      = 0;

//    fhasElect                      = 0;

      favg_missing_energy         = -99999;
      favg_missing_numelec         = -99999;  
      
      fsel_run                    = -99999;
      fsel_subrun                 = -99999;
      fsel_event                  = -99999;
      fsel_evttime                = -99999;
      fsel_endhitkey              = -99999;
      fsel_endwire                = -99999;
      fsel_endchno                = -99999;
      fsel_endtpcno               = -99999;
      fsel_endptime               = -99999;
      fsel_endhitchrg             = -99999;
      fsel_endhitx                = -99999;
      fsel_endhity                = -99999;
      fsel_endhitz                = -99999;
      fsel_ccrosser               = -99999;
      fsel_endhitmult             = -99999;
      fsel_endhitsigptime         = -99999;
      fsel_endhitsigchrg          = -99999;
      fsel_endhitsigpamp           = -99999;
      fsel_endhitdof              = -99999;
      fsel_endhitgof              = -99999;

      fsel_dist_hit_end           = -99999;
      fsel_dist_times             = -99999;
      fsel_trkstartx              = -99999;
      fsel_trkstarty              = -99999;
      fsel_trkstartz              = -99999;
      fsel_trkendx                = -99999;
      fsel_trkendy                = -99999;
      fsel_trkendz                = -99999;
      fsel_trkstartcosx           = -99999;
      fsel_trkstartcosy           = -99999;
      fsel_trkstartcosz           = -99999; 
      fsel_trkendcosx             = -99999;
      fsel_trkendcosy             = -99999;
      fsel_trkendcosz             = -99999;
      fsel_trackthetaxz           = -99999;
      fsel_trackthetayz           = -99999;
      fsel_trklen                 = -99999;
      fsel_trktheta               = -99999;
      fsel_trkphi                 = -99999;
      fsel_trkID                  = -99999;
//      fsel_PHratio                = -99999;
//      fsel_PH                     = -99999;
      fsel_minhitptime            = -99999;
      fsel_maxhitptime            = -99999;
      fsel_ncolhits               = -99999;
      fsel_nearhitcount           = -99999;
      fsel_CorrWirePtime          = -99999;
      fsel_trkrecotime            = -99999;
      fMichelcountselana          = -99999;

      fselshwr_key                = -99999;
      fselshwr_ID                 = -99999;
      fselshwr_length             = -99999;
      fselshwr_startx             = -99999;
      fselshwr_starty             = -99999;
      fselshwr_startz             = -99999;
      fselshwr_bestplane          = -99999;
      fselshwr_startdcosx         = -99999;
      fselshwr_startdcosy         = -99999;
      fselshwr_startdcosz         = -99999;
      fselshwr_openangle          = -99999;
//      fselall_shwrEnergy          = -99999;
//      fselcol_shwrEnergy          = -99999;
      fselshwr_dist               = -99999;
//      fselshwr_dEdx               = -99999;
//      fselshwr_energy             = -99999;
//      fselshwr_mipenergy          = -99999;


     // mc truth for candidate mu
      fmcsel_trkid                =-99999;
      fmcsel_vx                   =-99999;
      fmcsel_vy                   =-99999;
      fmcsel_vz                   =-99999;
      fmcsel_t                    =-99999;
      fmcsel_endx                 =-99999;
      fmcsel_endy                 =-99999;
      fmcsel_endz                 =-99999;
      fmcsel_endt                 =-99999;
      fmcsel_px                   =-99999;
      fmcsel_py                   =-99999;
      fmcsel_pz                   =-99999;
      fmcsel_momentum             =-99999;
      fmcsel_energy               =-99999;
      fmcsel_endpx                =-99999;
      fmcsel_endpy                =-99999;
      fmcsel_endpz                =-99999;
      fmcsel_endenergy            =-99999;
      fmcsel_pathlen              =-99999;
      fmcsel_length               =-99999;
      fmcsel_vxdrifted            =-99999;
      fmcsel_vydrifted            =-99999;
      fmcsel_vzdrifted            =-99999;
      fmcsel_tdrifted             =-99999;
      fmcsel_endxdrifted          =-99999;
      fmcsel_endydrifted          =-99999;
      fmcsel_endzdrifted          =-99999;
      fmcsel_endtdrifted          =-99999;
      fmcsel_pxdrifted            =-99999;
      fmcsel_pydrifted            =-99999;
      fmcsel_pzdrifted            =-99999;
      fmcsel_momentumdrifted      =-99999;
      fmcsel_energydrifted        =-99999;
      fmcsel_endpxdrifted         =-99999;
      fmcsel_endpydrifted         =-99999;
      fmcsel_endpzdrifted         =-99999;
      fmcsel_endenergydrifted     =-99999;
      fmcsel_pathlendrifted       =-99999;
      fmcsel_endprocess           =-99999;
      fmcsel_theta                =-99999;
      fmcsel_phi                  =-99999;
      fmcsel_pdg                  =-99999;
      fmcsel_status_code          =-99999;
      fmcsel_mass                 =-99999;
      fmcsel_ND                   =-99999;
      fmcsel_mother               =-99999;
      fmcsel_origin               =-99999;
      fmcsel_process              =-99999;
      fmcsel_rescatter            =-99999;

    // mc truth for true michels
      fmcd_trkid               =-99999;
      fmcd_vx                  =-99999;
      fmcd_vy                  =-99999;
      fmcd_vz                  =-99999;
      fmcd_t                   =-99999;
      fmcd_endx                =-99999;
      fmcd_endy                =-99999;
      fmcd_endz                =-99999;
      fmcd_endt                =-99999;
      fmcd_px                  =-99999;
      fmcd_py                  =-99999;
      fmcd_pz                  =-99999;
      fmcd_momentum            =-99999;
      fmcd_energy              =-99999;
      fmcd_trueselhitsE        =-99999;
      ftrueEdepo               =-99999;
      fmcd_endpx               =-99999;
      fmcd_endpy               =-99999;
      fmcd_endpz               =-99999;
      fmcd_endenergy           =-99999;
      fmcd_pathlen             =-99999;
      fmcd_vxdrifted           =-99999;
      fmcd_vydrifted           =-99999;
      fmcd_vzdrifted           =-99999;
      fmcd_tdrifted            =-99999;
      fmcd_endxdrifted         =-99999;
      fmcd_endydrifted         =-99999;
      fmcd_endzdrifted         =-99999;
      fmcd_endtdrifted         =-99999;
      fmcd_pxdrifted           =-99999;
      fmcd_pydrifted           =-99999;
      fmcd_pzdrifted           =-99999;
      fmcd_momentumdrifted     =-99999;
      fmcd_energydrifted       =-99999;
      fmcd_endpxdrifted        =-99999;
      fmcd_endpydrifted        =-99999;
      fmcd_endpzdrifted        =-99999;
      fmcd_endenergydrifted    =-99999;
      fmcd_pathlendrifted      =-99999;
      fmcd_endprocess          =-99999;
      fmcd_theta               =-99999;
      fmcd_phi                 =-99999;
      fmcd_pdg                 =-99999;
      fmcd_status_code         =-99999;
      fmcd_mass                =-99999;
      fmcd_ND                  =-99999;
      fmcd_mother              =-99999;
      fmcd_origin              =-99999;
      fmcd_process             =-99999;
      fmcd_rescatter           =-99999;
     
/*    fflash_reco_time_diff[i]       =-99999;
      fflash_time_sel[i]             =-99999;
      fflash_pe_sel[i]               =-99999;
      fflash_ycenter_sel[i]          =-99999;
      fflash_zcenter_sel[i]          =-99999;
      fflash_ywidth_sel[i]           =-99999;
      fflash_zwidth_sel[i]           =-99999;
      fflash_timewidth_sel[i]        =-99999;
      fflash_abstime_sel[i]          =-99999;
      fflash_frame_sel[i]            =-99999;
      fflash_time_wrttrigger_sel[i]  =-99999;
*/ 
      //MC truth for nearby michel hits
      fmchits_trkid               =-99999;
      fmchits_vx                  =-99999;
      fmchits_vy                  =-99999;
      fmchits_vz                  =-99999;
      fmchits_t                   =-99999;
      fmchits_endx                =-99999;
      fmchits_endy                =-99999;
      fmchits_endz                =-99999;
      fmchits_endt                =-99999;
      fmchits_px                  =-99999;
      fmchits_py                  =-99999;
      fmchits_pz                  =-99999;
      fmchits_momentum            =-99999;
      fmchits_energy              =-99999;
      fmchits_endpx               =-99999;
      fmchits_endpy               =-99999;
      fmchits_endpz               =-99999;
      fmchits_endenergy           =-99999;
      fmchits_pathlen             =-99999;
      fmchits_vxdrifted           =-99999;
      fmchits_vydrifted           =-99999;
      fmchits_vzdrifted           =-99999;
      fmchits_tdrifted            =-99999;
      fmchits_endxdrifted         =-99999;
      fmchits_endydrifted         =-99999;
      fmchits_endzdrifted         =-99999;
      fmchits_endtdrifted         =-99999;
      fmchits_pxdrifted           =-99999;
      fmchits_pydrifted           =-99999;
      fmchits_pzdrifted           =-99999;
      fmchits_momentumdrifted     =-99999;
      fmchits_energydrifted       =-99999;
      fmchits_endpxdrifted        =-99999;
      fmchits_endpydrifted        =-99999;
      fmchits_endpzdrifted        =-99999;
      fmchits_endenergydrifted    =-99999;
      fmchits_pathlendrifted      =-99999;
      fmchits_endprocess          =-99999;
      fmchits_theta               =-99999;
      fmchits_phi                 =-99999;
      fmchits_pdg                 =-99999;
      fmchits_status_code         =-99999;
      fmchits_mass                =-99999;
      fmchits_ND                  =-99999;
      fmchits_mother              =-99999;
      fmchits_origin              =-99999;
      fmchits_process             =-99999;
      fmchits_rescatter           =-99999;

      //MC truth for michel cone hits
      fmcconehits_trkid               =-99999;
      fmcconehits_vx                  =-99999;
      fmcconehits_vy                  =-99999;
      fmcconehits_vz                  =-99999;
      fmcconehits_t                   =-99999;
      fmcconehits_endx                =-99999;
      fmcconehits_endy                =-99999;
      fmcconehits_endz                =-99999;
      fmcconehits_endt                =-99999;
      fmcconehits_px                  =-99999;
      fmcconehits_py                  =-99999;
      fmcconehits_pz                  =-99999;
      fmcconehits_momentum            =-99999;
      fmcconehits_energy              =-99999;
      fmcconehits_endpx               =-99999;
      fmcconehits_endpy               =-99999;
      fmcconehits_endpz               =-99999;
      fmcconehits_endenergy           =-99999;
      fmcconehits_pathlen             =-99999;
      fmcconehits_vxdrifted           =-99999;
      fmcconehits_vydrifted           =-99999;
      fmcconehits_vzdrifted           =-99999;
      fmcconehits_tdrifted            =-99999;
      fmcconehits_endxdrifted         =-99999;
      fmcconehits_endydrifted         =-99999;
      fmcconehits_endzdrifted         =-99999;
      fmcconehits_endtdrifted         =-99999;
      fmcconehits_pxdrifted           =-99999;
      fmcconehits_pydrifted           =-99999;
      fmcconehits_pzdrifted           =-99999;
      fmcconehits_momentumdrifted     =-99999;
      fmcconehits_energydrifted       =-99999;
      fmcconehits_endpxdrifted        =-99999;
      fmcconehits_endpydrifted        =-99999;
      fmcconehits_endpzdrifted        =-99999;
      fmcconehits_endenergydrifted    =-99999;
      fmcconehits_pathlendrifted      =-99999;
      fmcconehits_endprocess          =-99999;
      fmcconehits_theta               =-99999;
      fmcconehits_phi                 =-99999;
      fmcconehits_pdg                 =-99999;
      fmcconehits_status_code         =-99999;
      fmcconehits_mass                =-99999;
      fmcconehits_ND                  =-99999;
      fmcconehits_mother              =-99999;
      fmcconehits_origin              =-99999;
      fmcconehits_process             =-99999;
      fmcconehits_rescatter           =-99999;

      //MC truth for showers
      fmcshwr_trkid               =-99999;
      fmcshwr_vx                  =-99999;
      fmcshwr_vy                  =-99999;
      fmcshwr_vz                  =-99999;
      fmcshwr_t                   =-99999;
      fmcshwr_endx                =-99999;
      fmcshwr_endy                =-99999;
      fmcshwr_endz                =-99999;
      fmcshwr_endt                =-99999;
      fmcshwr_px                  =-99999;
      fmcshwr_py                  =-99999;
      fmcshwr_pz                  =-99999;
      fmcshwr_momentum            =-99999;
      fmcshwr_energy              =-99999;
      fmcshwr_endpx               =-99999;
      fmcshwr_endpy               =-99999;
      fmcshwr_endpz               =-99999;
      fmcshwr_endenergy           =-99999;
      fmcshwr_pathlen             =-99999;
      fmcshwr_vxdrifted           =-99999;
      fmcshwr_vydrifted           =-99999;
      fmcshwr_vzdrifted           =-99999;
      fmcshwr_tdrifted            =-99999;
      fmcshwr_endxdrifted         =-99999;
      fmcshwr_endydrifted         =-99999;
      fmcshwr_endzdrifted         =-99999;
      fmcshwr_endtdrifted         =-99999;
      fmcshwr_pxdrifted           =-99999;
      fmcshwr_pydrifted           =-99999;
      fmcshwr_pzdrifted           =-99999;
      fmcshwr_momentumdrifted     =-99999;
      fmcshwr_energydrifted       =-99999;
      fmcshwr_endpxdrifted        =-99999;
      fmcshwr_endpydrifted        =-99999;
      fmcshwr_endpzdrifted        =-99999;
      fmcshwr_endenergydrifted    =-99999;
      fmcshwr_pathlendrifted      =-99999;
      fmcshwr_endprocess          =-99999;
      fmcshwr_theta               =-99999;
      fmcshwr_phi                 =-99999;
      fmcshwr_pdg                 =-99999;
      fmcshwr_status_code         =-99999;
      fmcshwr_mass                =-99999;
      fmcshwr_ND                  =-99999;
      fmcshwr_mother              =-99999;
      fmcshwr_origin              =-99999;
      fmcshwr_process             =-99999;
      fmcshwr_rescatter           =-99999;

      fflashsel_track_dist_int    =-99999;	       
      fflashsel_time_int          =-99999;
      fflashsel_pe_int            =-99999;
      fflashsel_ycenter_int       =-99999;
      fflashsel_zcenter_int       =-99999;
      fflashsel_ywidth_int	  =-99999;
      fflashsel_zwidth_int        =-99999;
      fflashsel_timewidth_int     =-99999;
      fflashsel_abstime_int	  =-99999;
      fflashsel_frame_int         =-99999;
      
      ftotintwf                   =-99999;

    for(int i=0; i<kMaxFlashes; i++)
    {
      //external flashes
      fflash_time_ext[i]             = -99999;
      fflash_pe_ext[i]               = -99999;
      fflash_ycenter_ext[i]          = -99999;
      fflash_zcenter_ext[i]          = -99999;
      fflash_ywidth_ext[i]           = -99999;
      fflash_zwidth_ext[i]           = -99999;
      fflash_timewidth_ext[i]        = -99999;
      fflash_abstime_ext[i]          = -99999;
      fflash_frame_ext[i]            = -99999;

      //internal flashes
      fflash_time_int[i]             = -99999;
      fflash_pe_int[i]               = -99999;
      fflash_ycenter_int[i]          = -99999;
      fflash_zcenter_int[i]          = -99999;
      fflash_ywidth_int[i]           = -99999;
      fflash_zwidth_int[i]           = -99999;
      fflash_timewidth_int[i]        = -99999;
      fflash_abstime_int[i]          = -99999;
      fflash_frame_int[i]            = -99999;
    }  
    
        

    for(int i=0; i<kMaxTracks; i++)
    {
      fpfpsana[i]                    = -99999;
      ft0sana[i]                     = -99999;
      fstartinboundsana[i]           = -99999;
      fendinFVsana[i]               = -99999;
      fccrosserana[i]                = -99999;
      felecdivstopzana[i]            = -99999;
      fdistanceana[i]                = -99999;
      fbrokencountana[i]             = -99999;
      fminhitptimeana[i]             = -99999;
      fmaxhitptimeana[i]             = -99999;      
//      fmuonendyana[i]                = -99999;      
//      fmuonendzana[i]                = -99999;      
      ftrklengthana[i]               = -99999;
//      ftrkdistcollhitsana[i]             = -99999;
//      fPHtestana[i]                  = -99999;
      fnearhitcountana[i]            = -99999;
      fnshwrdistana[i]               = -99999;
      fshwrdistana[i]                = -99999;

      fMichelcountpfpana[i]          = -99999;
      fMichelcountt0ana[i]           = -99999;
      fMichelcountstartinboundana[i] = -99999;
      fMichelcountendinFVana[i]      = -99999;
      fMichelcountccrosserana[i]     = -99999;
      fMichelcountelecdivstopzana[i] = -99999;
//      fMichelcountdistana[i]         = -99999;
      fMichelcountbrokencountana[i]  = -99999;
      fMichelcountminhitptimeana[i]  = -99999;
      fMichelcountmaxhitptimeana[i]  = -99999;
//      fMichelcountmuonendyana[i]     = -99999;
//      fMichelcountmuonendzana[i]     = -99999;
      fMichelcountlenana[i]          = -99999;
//      fMichelcountcollana[i]         = -99999;
//      fMichelcountPHtestana[i]       = -99999;
//      fMichelcountdistcollana[i]      = -99999;
      fMichelcountnearhitana[i]      = -99999;
      fMichelcountshwrdistana[i]     = -99999;

      ftrueEpfpana[i]          = -99999;
      ftrueEt0ana[i]           = -99999;
      ftrueEstartinboundana[i] = -99999;
      ftrueEendinFVana[i]      = -99999;
      ftrueEccrosserana[i]     = -99999;
      ftrueEelecdivstopzana[i] = -99999;
//      ftrueEdistana[i]         = -99999;
      ftrueEbrokencountana[i]  = -99999;
      ftrueEminhitptimeana[i]  = -99999;
      ftrueEmaxhitptimeana[i]  = -99999;
//      fMichelcountmuonendyana[i]     = -99999;
//      fMichelcountmuonendzana[i]     = -99999;
      ftrueElenana[i]          = -99999;
//      fMichelcountcollana[i]         = -99999;
//      fMichelcountPHtestana[i]       = -99999;
//      fMichelcountdistcollana[i]      = -99999;
      ftrueEnearhitana[i]	= -99999;
      ftrueEshwrdistana[i]     = -99999;
 

    }
      
   
    for(int i=0; i<kMaxHits; i++)
    {
 			   

      fhits_key[i]       = -99999;
      fhits_charge[i]    = -99999;
      fhits_wire[i]      = -99999;
      fhits_peakT[i]     = -99999;
      fhits_TPC[i]       = -99999;
      fhits_chno[i]      = -99999;
      fhits_xpos[i]      = -99999;
      fhits_ypos[i]      = -99999;
      fhits_zpos[i]      = -99999;
      fhits_mult[i]      = -99999;
      fhits_sigptime[i]  = -99999;
      fhits_sigchrg[i]   = -99999;
      fhits_sigpamp[i]    = -99999;
      fhits_dof[i]       = -99999;
      fhits_gof[i]       = -99999;
      fhits_ptplusRMS[i]     = -99999;
      fhits_ptminusRMS[i]    = -99999;
      fhits_cnnMichel[i]     = -99999;
      fhits_cnnEM[i]         = -99999;
      fhits_cnnTrack[i]      = -99999;

      ftrkdqdxU[i]               = -99999;
      ftrkdedxU[i]               = -99999;
      ftrkresrangeU[i]               = -99999;
      ftrkhitxU[i]               = -99999;
      ftrkhityU[i]               = -99999;
      ftrkhitzU[i]               = -99999;
      ftrkpitchU[i]               = -99999;
      ftrkdqdxV[i]               = -99999;
      ftrkdedxV[i]               = -99999;
      ftrkresrangeV[i]               = -99999;
      ftrkhitxV[i]               = -99999;
      ftrkhityV[i]               = -99999;
      ftrkhitzV[i]               = -99999;
      ftrkpitchV[i]               = -99999;
      ftrkdqdxY[i]               = -99999;
      ftrkdedxY[i]               = -99999;
      ftrkresrangeY[i]               = -99999;
      ftrkhitxY[i]		= -99999;
      ftrkhityY[i]		 = -99999;
      ftrkhitzY[i]               = -99999;
      ftrkpitchY[i]               = -99999;

      fnearhits_key[i]               = -99999;
      fnearhits_chno[i]              = -99999;
      fnearhits_peakT[i]             = -99999;
      fnearhits_charge[i]            = -99999;
      fnearhits_wire[i]              = -99999;
      fnearhits_TPC[i]               = -99999;
      fnearhits_plane[i]             = -99999;
      fnearhits_xpos[i]              = -99999;
      fnearhits_ypos[i]              = -99999;
      fnearhits_zpos[i]              = -99999;
      fnearhits_mult[i]              = -99999;
      fnearhits_sigptime[i]          = -99999;
      fnearhits_sigchrg[i]           = -99999;
      fnearhits_sigpamp[i]            = -99999;
      fnearhits_dof[i]               = -99999;
      fnearhits_gof[i]               = -99999;
      fnearhits_ptplusRMS[i]         = -99999;
      fnearhits_ptminusRMS[i]        = -99999;
      fnearhits_cnnMichel[i]         = -99999;
      fnearhits_cnnEM[i]             = -99999;
      fnearhits_cnnTrack[i]          = -99999;

      ftrueparhitsall_key[i]                  = -99999;
      ftrueparhitsall_chno[i]                 = -99999;
      ftrueparhitsall_peakT[i]                = -99999;
      ftrueparhitsall_charge[i]               = -99999;
      ftrueparhitsall_wire[i]                 = -99999;
      ftrueparhitsall_TPC[i]                  = -99999;
      ftrueparhitsall_plane[i]                = -99999;
      ftrueparhitsall_xpos[i]                 = -99999;
      ftrueparhitsall_ypos[i]                 = -99999;
      ftrueparhitsall_zpos[i]                 = -99999;
      ftrueparhitsall_mult[i]              = -99999;
      ftrueparhitsall_sigptime[i]          = -99999;
      ftrueparhitsall_sigchrg[i]           = -99999;
      ftrueparhitsall_sigpamp[i]            = -99999;
      ftrueparhitsall_dof[i]		   = -99999;
      ftrueparhitsall_gof[i]               = -99999;

      ftrueMiEFrac[i]                         = -99999;
      ftrueparhitscol_key[i]                  = -99999;
      ftrueparhitscol_chno[i]                 = -99999;
      ftrueparhitscol_peakT[i]                = -99999;
      ftrueparhitscol_charge[i]               = -99999;
      ftrueparhitscol_wire[i]                 = -99999;
      ftrueparhitscol_TPC[i]                  = -99999;
      ftrueparhitscol_plane[i]                = -99999;
      ftrueparhitscol_xpos[i]                 = -99999;
      ftrueparhitscol_ypos[i]                 = -99999;
      ftrueparhitscol_zpos[i]                 = -99999;
      ftrueparhitscol_mult[i]              = -99999;
      ftrueparhitscol_sigptime[i]          = -99999;
      ftrueparhitscol_sigchrg[i]           = -99999;
      ftrueparhitscol_sigpamp[i]            = -99999;
      ftrueparhitscol_dof[i]		   = -99999;
      ftrueparhitscol_gof[i]               = -99999;
      ftrueparhitscol_angledeg[i]             = -99999;
      ftrueparhitscol_maghitveccostheta[i]    = -99999;
      ftrueparhitscol_distance[i]             = -99999;

        fshwrhits_chno[i]              = -99999;
        fshwrhits_peakT[i]             = -99999;
        fshwrhits_charge[i]            = -99999;
        fshwrhits_wire[i] 	       = -99999;
        fshwrhits_TPC[i]               = -99999;
        fshwrhits_plane[i]	       = -99999;
        fshwrhits_xpos[i]              = -99999;
        fshwrhits_ypos[i] 	       = -99999;
        fshwrhits_zpos[i] 	       = -99999;
        fshwrhits_mult[i]              = -99999;
        fshwrhits_sigptime[i]          = -99999;
        fshwrhits_sigchrg[i]           = -99999;
        fshwrhits_sigpamp[i]            = -99999;
        fshwrhits_dof[i]	       = -99999;
        fshwrhits_gof[i]               = -99999;

      fshwrallhits_chno[i]           = -99999;
      fshwrallhits_peakT[i]          = -99999;
      fshwrallhits_charge[i]         = -99999;
      fshwrallhits_wire[i] 	     = -99999;
      fshwrallhits_TPC[i]            = -99999;
      fshwrallhits_plane[i]	     = -99999;
      fshwrallhits_xpos[i]           = -99999;
      fshwrallhits_ypos[i] 	     = -99999;
      fshwrallhits_zpos[i] 	     = -99999;
      fshwrallhits_mult[i]           = -99999;
      fshwrallhits_sigptime[i]       = -99999;
      fshwrallhits_sigchrg[i]        = -99999;
      fshwrallhits_sigpamp[i]         = -99999;
      fshwrallhits_dof[i]	     = -99999;
      fshwrallhits_gof[i]            = -99999;
    
      fmhits_key[i]                  = -99999;
      fmhits_chno[i]                 = -99999;
      fmhits_peakT[i]                = -99999;
      fmhits_charge[i]               = -99999;
      fmhits_wire[i]                 = -99999;
      fmhits_TPC[i]                  = -99999;
      fmhits_plane[i]                = -99999;
      fmhits_xpos[i]                 = -99999;
      fmhits_ypos[i]                 = -99999;
      fmhits_zpos[i]                 = -99999;
      fmhits_mult[i]                 = -99999;
      fmhits_sigptime[i]             = -99999;
      fmhits_sigchrg[i]              = -99999;
      fmhits_sigpamp[i]               = -99999;
      fmhits_dof[i]	             = -99999;
      fmhits_gof[i]	             = -99999;
      fmhits_angledeg[i]             = -99999;
      fmhits_maghitveccostheta[i]    = -99999;
      fmhits_distance[i]             = -99999;
      fmhits_longtrk[i]              = -99999;
      fmhits_sametrk[i]              = -99999;
      fmhits_corrhit[i]              = -99999;
      fmhits_ptplusRMS[i]            = -99999;
      fmhits_ptminusRMS[i]           = -99999;
      fmhits_cnnMichel[i]            = -99999;
      fmhits_cnnEM[i]                = -99999;
      fmhits_cnnTrack[i]             = -99999;
    }
    

    
    totflash = 0;
    for(int i = 0; i<kMaxFlashes; i++)
    {
	flash_distana[i] = -99999;
	flash_peana[i]   = -99999;
    }
    
 
   for(int i = 0; i<kMaxWf; i++)
    {
	fwftimeint[i]                 = -99999;
	fwfchan[i]                    = -99999;
	fwfsel_endhitx[i]             = -99999;
	fwfsel_endhity[i]             = -99999;
	fwfsel_endhitz[i]             = -99999;
	fwfsel_endhitkey[i]           = -99999;
	fwfsel_endwire[i]             = -99999;
	fwfsel_endchno[i]             = -99999;
	fwfsel_endtpcno[i]            = -99999;
	fwfsel_endhitchrg[i]          = -99999;
	fwfsel_endptime[i]            = -99999;
        fwfsel_endhitmult[i]          = -99999;
        fwfsel_endhitsigptime[i]      = -99999;
        fwfsel_endhitsigchrg[i]       = -99999;
        fwfsel_endhitsigpamp[i]       = -99999;
        fwfsel_endhitdof[i]	      = -99999;
        fwfsel_endhitgof[i]	      = -99999;
	fwftimeext[i]                 = -99999;
	ftrkrecotime[i]               = -99999;
	fwftime[i]                    = -99999;
	fwftracktimediff[i]           = -99999;            
	
    }

    fmichel_conesize = -99999;
    fmichel_wcount = -99999;
    for(int i = 0; i<kMaxct; i++)
    {
      
      fmichel_zpos[i]          = -99999;
      fmichel_ypos[i]          = -99999;
      fmichel_xpos[i]          = -99999;
      fmichel_chrg[i]          = -99999;
      fmichel_chno[i]          = -99999;
      fmichel_key[i]          = -99999;
      fmichel_wire[i]          = -99999;
      fmichel_chargehit[i]     = -99999;
      fmichel_tpc[i]          = -99999;
      fmichel_ptime[i]          = -99999;
      fmichel_angledeg[i]	   = -99999;
      fmichel_maghitveccostheta[i]	   = -99999;
      fmichel_distance[i]          = -99999;
      fmichel_mult[i]          = -99999;
      fmichel_sigptime[i]	   = -99999;
      fmichel_sigchrg[i]	  = -99999;
      fmichel_sigpamp[i]          = -99999;
      fmichel_dof[i]	      = -99999;
      fmichel_gof[i]	      = -99999;
      fmichel_ptminusRMS[i]         = -99999;
      fmichel_ptplusRMS[i]	    = -99999;
      fmichel_status[i]             = -99999;
      fmichel_cnnMichel[i]          = -99999;
      fmichel_cnnEM[i]              = -99999;
      fmichel_cnnTrack[i]           = -99999;
    }
    
  }  
    //////////////////////// End of definition ///////////////	
	  
    DEFINE_ART_MODULE(MichelStudy)
  }



