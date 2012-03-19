 
#include "PhysicsTools/SelectorUtils/interface/WPlusJetsEventSelector.h"
#include "PhysicsTools/SelectorUtils/interface/SimpleCutBasedElectronIDSelectionFunctor.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"

#include <iostream>
#include <string>                                                                                    
#include "TRegexp.h" 

using namespace std;

WPlusJetsEventSelector::WPlusJetsEventSelector( edm::ParameterSet const & params ) :
  EventSelector(),
  muonTag_         (params.getParameter<edm::InputTag>("muonSrc") ),
  electronTag_     (params.getParameter<edm::InputTag>("electronSrc") ),
  jetTag_          (params.getParameter<edm::InputTag>("jetSrc") ),
  metTag_          (params.getParameter<edm::InputTag>("metSrc") ),  
  trigTag_         (params.getParameter<edm::InputTag>("trigSrc") ),
  muTrig_          (params.getParameter<std::string>("muTrig")),
  eleTrig_         (params.getParameter<std::string>("eleTrig")),
  pvSelector_      (params.getParameter<edm::ParameterSet>("pvSelector") ),
  muonIdTight_     (params.getParameter<edm::ParameterSet>("muonIdTight") ),
  electronIdTight_ (params.getParameter<edm::ParameterSet>("electronIdTight") ),
  muonIdLoose_     (params.getParameter<edm::ParameterSet>("muonIdLoose") ),
  electronIdLoose_ (params.getParameter<edm::ParameterSet>("electronIdLoose") ),
  jetIdLoose_      (params.getParameter<edm::ParameterSet>("jetIdLoose") ),
  pfjetIdLoose_    (params.getParameter<edm::ParameterSet>("pfjetIdLoose") ),
  jecUnc_          (0),
  JECUncertaintyFile_(params.getParameter<std::string>("JECUncertaintyFile")),
  fancyJES_        ("none"),
  flatAdditionalUncer_(params.getParameter<double>("flatAdditionalUncer")),
  minJets_         (params.getParameter<int> ("minJets") ),
  muJetDRJets_     (params.getParameter<double>("muJetDRJets")),
  muJetDRMuon_     (params.getParameter<double>("muJetDRMuon")),
  eleJetDR_        (params.getParameter<double>("eleJetDR")),
  muPlusJets_      (params.getParameter<bool>("muPlusJets") ),
  ePlusJets_       (params.getParameter<bool>("ePlusJets") ),
  isHWW_           (params.getParameter<bool>("isHWW") ),
  isLoose_           (params.getParameter<bool>("isLoose") ),
  muPtMin_         (params.getParameter<double>("muPtMin")), 
  muEtaMax_        (params.getParameter<double>("muEtaMax")), 
  eleEtMin_        (params.getParameter<double>("eleEtMin")), 
  eleEtaMax_       (params.getParameter<double>("eleEtaMax")), 
  muPtMinLoose_    (params.getParameter<double>("muPtMinLoose")), 
  muEtaMaxLoose_   (params.getParameter<double>("muEtaMaxLoose")), 
  eleEtMinLoose_   (params.getParameter<double>("eleEtMinLoose")), 
  eleEtaMaxLoose_  (params.getParameter<double>("eleEtaMaxLoose")), 
  jetPtMin_        (params.getParameter<double>("jetPtMin")), 
  jetEtaMax_       (params.getParameter<double>("jetEtaMax")), 
  jetScale_        (params.getParameter<double>("jetScale")),
  jerFactor_        (params.getParameter<double>("jerFactor")),
  metMin_          (params.getParameter<double>("metMin"))
{
  // make the bitset
  push_back( "Inclusive"      );
  push_back( "Trigger"        );
  push_back( "PV"             );
  push_back( ">= 1 Lepton"    );
  push_back( "== 1 Tight Lepton"    );
  push_back( "== 1 Tight Lepton, Mu Veto");
  push_back( "== 1 Lepton"    );
  push_back( "MET Cut"        );
  push_back( "Z Veto"         );
  push_back( "Conversion Veto");
  push_back( "Cosmic Veto"    );
  push_back( ">=1 Jets"       );
  push_back( ">=2 Jets"       );
  push_back( ">=3 Jets"       );
  push_back( ">=4 Jets"       );
  push_back( ">=5 Jets"       );


  // turn (almost) everything on by default
  set( "Inclusive"      );
  set( "Trigger"        );
  set( "PV"             );
  set( ">= 1 Lepton"    );
  set( "== 1 Tight Lepton"    );
  set( "== 1 Tight Lepton, Mu Veto");
  set( "== 1 Lepton"    );
  set( "MET Cut"        );
  set( "Z Veto"         );
  set( "Conversion Veto");
  set( "Cosmic Veto"    );
  set( ">=1 Jets", minJets_ >= 1);
  set( ">=2 Jets", minJets_ >= 2);
  set( ">=3 Jets", minJets_ >= 3);
  set( ">=4 Jets", minJets_ >= 4);
  set( ">=5 Jets", minJets_ >= 5); 


  inclusiveIndex_ = index_type(&bits_, std::string("Inclusive"      ));
  triggerIndex_ = index_type(&bits_, std::string("Trigger"        ));
  pvIndex_ = index_type(&bits_, std::string("PV"             ));
  lep1Index_ = index_type(&bits_, std::string(">= 1 Lepton"    ));
  lep2Index_ = index_type(&bits_, std::string("== 1 Tight Lepton"    ));
  lep3Index_ = index_type(&bits_, std::string("== 1 Tight Lepton, Mu Veto"));
  lep4Index_ = index_type(&bits_, std::string("== 1 Lepton"    ));
  metIndex_ = index_type(&bits_, std::string("MET Cut"        ));
  zvetoIndex_ = index_type(&bits_, std::string("Z Veto"         ));
  conversionIndex_ = index_type(&bits_, std::string("Conversion Veto"));
  cosmicIndex_ = index_type(&bits_, std::string("Cosmic Veto"    ));
  jet1Index_ = index_type(&bits_, std::string(">=1 Jets"));
  jet2Index_ = index_type(&bits_, std::string(">=2 Jets"));
  jet3Index_ = index_type(&bits_, std::string(">=3 Jets"));
  jet4Index_ = index_type(&bits_, std::string(">=4 Jets"));
  jet5Index_ = index_type(&bits_, std::string(">=5 Jets")); 

  if ( params.exists("cutsToIgnore") )
    setIgnoredCuts( params.getParameter<std::vector<std::string> >("cutsToIgnore") );
  
  // See if the config file requests a valid form of fancy JES.
  // Check to see that the request is for a valid type
  if ( params.exists("fancyJES") ){
    std::string configFancyJES = params.getParameter< std::string >("fancyJES");
    if (configFancyJES == "up" ||
        configFancyJES == "down" ||
        configFancyJES == "none") {
      fancyJES_ = configFancyJES;      
      std::cout << "FANCY JES ---- using setting      " << fancyJES_ << std::endl;      
    } else {
      std::cout << "FANCY JES ---- you requested setting     "  << fancyJES_ << std::endl
                << "   but I don't know how to use that seting, exiting" << std:: endl;

      exit(33);
      
    }
    if (fancyJES_ == "up" || fancyJES_ == "down") {
      jecUnc_ = new JetCorrectionUncertainty(JECUncertaintyFile_);
    }
  }
  
	

  retInternal_ = getBitTemplate();
}

bool WPlusJetsEventSelector::operator() ( edm::EventBase const & event, pat::strbitset & ret)
{

  ret.set(false);

  selectedJets_.clear();
  cleanedJets_.clear();
  temporaryMuons_.clear();
  selectedMuons_.clear();
  selectedElectrons_.clear();
  looseMuons_.clear();
  looseElectrons_.clear();
  selectedMETs_.clear();


  passCut( ret, inclusiveIndex_);


  bool passTrig = false;
  if (!ignoreCut(triggerIndex_) ) {

    edm::Handle<pat::TriggerEvent> triggerEvent;
    event.getByLabel(trigTag_, triggerEvent);

    pat::TriggerEvent const * trig = &*triggerEvent;

    if ( trig->wasRun() && trig->wasAccept() ) {

                                                                                            
      //match muon trigger names to our wild card expression                                      
      TRegexp  muTrigRegexp(muTrig_);                                                             
      bool matchMuTrigName = false;                                                                 
      for(pat::TriggerPathCollection::const_iterator iPath = trig->paths()->begin(); iPath != trig->paths()->end(); ++iPath){ 
	    TString thisTrigPath = iPath->name();                                                     
	  matchMuTrigName =  thisTrigPath.Contains(muTrigRegexp);                                     
	  if(matchMuTrigName == true){                                                                
	    pat::TriggerPath const * muPath = trig->path(iPath->name());                            
	    if ( muPlusJets_ && muPath != 0 && muPath->wasAccept() ) {                             
              passTrig = true;                                                                    
            }                                                                                     
	  }                                                                                         
    }// full Trigger path collection     

     //match electron trigger names to our wild card expression 
     TRegexp  eleTrigRegexp(eleTrig_);
     bool matchElTrigName = false;
     for(pat::TriggerPathCollection::const_iterator iPath = trig->paths()->begin(); iPath != trig->paths()->end(); ++iPath){
	TString thisTrigPath = iPath->name();
	matchElTrigName =  thisTrigPath.Contains(eleTrigRegexp);
	if(matchElTrigName == true){
	  pat::TriggerPath const * elePath = trig->path(iPath->name());
	  if ( ePlusJets_ && elePath != 0 && elePath->wasAccept() ) {
	    passTrig = true;
	  }
	}
     }//end loop over full trigger paths collection
    }
  }
  
  if ( ignoreCut(triggerIndex_) || 
       passTrig ) {
    passCut(ret, triggerIndex_);


    bool passPV = false;

    // This way is all broken don't do it
    // passPV = pvSelector_( event );

    passPV = pvSelector_.passPVSelection(event);
    
    if ( ignoreCut(pvIndex_) || passPV ) {
      passCut(ret, pvIndex_);
  
      edm::Handle< vector< pat::Electron > > electronHandle;
      event.getByLabel (electronTag_, electronHandle);
  
      edm::Handle< vector< pat::Muon > > muonHandle;
      event.getByLabel (muonTag_, muonHandle);

      edm::Handle< vector< pat::Jet > > jetHandle;

      edm::Handle< edm::OwnVector<reco::Candidate> > jetClonesHandle ;

      edm::Handle< vector< pat::MET > > metHandle;
      event.getByLabel (metTag_, metHandle);

      int nElectrons = 0;
      for ( std::vector<pat::Electron>::const_iterator electronBegin = electronHandle->begin(),
              electronEnd = electronHandle->end(), ielectron = electronBegin;
            ielectron != electronEnd; ++ielectron ) {
        ++nElectrons;
        // Tight cuts

	if(isHWW_){
	  bool passIDLoose = false;
	  bool passID = false;
	  
	  // Calculate 
	  reco::TrackBase::Point beamPoint(0,0,0);
	  reco::BeamSpot beamSpot;
	  edm::Handle<reco::BeamSpot> beamSpotHandle;
	  event.getByLabel(edm::InputTag("offlineBeamSpot"), beamSpotHandle);
	  if( beamSpotHandle.isValid() ){
	    beamSpot = *beamSpotHandle;
	  } else{
	    edm::LogError("DataNotAvailable")
	      << "No beam spot available from EventSetup, not adding high level selection \n";
	  }
	  beamPoint = reco::TrackBase::Point ( beamSpot.x0(), beamSpot.y0(), beamSpot.z0() );
	  double d0bs = ielectron->gsfTrack()->dxy(beamPoint);
	  

	  SimpleCutBasedElectronIDSelectionFunctor patSele95(SimpleCutBasedElectronIDSelectionFunctor::trkIso95);
	  SimpleCutBasedElectronIDSelectionFunctor patSele70(SimpleCutBasedElectronIDSelectionFunctor::trkIso70);
	  SimpleCutBasedElectronIDSelectionFunctor patSele95comb(SimpleCutBasedElectronIDSelectionFunctor::cIso95);
	  SimpleCutBasedElectronIDSelectionFunctor patSele70comb(SimpleCutBasedElectronIDSelectionFunctor::cIso70);
	  
	  if(isLoose_){
	    passIDLoose = patSele95(*ielectron,event);
	    passID = patSele70(*ielectron,event);
	  }
	  else{
	    passIDLoose = patSele95comb(*ielectron,event);
	    passID = patSele70comb(*ielectron,event);
	  }
	  if ( ielectron->et() > eleEtMin_ && fabs(ielectron->eta()) < eleEtaMax_ && 
	       (fabs(ielectron->superCluster()->eta()) < 1.4442 || fabs(ielectron->superCluster()->eta()) > 1.5660) &&
	       passID && fabs(d0bs) < 0.02) {
	    selectedElectrons_.push_back( reco::ShallowClonePtrCandidate( edm::Ptr<pat::Electron>( electronHandle, ielectron - electronBegin ) ) );
	  } else {
	    // Loose cuts
	    if ( ielectron->et() > eleEtMinLoose_ && fabs(ielectron->eta()) < eleEtaMaxLoose_ && 
		 (fabs(ielectron->superCluster()->eta()) < 1.4442 || fabs(ielectron->superCluster()->eta()) > 1.5660) &&
		 passIDLoose) {
	      looseElectrons_.push_back( reco::ShallowClonePtrCandidate( edm::Ptr<pat::Electron>( electronHandle, ielectron - electronBegin ) ) );
	    }
	  }
	}
	else {
	  if ( ielectron->et() > eleEtMin_ && fabs(ielectron->eta()) < eleEtaMax_ && 
	       electronIdTight_(*ielectron) &&
	       ielectron->electronID( "eidRobustTight" ) > 0  ) {
	    selectedElectrons_.push_back( reco::ShallowClonePtrCandidate( edm::Ptr<pat::Electron>( electronHandle, ielectron - electronBegin ) ) );
	  } else {
	    // Loose cuts
	    if ( ielectron->et() > eleEtMinLoose_ && fabs(ielectron->eta()) < eleEtaMaxLoose_ && 
		 electronIdLoose_(*ielectron) ) {
	      looseElectrons_.push_back( reco::ShallowClonePtrCandidate( edm::Ptr<pat::Electron>( electronHandle, ielectron - electronBegin ) ) );
	    }
	  }	  
	}
      }
      
      int nMuons = -1;
      
      for ( std::vector<pat::Muon>::const_iterator muonBegin = muonHandle->begin(),
              muonEnd = muonHandle->end(), imuon = muonBegin;
            imuon != muonEnd; ++imuon ) {
        nMuons++;
        if ( !imuon->isGlobalMuon() ) continue;



        // Tight cuts
        bool passTight = muonIdTight_(*imuon,event) && imuon->isTrackerMuon() ;
        if (  imuon->pt() > muPtMin_ && fabs(imuon->eta()) < muEtaMax_ && 
              passTight) {

          // dereference the iterartor
          // and get the address of the object it points to
          //const pat::Muon *tempMuonObj = &(*imuon);

          //cout <<"===================== New event ==========================" <<endl;
          
          //cout << "Pushing on a muon object with pt = " << tempMuonObj->pt() << endl;

          //cout << "Muon object has index nMuons = " << nMuons << endl;
          //cout << "Calculating offset imuon-muonBegin = " << (unsigned) (imuon-muonBegin) << endl;
          
          edm::Ptr<pat::Muon> testMuonPtr(muonHandle, nMuons);     
                                      
          //cout << "Trying to test that muon object to see if we can use it's ptr" << endl;
      
          //cout << "Pt object ptr = " << testMuonPtr->pt() << endl;
           

        
          temporaryMuons_.push_back(testMuonPtr);
          
        } else {
          // Loose cuts
          if ( imuon->pt() > muPtMinLoose_ && fabs(imuon->eta()) < muEtaMaxLoose_ && 
               muonIdLoose_(*imuon,event) ) {
            looseMuons_.push_back( reco::ShallowClonePtrCandidate( edm::Ptr<pat::Muon>( muonHandle, imuon - muonBegin ) ) );
          }
        }

        
      }// end for each muon
      
      met_ = reco::ShallowClonePtrCandidate( edm::Ptr<pat::MET>( metHandle, 0),
                                             metHandle->at(0).charge(),
                                             metHandle->at(0).p4() );

      double deltaPx = 0.0;
      double deltaPy = 0.0;
      double deltaSumEt = 0.0;

      event.getByLabel (jetTag_, jetHandle);
      pat::strbitset ret1 = jetIdLoose_.getBitTemplate();
      pat::strbitset ret2 = pfjetIdLoose_.getBitTemplate();


      bool jetDebug = false;

      if (jetDebug) cout << "====================NEW LOOP OVER JETS=====================" << endl;
      
      for ( std::vector<pat::Jet>::const_iterator jetBegin = jetHandle->begin(),jetEnd = jetHandle->end(), ijet = jetBegin;
            ijet != jetEnd; ++ijet ) {

        double ptScale = 1.0; // start at 1.0, we multiply in JES later
        //Now let's evaluate JER
        if (jerFactor_ > 0.) {
          const reco::GenJet *myGenJet = ijet->genJet();
          if (myGenJet) { //Make sure we actually have one of these
            if (myGenJet->pt() > 15) { //That's what the Twiki says...
              double deltaPt = (ijet->pt() - myGenJet->pt()) * jerFactor_;
              ptScale *= max(0.0, (ijet->pt() + deltaPt) / ijet->pt());  
            }
            
          }
        }

        // Decide how you want to scale the JES
        // none = if you're not doing it the fancy way, do it the plain way
        // up = add the uncertainty to 1.0
        // down = subtract the uncertainty from 1.0
        if (fancyJES_ == "none") {
          ptScale *= jetScale_;          
        } else {
          jecUnc_->setJetEta(ijet->eta());
          jecUnc_->setJetPt(ijet->pt());
          
          double uncert = jecUnc_->getUncertainty(true);
	  if (jetDebug) cout << "  *** TOTAL = " << uncert << endl << endl << endl;
          
          if (fancyJES_ == "up") {
            ptScale *= (1 + uncert);          
          } else if (fancyJES_ == "down") {
            ptScale *= (1 - uncert);
          }
        }

	pat::Jet iRawJet = ijet->correctedJet("Uncorrected");
	deltaPx += ijet->px() * (1-ptScale);
	deltaPy += ijet->py() * (1-ptScale);
	deltaSumEt += ijet->et() *(1-ptScale);      

        reco::ShallowClonePtrCandidate scaledJet ( reco::ShallowClonePtrCandidate( edm::Ptr<pat::Jet>( jetHandle, ijet - jetBegin ),
                                                                                   ijet->charge(),
                                                                                   ijet->p4() * ptScale ) );


    
        bool passJetID = false;
        if ( ijet->isCaloJet() || ijet->isJPTJet() ) passJetID = jetIdLoose_(*ijet, ret1);
        else passJetID = pfjetIdLoose_(*ijet, ret2);
	if ( scaledJet.pt() > jetPtMin_ && fabs(scaledJet.eta()) < jetEtaMax_ && passJetID ) {
          selectedJets_.push_back( scaledJet );

          if ( muPlusJets_ ) {
        
            //Remove some jets
            bool indeltaR = false;
	    
            for( std::vector< edm::Ptr<pat::Muon> >::const_iterator muonBegin = temporaryMuons_.begin(),
                   muonEnd = temporaryMuons_.end(), imuon = muonBegin;
                 imuon != muonEnd; ++imuon ) {
              
              //const pat::Muon * iMuonPtr = (*imuon);
              if( reco::deltaR( (*imuon)->eta(), (*imuon)->phi(), scaledJet.eta(), scaledJet.phi() ) < muJetDRJets_ )
                {  indeltaR = true; }
            }

            if( !indeltaR ) {
              cleanedJets_.push_back( scaledJet );
            }// end if jet is not within dR of a muon
          }// end if mu+jets
          else {
            //Remove some jets
            bool indeltaR = false;
            for( std::vector<reco::ShallowClonePtrCandidate>::const_iterator electronBegin = selectedElectrons_.begin(),
                   electronEnd = selectedElectrons_.end(), ielectron = electronBegin;
                 ielectron != electronEnd; ++ielectron ) {
              if( reco::deltaR( ielectron->eta(), ielectron->phi(), scaledJet.eta(), scaledJet.phi() ) < eleJetDR_ )
                {  indeltaR = true; }
            }
            if( !indeltaR ) {
              cleanedJets_.push_back( scaledJet );
            }// end if jet is not within dR of an electron
          }// end if e+jets
        }// end if pass id and kin cuts
      }// end loop over jets


      double corrPx = met_.px() + deltaPx;
      double corrPy = met_.py() + deltaPy;
      reco::Particle::LorentzVector corrMetLV (corrPx,
					       corrPy,
					       0,
					       sqrt(corrPx*corrPx + corrPy*corrPy)
					       );
      
      met_ = reco::ShallowClonePtrCandidate( edm::Ptr<pat::MET>( metHandle, 0),
					     metHandle->at(0).charge(),
					     corrMetLV);

      for( std::vector< edm::Ptr<pat::Muon> >::const_iterator muonBegin = temporaryMuons_.begin(),
             muonEnd = temporaryMuons_.end(), imuon = muonBegin;
           imuon != muonEnd; ++imuon ) {

        //const pat::Muon * iMuonPtr = (*imuon);

        //Now, check that the muon isn't within muJetDRMuon_ of any jet
        bool inDeltaR_final = false;
        for (std::vector<reco::ShallowClonePtrCandidate>::const_iterator iJet = cleanedJets_.begin();
             iJet != cleanedJets_.end(); ++iJet) {
          if ( reco::deltaR((*imuon)->eta(), (*imuon)->phi(), iJet->eta(), iJet->phi()) < muJetDRMuon_ ) inDeltaR_final = true;
        }
        if (  !inDeltaR_final  ){

          //cout << " You have successfully gotten the originalObjectRef from an object with Pt = " << (*imuon)->pt() <<  endl;
          
          
          selectedMuons_.push_back( reco::ShallowClonePtrCandidate((*imuon)) );
        } 
        else {
          
          // Loose cuts
          // **imuon is a double dereference
          // de-reference the iterator
          // de-reference the Ptr
          // end result is a pat muon
          
          if ( (*imuon)->pt() > muPtMinLoose_ && fabs((*imuon)->eta()) < muEtaMaxLoose_ &&
               muonIdLoose_( (**imuon), event) ) {
            // put me back in soon
            looseMuons_.push_back( reco::ShallowClonePtrCandidate((*imuon)) );
          }
        }
      }
      

      int nleptons = 0;
      bool passlepton = false;
      if ( muPlusJets_ )
        nleptons += selectedMuons_.size();
      
      if ( ePlusJets_ ) 
        nleptons += selectedElectrons_.size();

      if ( ignoreCut(lep1Index_) || 
           ( nleptons > 0 ) ){
        passCut( ret, lep1Index_);

        if ( ignoreCut(lep2Index_) || 
             ( nleptons == 1 ) ){
          passCut( ret, lep2Index_);

          bool oneMuon = 
            ( selectedMuons_.size() == 1 && 
              looseMuons_.size() + selectedElectrons_.size() + looseElectrons_.size() == 0 
              );
          bool oneElectron = 
            ( selectedElectrons_.size() == 1 &&
              selectedMuons_.size() == 0 
              );

          bool oneMuonMuVeto = 
            ( selectedMuons_.size() == 1 && 
              looseMuons_.size() == 0 
              );
	  if(isHWW_){
	    if ( ignoreCut(lep3Index_) || 
		 (ePlusJets_ && selectedMuons_.size() == 0 && looseMuons_.size() == 0) ||
		 (muPlusJets_ && oneMuonMuVeto)
		 ) {
	      passCut(ret, lep3Index_);
	      
	      if ( ignoreCut(lep4Index_) || 
		   ( (muPlusJets_ &&  selectedElectrons_.size() == 0 && looseElectrons_.size() == 0) ^ 
		     (ePlusJets_ && selectedElectrons_.size() == 1 && looseElectrons_.size() == 0 )  )
		   //( (muPlusJets_ && oneMuon) ^ (ePlusJets_ && oneElectron )  )
		   
		   ) {
		passCut(ret, lep4Index_);   
		passlepton = true;
	      }	      
	    }
	  }
	  else {
	    if ( ignoreCut(lep3Index_) || 
		 ePlusJets_ ||
		 (muPlusJets_ && oneMuonMuVeto)
		 ) {
	      passCut(ret, lep3Index_);
	      
	      if ( ignoreCut(lep4Index_) || 
		   ( (muPlusJets_ && oneMuon) ^ (ePlusJets_ && oneElectron )  )
		   ) {
		passCut(ret, lep4Index_);	  
		passlepton = true;
	      }
	    }
	  }
	  
	  if(passlepton){
	    bool metCut = met_.pt() > metMin_;
	    if ( ignoreCut(metIndex_) ||
		 metCut ) {
	      passCut( ret, metIndex_ );
	      

	      bool zVeto = true;
	      if ( selectedMuons_.size() == 2 ) {
	      }
	      if ( selectedElectrons_.size() == 2 ) {
	      }
	      if ( ignoreCut(zvetoIndex_) ||
		   zVeto ){
		passCut(ret, zvetoIndex_);
		
		
		bool conversionVeto = true;
		if ( ignoreCut(conversionIndex_) ||
		     conversionVeto ) {
		  passCut(ret,conversionIndex_);
		  
		  
		  
		  bool cosmicVeto = true;
		  if ( ignoreCut(cosmicIndex_) ||
		       cosmicVeto ) {
		    passCut(ret,cosmicIndex_);
		    
		    if ( ignoreCut(jet1Index_) ||
			 static_cast<int>(cleanedJets_.size()) >=  1 ){
		      passCut(ret,jet1Index_);  
		    } // end if >=1 tight jets
		    
		    if ( ignoreCut(jet2Index_) ||
			 static_cast<int>(cleanedJets_.size()) >=  2 ){
		      passCut(ret,jet2Index_);  
		    } // end if >=2 tight jets
		    
		    if ( ignoreCut(jet3Index_) ||
			 static_cast<int>(cleanedJets_.size()) >=  3 ){
		      passCut(ret,jet3Index_);  
		    } // end if >=3 tight jets
		    
		    if ( ignoreCut(jet4Index_) ||
			 static_cast<int>(cleanedJets_.size()) >=  4 ){
		      passCut(ret,jet4Index_);  
		    } // end if >=4 tight jets
		    
		    if ( ignoreCut(jet5Index_) ||
			 static_cast<int>(cleanedJets_.size()) >=  5 ){
		      passCut(ret,jet5Index_);  
		    } // end if >=5 tight jets
		    
		    
		    
		  } // end if cosmic veto
		  
		} // end if conversion veto
		
	      } // end if z veto
	      
	    } // end if met cut
	    
	  } // end if == 1 lepton
	  
        } // end if == 1 tight lepton
	
      } // end if >= 1 lepton

    } // end if PV
    
  } // end if trigger


  setIgnored(ret);

  return (bool)ret;
}


double WPlusJetsEventSelector::getPtEtaJESUncert ( pat::Jet  anyJet ) {

  if (jecUnc_) {
    jecUnc_->setJetEta(anyJet.eta());
    jecUnc_->setJetPt(anyJet.pt());
    double uncert = jecUnc_->getUncertainty(true);


    double softwareUncert = 0.015;
    double pileUpUncert = (0.75 * 0.8 * 2.2) / anyJet.pt();
    // method defined with british spelling only
    int jetFlavor = anyJet.partonFlavour();
    double bjetUncert = 0.0;


    if ( abs(jetFlavor) == 5 ) {

      if ( anyJet.pt() > 50 && anyJet.pt() < 200 && fabs(anyJet.eta()) < 2.0) {
        bjetUncert = 0.02;
      } else {
        bjetUncert = 0.03;
      }
    }
          

    if (flatAdditionalUncer_ > 0) {
      uncert = sqrt(uncert*uncert
                    + flatAdditionalUncer_*flatAdditionalUncer_
                    + softwareUncert*softwareUncert
                    + pileUpUncert*pileUpUncert
                    + bjetUncert*bjetUncert);
    } else {

      uncert = sqrt (uncert*uncert
                     + softwareUncert*softwareUncert
                     + pileUpUncert*pileUpUncert
                     + bjetUncert*bjetUncert);
    }

    
    
    return uncert;
  } else {
    return -1.;
  }
  
}
