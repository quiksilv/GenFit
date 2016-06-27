#include <ConstField.h>
#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFitterRefTrack.h>
#include <KalmanFittedStateOnPlane.h>
#include <KalmanFitterInfo.h>
#include <StateOnPlane.h>
#include <Track.h>
#include <TrackCand.h>
#include <TrackPoint.h>

#include <MeasurementProducer.h>
#include <MeasurementFactory.h>

#include "mySpacepointDetectorHit.h"
#include "myProlateSpacepointMeasurement.h"

#include "FullMeasurement.h"

#include <MaterialEffects.h>
#include <RKTrackRep.h>
#include <TGeoMaterialInterface.h>

#include <EventDisplay.h>

#include <HelixTrackModel.h>
#include <MeasurementCreator.h>

#include <TDatabasePDG.h>
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TRandom.h>
#include <TVector3.h>
#include <vector>

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include "TDatabasePDG.h"
#include <TMath.h>


//require one argument, the number of measurements
int main(int argc, char** argv) {


//  gROOT->SetStyle("Plain");
//  gStyle->SetPalette(1);
//  gStyle->SetOptFit(1111);
  const unsigned int nMeasurements = 150;
  double resolution = 0.02;
  const double momentum = 0.1045;
  //TH1D *hmomRes = new TH1D("hmomRes","mom res",500,-20*resolution*momentum/nMeasurements,20*resolution*momentum/nMeasurements);
  TH1D *hmomRes = new TH1D("hmomRes","mom res",500,-1*resolution,1*resolution);


  TFile *ftr3 = new TFile("analysis.root");
  TTree* cdc = (TTree*)ftr3->Get("cdc");
  Float_t fX;
  Float_t fY;
  Float_t fZ;
  Float_t fPhi;
  Float_t fTheta;
  Float_t fMag;
  cdc->SetBranchAddress("x", &fX);
  cdc->SetBranchAddress("y", &fY);
  cdc->SetBranchAddress("z", &fZ);
  cdc->SetBranchAddress("phi", &fPhi);
  cdc->SetBranchAddress("theta", &fTheta);
  cdc->SetBranchAddress("mag", &fMag);

  double resolutionWire = 0.3;
  const double momSmear = 3. /180.*TMath::Pi();     // rad
  const double momMagSmear = 0.1;   // relative
  // smeared start values (would come from the pattern recognition)
  //const bool smearPosMom = true;   // init the Reps with smeared pos and mom
  //const double posSmear = 0.1;     // cm

  //total number of measurements, or rather the ending index
  //unsigned int nMeasurements = atoi(argv[1]);
  //unsigned int nMeasurements = 150;
  // init MeasurementCreator
  genfit::MeasurementCreator measurementCreator;


  // init geometry and mag. field
  new TGeoManager("Geometry", "Geane geometry");
//  TGeoManager::Import("genfitGeom.root");
  TGeoManager::Import("analysis.root");
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(10.,0., 0.)); // 15 kGauss
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
  genfit::MaterialEffects::getInstance()->setEnergyLossBrems(false);
//  genfit::MaterialEffects::getInstance()->setNoEffects();

  // init event display
  genfit::EventDisplay* display = genfit::EventDisplay::getInstance();

  // init fitter
  genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack(10);
  fitter->setMinIterations(3);
//  fitter->setRelChi2Change(0.2);
//  fitter->setDeltaChi2Ref(1);
  fitter->setDebugLvl(1);

  TClonesArray myDetectorHitArray("genfit::mySpacepointDetectorHit");

  // init the factory
  int myDetId(1);
  genfit::MeasurementFactory<genfit::AbsMeasurement> factory;
  genfit::MeasurementProducer<genfit::mySpacepointDetectorHit, genfit::myProlateSpacepointMeasurement> myProducer(&myDetectorHitArray);
  factory.addProducer(myDetId, &myProducer);

  // main loop
  for (unsigned int iEvent=0; iEvent<1; ++iEvent){
    gRandom->SetSeed(15);

    myDetectorHitArray.Clear();

    // true start values
    TVector3 pos(0, 0, 0);
    TVector3 mom(1.,0.,0.);
    cdc->GetEntry(0);
    pos.SetX(gRandom->Gaus(fX/10, resolutionWire));
    pos.SetY(gRandom->Gaus(fY/10, resolution));
    pos.SetZ(gRandom->Gaus(fZ/10, resolution));
    mom.SetPhi(gRandom->Gaus(fPhi, momSmear));
    mom.SetTheta(gRandom->Gaus(fTheta, momSmear));
    mom.SetMag(gRandom->Gaus(fMag, momMagSmear*mom.Mag() ) ); //relative error = sigma/momentum

    TVector3 myDir(mom);
    myDir.SetMag(1);
    TVector3 posInit = pos - 0.5*myDir;
    const int pdg = 11;               // particle pdg code

    // covariance
    TMatrixDSym cov(3);
    for (int i = 0; i < 3; ++i)
      cov(0,0) = resolutionWire*resolutionWire;
      cov(1,1) = resolution*resolution;
      cov(2,2) = resolution*resolution;
    // initial guess for cov
    TMatrixDSym covSeed(6);
    for (int i = 0; i < 3; ++i)
      covSeed(i,i) = resolution*resolution;
    for (int i = 3; i < 6; ++i)
      covSeed(i,i) = pow(resolution / nMeasurements / sqrt(3), 2);

    // trackrep
    genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);
//    rep->setPropDir(1);

genfit::MeasuredStateOnPlane stateRef(rep);
rep->setPosMomCov(stateRef, pos, mom, covSeed);
// remember original initial state
const genfit::StateOnPlane stateRefOrig(stateRef);

    // smeared start state
    genfit::MeasuredStateOnPlane stateSmeared(rep);
    rep->setPosMomCov(stateSmeared, pos, mom, covSeed);
    // create track
    TVectorD seedState(6);
    TMatrixDSym seedCov(6);
    rep->get6DStateCov(stateSmeared, seedState, seedCov);
    genfit::Track fitTrack(rep, seedState, seedCov);
    genfit::Track secondTrack(rep->clone(), seedState, seedCov);

    std::vector<genfit::eMeasurementType> measurementTypes;
    for (unsigned int i=0; i<nMeasurements; ++i)
      measurementTypes.push_back(genfit::ProlateSpacepoint);
      //measurementTypes.push_back(genfit::Pixel);
    for (unsigned int i=0; i<nMeasurements; ++i) {
      cdc->GetEntry(i);
      TVector3 currentPos;
      currentPos.SetX(gRandom->Gaus(fX/10, resolutionWire));
      currentPos.SetY(gRandom->Gaus(fY/10, resolution));
      currentPos.SetZ(gRandom->Gaus(fZ/10, resolution));

      std::vector<genfit::AbsMeasurement*> measurements = measurementCreator.create(measurementTypes[i], currentPos, cov);
      genfit::TrackPoint* trackPoint = new genfit::TrackPoint(measurements, &fitTrack);
      if(i<87) {
        fitTrack.insertPoint(trackPoint);
      } else {
        secondTrack.insertPoint(trackPoint);
      }
    }

    //check
    assert(fitTrack.checkConsistency());
    assert(secondTrack.checkConsistency());

    // do the fit
    try{
      fitter->processTrack(&fitTrack);
      //fitter->processTrack(&secondTrack, false);
      fitter->processTrack(&secondTrack);
    }
    catch(genfit::Exception& e){
      std::cerr << e.what();
      std::cerr << "Exception, next track" << std::endl;
      continue;
    }
    bool fullMeasurement = false;
    if (fullMeasurement) {
      genfit::FullMeasurement* fullM = new genfit::FullMeasurement(secondTrack.getFittedState());
      fitTrack.insertPoint(new genfit::TrackPoint(fullM, &fitTrack));
    }
    else
      fitTrack.mergeTrack(&secondTrack);

    rep->setPropDir(1);

    fitter->processTrack(&fitTrack);
    if (iEvent < 1000) {
      // add track to event display
//      display->addEvent(&fitTrack);
        std::vector<genfit::Track*> event;
        event.push_back(&fitTrack);
        display->addEvent(event);
    }

      genfit::TrackPoint* tp = fitTrack.getPointWithMeasurementAndFitterInfo(0, rep);
      if (tp == NULL) {
        std::cout << "Track has no TrackPoint with fitterInfo! \n";
        continue;
      }
      genfit::KalmanFittedStateOnPlane kfsop(*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getBackwardUpdate()));
      if (true) {
        std::cout << "state before extrapolating back to reference plane \n";
        kfsop.Print();
      }

      // extrapolate back to reference plane.
      try{
        rep->extrapolateToPlane(kfsop, stateRefOrig.getPlane());;
      }
      catch(genfit::Exception& e){
        std::cerr<<"Exception, next track"<<std::endl;
        std::cerr << e.what();
        continue;
      }
      // calculate pulls
      const TVectorD& referenceState = stateRefOrig.getState();

      const TVectorD& _state = kfsop.getState();
      const TMatrixDSym& _cov = kfsop.getCov();

      double _pval = fitter->getPVal(&fitTrack, rep);
      const double charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge()/(3.);
      Double_t entry = charge/_state[0] - momentum;
      hmomRes->Fill(entry);

  }// end loop over events

  delete fitter;

  TCanvas* c1 = new TCanvas();
  hmomRes->Fit("gaus", "SMEQ");
  hmomRes->Draw();
  c1->Update();
  c1->Write();
  // open event display
  display->open();

}


