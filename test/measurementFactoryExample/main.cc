#include <ConstField.h>
#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFitterRefTrack.h>
#include <StateOnPlane.h>
#include <Track.h>
#include <TrackCand.h>
#include <TrackPoint.h>

#include <MeasurementProducer.h>
#include <MeasurementFactory.h>

#include "mySpacepointDetectorHit.h"
#include "mySpacepointMeasurement.h"

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

#include <TFile.h>
#include <TTree.h>
#include "TDatabasePDG.h"
#include <TMath.h>


//require one argument which is the number of measurements
int main(int argc, char* argv[]) {

  gRandom->SetSeed(14);

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

  double resolution = 0.02;
  const double momSmear = 3. /180.*TMath::Pi();     // rad
  const double momMagSmear = 0.1;   // relative
  // smeared start values (would come from the pattern recognition)
  const bool smearPosMom = true;   // init the Reps with smeared pos and mom
  const double posSmear = 0.1;     // cm

  // init MeasurementCreator
  genfit::MeasurementCreator measurementCreator;


  // init geometry and mag. field
  new TGeoManager("Geometry", "Geane geometry");
//  TGeoManager::Import("genfitGeom.root");
  TGeoManager::Import("analysis.root");
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0., 10.)); // 15 kGauss
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());


  // init event display
  genfit::EventDisplay* display = genfit::EventDisplay::getInstance();


  // init fitter
  genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();


  TClonesArray myDetectorHitArray("genfit::mySpacepointDetectorHit");

  // init the factory
  int myDetId(1);
  genfit::MeasurementFactory<genfit::AbsMeasurement> factory;
  //genfit::MeasurementProducer<genfit::mySpacepointDetectorHit, genfit::mySpacepointMeasurement> myProducer(&myDetectorHitArray);
  genfit::MeasurementProducer<genfit::mySpacepointDetectorHit, genfit::mySpacepointMeasurement> myProducer(&myDetectorHitArray);
  factory.addProducer(myDetId, &myProducer);


  // main loop
  for (unsigned int iEvent=0; iEvent<1; ++iEvent){

    myDetectorHitArray.Clear();

    //TrackCand
    genfit::TrackCand myCand;

    // true start values
    TVector3 pos(0, 0, 0);
    TVector3 mom(1.,0,0);
    cdc->GetEntry(0);
    pos.SetX(gRandom->Gaus(fX/10, 0.3));
    pos.SetY(gRandom->Gaus(fY/10, resolution));
    pos.SetZ(gRandom->Gaus(fZ/10, resolution));
    mom.SetPhi(gRandom->Gaus(fPhi, momSmear));
    mom.SetTheta(gRandom->Gaus(fTheta, momSmear));
    mom.SetMag(gRandom->Gaus(fMag, momMagSmear*mom.Mag() ) ); //relative error = sigma/momentum

    // helix track model
    const int pdg = 11;               // particle pdg code

//    unsigned int nMeasurements = gRandom->Uniform(5, 15);
    unsigned int nMeasurements = atoi(argv[1]);

    // covariance
//    double resolution = 0.02;
    TMatrixDSym cov(3);
    for (int i = 0; i < 3; ++i)
      cov(i,i) = resolution*resolution;

    std::vector<genfit::eMeasurementType> measurementTypes;
    for (unsigned int i=0; i<nMeasurements; ++i) {
      cdc->GetEntry(i);
      TVector3 currentPos;
      currentPos.SetX(gRandom->Gaus(fX/10, 0.3));
      currentPos.SetY(gRandom->Gaus(fY/10, resolution));
      currentPos.SetZ(gRandom->Gaus(fZ/10, resolution));

      // Fill the TClonesArray and the TrackCand
      // In a real experiment, you detector code would deliver mySpacepointDetectorHits and fill the TClonesArray.
      // The patternRecognition would create the TrackCand.
      new(myDetectorHitArray[i]) genfit::mySpacepointDetectorHit(currentPos, cov);
      myCand.addHit(myDetId, i);
    }



    TVector3 posM(pos);
    TVector3 momM(mom);
    if (smearPosMom) {
      posM.SetX(gRandom->Gaus(posM.X(),posSmear));
      posM.SetY(gRandom->Gaus(posM.Y(),posSmear));
      posM.SetZ(gRandom->Gaus(posM.Z(),posSmear));

      momM.SetPhi(gRandom->Gaus(mom.Phi(),momSmear));
      momM.SetTheta(gRandom->Gaus(mom.Theta(),momSmear));
      momM.SetMag(gRandom->Gaus(mom.Mag(), momMagSmear*mom.Mag()));
    }

    // initial guess for cov
    TMatrixDSym covSeed(6);
    for (int i = 0; i < 3; ++i)
      covSeed(i,i) = resolution*resolution;
    for (int i = 3; i < 6; ++i)
      covSeed(i,i) = pow(resolution / nMeasurements / sqrt(3), 2);


    // set start values and pdg to cand
    myCand.setPosMomSeedAndPdgCode(posM, momM, pdg);
    myCand.setCovSeed(covSeed);


    // create track
    genfit::Track fitTrack(myCand, factory, new genfit::RKTrackRep(pdg));


    // do the fit
    try{
      fitter->processTrack(&fitTrack);
    }
    catch(genfit::Exception& e){
      std::cerr << e.what();
      std::cerr << "Exception, next track" << std::endl;
      continue;
    }

    //check
    assert(fitTrack.checkConsistency());


    if (iEvent < 1000) {
      // add track to event display
      display->addEvent(&fitTrack);
    }


  }// end loop over events

  delete fitter;

  // open event display
  display->open();

}


