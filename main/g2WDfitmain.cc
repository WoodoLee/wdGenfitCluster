#include "TDatabasePDG.h"
#include "TPie.h"
#include "TPieSlice.h"
#include "TLegend.h"
#include "TFile.h"
#include "TH1.h"
#include "TRandom3.h"
#include "TNtuple.h"
#include "TSystem.h"
#include "TEfficiency.h"
#include "TTree.h"
#include <AbsFinitePlane.h>
#include <AbsFitterInfo.h>
#include <AbsKalmanFitter.h>
#include <AbsMeasurement.h>
#include <AbsTrackRep.h>
#include <ConstField.h>
#include <DAF.h>
#include <DetPlane.h>
#include <EventDisplay.h>
#include <Exception.h>
#include <FieldManager.h>
#include <Fit/Fitter.h>
#include <FullMeasurement.h>
#include <GFGbl.h>
#include <HelixTrackModel.h>
#include <KalmanFitStatus.h>
#include <KalmanFittedStateOnPlane.h>
#include <KalmanFitter.h>
#include <KalmanFitterInfo.h>
#include <KalmanFitterRefTrack.h>
#include <MaterialEffects.h>
#include <Math/Functor.h>
#include <Math/Vector3D.h>
#include <MeasuredStateOnPlane.h>
#include <MeasurementCreator.h>
#include <MeasurementOnPlane.h>
#include <PlanarMeasurement.h>
#include <ProlateSpacepointMeasurement.h>
#include <RKTrackRep.h>
#include <RectangularFinitePlane.h>
#include <ReferenceStateOnPlane.h>
#include <SharedPlanePtr.h>
#include <SpacepointMeasurement.h>
#include <StateOnPlane.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TChain.h>
#include <TDatabasePDG.h>
#include <TEllipse.h>
#include <TEveManager.h>
#include <TF1.h>
#include <TF2.h>
#include <TFrame.h>
#include <TGeoManager.h>
#include <TGeoMaterialInterface.h>
#include <TGraph2D.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TH2I.h>
#include <TH3F.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMarker.h>
#include <TMath.h>
#include <TPad.h>
#include <TPie.h>
#include <TPolyLine.h>
#include <TPolyLine3D.h>
#include <TPolyMarker.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TRandom2.h>
#include <TStyle.h>
#include <TVector3.h>
#include <TVirtualPad.h>
#include <Tools.h>
#include <Track.h>
#include <TrackCand.h>
#include <TrackCandHit.h>
#include <TrackPoint.h>
#include <WireMeasurement.h>
#include <WirePointMeasurement.h>
#include <cassert>
#include <iostream>
#include <memory>
#include <vector>
#include <fstream>
#include "TAxis.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
using namespace ROOT::Math;
using namespace std;
using namespace genfit;
using namespace RooFit;
int main(int argc, char* argv[] )
{

  //Conditions
  double BZ = 30.;
  const int pdgWanted  = -11;  // particle pdg code -11 = positron
  //double pitch = 0.02;
  double pitch = 0.019;
  double resolution =  pitch ;
  gRandom->SetSeed(14);
  // init MeasurementCreator
  genfit::MeasurementCreator measurementCreator;
  // init geometry and mag. field
  TChain *Thit  = new TChain("Hit");
  TChain *Tfit  = new TChain("Fit");
  Thit -> Add(argv[1]);
  Tfit -> Add(argv[1]);
  //TString fileGeom = "/raid01/wdlee/data/g4data/sysError/dp/tdrValue/g2wdGeom.gdml";
  TString fileGeom = argv[3];
  new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import(fileGeom.Data());
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0., BZ)); // BZ kGauss
  genfit::MaterialEffects *mateff = genfit::MaterialEffects::getInstance();
  mateff->setEnergyLossBrems(true);
  mateff->setNoiseBrems(true);
  mateff->setEnergyLossBetheBloch(true);
  mateff->setNoiseBetheBloch(true);
  mateff->setNoiseCoulomb(true);
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
  // init event display
  //genfit::EventDisplay* display = genfit::EventDisplay::getInstance();
  // init fitter
  //genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();
  //genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitter();
  genfit::AbsKalmanFitter* fitter = new genfit::DAF();
  //genfit::AbsKalmanFitter* fitter = new genfit::DAF(true, 1e-3, 1e-3);
  //const genfit::eFitterType fitter = genfit::DafRef;
  //const genfit::eFitterType fitterId = genfit::DafSimple;

  //GenFit options
  const genfit::eMultipleMeasurementHandling mmHandling = genfit::weightedAverage;
  //const genfit::eMultipleMeasurementHandling mmHandling = genfit::unweightedClosestToReference;
  //const genfit::eMultipleMeasurementHandling mmHandling = genfit::unweightedClosestToPrediction;
  //const genfit::eMultipleMeasurementHandling mmHandling = genfit::weightedClosestToReference;
  //const genfit::eMultipleMeasurementHandling mmHandling = genfit::weightedClosestToPrediction;
  //const genfit::eMultipleMeasurementHandling mmHandling = genfit::unweightedClosestToReferenceWire;
  //const genfit::eMultipleMeasurementHandling mmHandling = genfit::unweightedClosestToPredictionWire;
  //const genfit::eMultipleMeasurementHandling mmHandling = genfit::weightedClosestToReferenceWire;
  //const genfit::eMultipleMeasurementHandling mmHandling = genfit::weightedClosestToPredictionWire;
  fitter->setMultipleMeasurementHandling(mmHandling);
  double dPVal(1.E-3);
  double  dRelChi2(0.2);
  //double  dChi2Ref(1.);
  double  nMinIter(3);
  double  nMaxIter(100);
  int nMaxFailed(-1);
  //bool refit(false);
  //bool  squareRootFormalism(false);

  fitter->setMinIterations(nMinIter);
  fitter->setMaxIterations(nMaxIter);
  fitter->setRelChi2Change(dRelChi2);
  fitter->setMaxFailedHits(nMaxFailed);

  TFile* f = new TFile(argv[2],"recreate");
  TTree *tr= new TTree("data", "data");
  
  double momInit, momFitted , momDiff, momRes;
  double firstHitX, firstHitY, firstHitZ, firstHitT; 
  tr -> Branch ("momInit",  &momInit , "momInit/D");
  tr -> Branch ("momFitted" ,  &momFitted  , "momFittied/D");
  tr -> Branch ("momDiff" ,  &momDiff  , "momDiff/D");
  tr -> Branch ("momRes" ,  &momRes  , "momRes/D");
  tr -> Branch ("firstHitX" ,  &firstHitX , "firstHitX/D");
  tr -> Branch ("firstHitY" ,  &firstHitY , "firstHitY/D");
  tr -> Branch ("firstHitZ" ,  &firstHitZ , "firstHitZ/D");
  tr -> Branch ("firstHitT" ,  &firstHitT , "firstHitT/D");
  
  //TFile* file = new TFile(fileData.Data());
  //std::cout << "\033[1;32m [Notice] File " << fileData.Data() << " is open.\033[0m" << std::endl;
  
  //TTree *treeHit = (TTree*)(file->Get("Hit"));
  int chitEventID;
  double chitTime;
  double chitPosX,chitPosY, chitPosZ;
  double cMomX, cMomY, cMomZ;


  
  Thit -> SetBranchAddress("chitEventID",&chitEventID);
  Thit -> SetBranchAddress("chitTime",&chitTime);
  Thit -> SetBranchAddress("chitPosX",&chitPosX);
  Thit -> SetBranchAddress("chitPosY",&chitPosY);
  Thit -> SetBranchAddress("chitPosZ",&chitPosZ);

  Tfit -> SetBranchAddress("cMomX",&cMomX);
  Tfit -> SetBranchAddress("cMomY",&cMomY);
  Tfit -> SetBranchAddress("cMomZ",&cMomZ);
  

  vector < vector <double> > hitT;
  vector < vector <double> > hitX;
  vector < vector <double> > hitY;
  vector < vector <double> > hitZ;

  hitT.clear();
  hitX.clear();
  hitY.clear();
  hitZ.clear();

  vector <double> hitTTemp;
  vector <double> hitXTemp;
  vector <double> hitYTemp;
  vector <double> hitZTemp;

  int preEid ;
  int postEid;
  //int eIDhit = treeHit -> GetEntries();
  int eIDhit = Thit -> GetEntries();
  std::cout << "\033[1;32m [Notice] " << eIDhit << " number of Entries. \033[0m" << std::endl;
  
  for (int i = 0; i < eIDhit; i++)
  {
    //std::cout << "\033[1;32m [Notice] event [ " << i << " ] is being fitted \033[0m" << std::endl;
    //treeHit -> GetEntry(i);
    Thit -> GetEntry(i);
    preEid = chitEventID;
    hitTTemp.push_back(chitTime);
    hitXTemp.push_back(chitPosX);
    hitYTemp.push_back(chitPosY);
    hitZTemp.push_back(chitPosZ);

    Thit -> GetEntry(i+1);
    postEid = chitEventID;
    //cout << "preEid = " << preEid << endl;
    //cout << "postEid = " << postEid << endl;
    //cout << "hitTime = " << hitTime << endl;
    if(preEid == eIDhit)
    {
      hitT.push_back(hitTTemp);
      hitX.push_back(hitXTemp);
      hitY.push_back(hitYTemp);
      hitZ.push_back(hitZTemp);

      hitTTemp.clear();
      hitXTemp.clear();
      hitYTemp.clear();
      hitZTemp.clear();

    }
    if(preEid != postEid)
    {
      hitT.push_back(hitTTemp);
      hitX.push_back(hitXTemp);
      hitY.push_back(hitYTemp);
      hitZ.push_back(hitZTemp);

      hitTTemp.clear();
      hitXTemp.clear();
      hitYTemp.clear();
      hitZTemp.clear();

    }
  }
  int trackN = hitX.size();
  
  // Looping over entries
  std::cout << "\033[1;31m [Notice] " << trackN << " number of Tracks. \033[0m" << std::endl;
  for ( int i = 0; i < trackN; i++ ) 
  {
    if(i % 1000 == 0)
    {
      std::cout << "\033[1;32m [Notice] track [ " << i << " ] is being fitted. \033[0m" << std::endl;
    }

    Tfit -> GetEntry(i);

    //Array for momentums
    //chekc the real Hit number ;
    //set the initial state(for Kalman filter step 0) Units: cm / GeV/c
    //std::cout << "\033[1;35m [Notice] Initial Condition [ " << i << " ] is being fitted. \033[0m" << std::endl;
    TVector3 pos( hitX[i][0] / 10.  , hitY[i][0] / 10.    , hitZ[i][0] / 10.  );
    TVector3 mom( cMomX / 1000. , cMomY / 1000.  , cMomZ / 1000. );
    
    firstHitX = hitX[i][0];
    firstHitY = hitY[i][0];
    firstHitZ = hitZ[i][0];
    firstHitT = hitT[i][0];
    //pos.Print();
    //mom.Print();
    //Calculate the Magnitude of Initial momentum - will be use in Momentum cutting
    double momMag = mom.Mag();
    //mom.Print();
    const int pdg = pdgWanted ;     // particle pdg
    //const double charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge()/(3.);
    //genfit::HelixTrackModel* helix = new genfit::HelixTrackModel(pos, mom, charge);
    //measurementCreator.setTrackModel(helix);
    //get the charge : why devide by 3? - Curious(from Dasilva) 
    //const double charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge()/(3.);
    //std::cout << "charge  = " << charge << std::endl;
    // smeared start values
    //const bool smearPosMom = true;     // init the Reps with smeared pos and mom
    //const bool smearPosMom = false;     // init the Reps with smeared pos and mom
    //const double posSmear = 0.019;     // cm
    //const double momSmear = 3. /180.*TMath::Pi();     // rad
    //const double momSmear = 0.019;     // rad
    //const double momMagSmear = 0.019;   // relative
    //hmomshit -> Fill( momMag* 1000. , countHits);
    TVector3 posM(pos);
    TVector3 momM(mom);
    //if (smearPosMom)
    //{
    //  posM.SetX(gRandom->Gaus(posM.X(),posSmear));
    //  posM.SetY(gRandom->Gaus(posM.Y(),posSmear));
    //  posM.SetZ(gRandom->Gaus(posM.Z(),posSmear));
    //  momM.SetX(gRandom->Gaus(momM.X(),momSmear));
    //  momM.SetY(gRandom->Gaus(momM.Y(),momSmear));
    //  momM.SetZ(gRandom->Gaus(momM.Z(),momSmear));
    //}
    // approximate covariance(Kalman filter step.1)

    int nHits = hitX[i].size();
    //if (nHits ==0) continue;

    TMatrixDSym covM(6);
    for (int k = 0; k < 3; ++k)
    {
      covM(k,k) = resolution * resolution;
      //covM(k,k) = pow(posSmear,2);
    }
    for (int k = 3; k < 6; ++k)
    {
      covM(k,k) = pow(resolution / nHits / sqrt(3), 2);
    }
    // trackrep
    genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);
    // smeared start state
    genfit::MeasuredStateOnPlane stateSmeared(rep);
    stateSmeared.setPosMomCov(posM, momM, covM);
    genfit::MeasuredStateOnPlane stateRef(rep);
    stateRef.setPosMomCov(pos, mom, covM);
    // remember original initial state
    const genfit::StateOnPlane stateRefOrig(stateRef);
    // create track
    TVectorD seedState(6);
    TMatrixDSym seedCov(6);
    stateSmeared.get6DStateCov(seedState, seedCov);
    genfit::Track fitTrack(rep, seedState, seedCov);
    // create random measurement types
    TMatrixDSym hitCov(2);
    hitCov.UnitMatrix();
    hitCov *=resolution *resolution;
    for (int j = 0; j < nHits ; j++)
    {
      //pick the hit positions
      TVector3 cPos(hitX[i][j]/ 10. , hitY[i][j] / 10.  , hitZ[i][j] / 10. );
      //double cPer = cPos.Perp() ; 
      if(cPos == TVector3(0.,0.,0.)) continue;
      double phi = cPos.Phi();
      //double phi = cPos.Theta();
      TVector3 o(18.425 * TMath::Cos(phi), 18.425 * TMath::Sin(phi), 0.);
      TVector3 u(0., 0., 1. );
      TVector3 v = o.Unit();
      TVectorD hitCoords(2);
      hitCoords[0] = cPos.Z();
      hitCoords[1] = cPos.Perp() - 18.425;
      //std::cout << "\033[1;32m [Notice] cPos [ " << j << " ]  \033[0m" << std::endl;
      //cPos.Print();
      genfit::PlanarMeasurement * measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, 0., j ,nullptr );
      measurement ->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(o,u,v) ), 0 );
      fitTrack.insertPoint(new genfit::TrackPoint(measurement,&fitTrack ));
    }
    
    // assert(fitTrack.checkConsistency());
    // do the fit
    fitTrack.checkConsistency();
    fitter->processTrack(&fitTrack);
    //display->addEvent(&fitTrack);
    if (i % 1000 == 0)
    {
      std::cout << "\033[1;32m [Notice] Fitting track [ " << i << " ] is done. \033[0m" << std::endl;
    }

    genfit::TrackPoint* tp = fitTrack.getPointWithMeasurementAndFitterInfo(0, rep);
    if (tp == nullptr) continue;
    //genfit::KalmanFittedStateOnPlane kfsop(*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getBackwardUpdate()));
    genfit::KalmanFittedStateOnPlane kfsop(*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getForwardUpdate()));
    const TVectorD& referenceState = stateRefOrig.getState();
    const TVectorD& state = kfsop.getState();
    const TMatrixDSym covf = kfsop.getCov();
    //genfit::FitStatus *fitstat = fitTrack->getFitStatus(rep);
    //bool fitconv = fitstat->isFitConverged();
    //double chi2 = fitstat->getChi2();
    //double NDF = fitstat->getNdf();
    //fittedTrack->_chi2ndf = chi2/NDF;
    
    //if(!fitconv)
    //{
    //  cout << "Track could not be fitted successfully" << endl;
    //  delete fittedTrack;
    //  delete fitTrack;
    //  continue;
    //}
//
//    //genfit::TrackPoint* tp = fitTrack->getPointWithMeasurementAndFitterInfo(0, rep);
//    //if(tp==NULL){
//    //  cout << "Track has no TrackPoint with fitterInfo" << endl;
//    //  delete fittedTrack;
//    //  delete fitTrack;
//    //  continue;
    //}







    //genfit::FitStatus *fitstat = fitTrack.getFitStatus(rep);
    //bool fitconv = fitstat->isFitConverged();
    //if(!fitconv)
    //{
    //  //cout << "Track could not be fitted successfully" << endl;
    //  continue;
    //}

    momInit = (1./referenceState[0]) *  1000.; 
    momFitted = (1./state[0]) * 1000.; 
    momDiff = momInit - momFitted; 
    momRes = (momInit - momFitted) / momInit;
    tr -> Fill();


  }
  f  -> cd();
  tr -> Write();   
  f  -> Close();
  
  delete fitter;

  
  //display->setOptions("ABDEFHMPT"); // G show geometry
  //display->open();


}
