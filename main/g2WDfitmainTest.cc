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
int main()
{
  //Conditions
  double BZ = 30.;
  const int pdgWanted  = -11;  // particle pdg code -11 = positron
  //double pitch = 0.02;
  double pitch = 0.019;
  vector<int> pvaltest;
  vector<int> negM;
  vector<double> Effv;
  vector<double> momFitting;
  vector<double> momInitial;
  vector<double> Effa;
  double resolution =  pitch ;
  gRandom->SetSeed(14);
  // init MeasurementCreator
  genfit::MeasurementCreator measurementCreator;
  // init geometry and mag. field
  new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("/home/wdlee/ssd1/work/data/genfit2g2Data/Detector.root");
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0., BZ)); // BZ kGauss
  genfit::MaterialEffects *mateff = genfit::MaterialEffects::getInstance();
  mateff->setEnergyLossBrems(true);
  mateff->setNoiseBrems(true);
  mateff->setEnergyLossBetheBloch(true);
  mateff->setNoiseBetheBloch(true);
  mateff->setNoiseCoulomb(true);
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
  // init event display
  genfit::EventDisplay* display = genfit::EventDisplay::getInstance();
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

  const bool debug =  true; //false;
  //if (debug)
  //fitter->setDebugLvl(10);
  // create histograms
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptFit(1111);
  int nbin1 =150; //500
  int nbin2 =200; //200
  double mRr = 1;
  double mRrC = 0.3;
  double mDr = 300;
  double mDrC = 300;

  TH1D *hmomB = new TH1D("All","Momentum Distribution before",nbin1, 0 , 300);
  TH1D *hmomA = new TH1D("All","Momentum Distribution After",nbin1, 0 , 300);
  TH1D *hmomE = new TH1D("All","Efficiency of convergence",nbin1, 0 , 300);

  TH1D *hmomResAll = new TH1D("All","Momentum Resolution",nbin1, -mRr , mRr);
  TH1D *hmomDiffAll = new TH1D("All","Momentum Difference",nbin1, -mDr , mDr);
  double momRange = 300.;
  int window = 100;
  double windowRange = momRange / window;
  
  TH1D *hmomRes[window];
  TH1D *hmomResCut[window];
  TH1D *hmomDiff[window];
  TH1D *hmomDiffCut[window];
  TH1D *hmomSig[window];
  TString name;
  TString title;
  
  for(int i=0; i < window; i++)
  {
    double rangeMin = i *  windowRange ;
    double rangeMax = (i+1) * windowRange ;
    title . Form("%i MeV < P_{Initial} < %i + 3 MeV/c", rangeMin);
    hmomRes[i]     = new TH1D(title . Data() , "momentumswum Resolution",nbin1, -mRr , mRr);
    hmomResCut[i]  = new TH1D(title . Data() , "Momentum Resolution Cut",nbin1, -mRr , mRr);
    hmomDiff[i]    = new TH1D(title . Data() , "Momentum Difference",nbin1, -mDr , mDr);
    hmomDiffCut[i] = new TH1D(title . Data() , "Momentum Difference Cut",nbin1, -mDr , mDr);
    hmomSig[i]     = new TH1D(title . Data() , "Momentum Difference / true Momentum",nbin1, -mDr , mDr);
  }


  TH1D *hmomResO200 = new TH1D("200 MeV/c < P_{Initial} < 275 MeV/c","Momentum Resolution",nbin1*5, -mRr , mRr);
  TH1D *hmomResCutO200 = new TH1D("200MeV/c P_{Initial} < 275 MeV/c","Momentum Resolution Cut < 0.3",nbin1, -mRrC , mRrC);
  TH1D *hmomDiffO200 = new TH1D("P_{Initial} > 200 MeV/c","Momentum Difference",750, -mDr , mDr);
  TH1D *hmomDiffCutO200 = new TH1D("P_{Initial} > 200 MeV/c","Momentum Difference Cut Resolution < 0.2",150, -mDrC , mDrC);
  TH1D *hmomResO200fit = new TH1D(" P_{Initial} > 200 MeV/c","Momentum Resolution",nbin1, -1 , 1);
  TH1D *hmomResCutO200fit = new TH1D("P_{Initial} > 200 MeV/c","Momentum Resolution Cut < 0.2",nbin1, -0.5 , 0.5);
  TH1D *hmomDiffO200fit = new TH1D("P_{Initial} > 200 MeV/c","Momentum Difference",nbin1, -10 , 10);
  TH1D *hmomDiffCutO200fit = new TH1D("P_{Initial} > 200 MeV/c","Momentum Difference Cut Resolution < 0.2",nbin1, -10 , 10);
  TH1D *hmoms = new TH1D(" p_{Fitting} ","  p_{Fitting} Distribution , ALL [MeV/c]",nbin1, -300 , 300);
  TH1D *hmomref = new TH1D(" p_{Initial} "," p_{Initial} Distribution , ALL [MeV/c]",nbin1, -300 , 300);
  TH1D *hmomsCut   = new TH1D(" p_{Fitting} ","  p_{Fitting} Distribution, p_{Initial} > 200 MeV/c  [MeV/c]",nbin1, -300 , 300);
  TH1D *hmomrefCut = new TH1D(" p_{Initial} "," p_{Initial} Distribution, p_{Initial} > 200 MeV/c [MeV/c]",nbin1, -300 , 300);
  TH1D *hmomNeg       = new TH1D(" p_{Fitting} ","  p_{Fitting} Distribution, p_{Fitting} < 0 MeV/c  [MeV/c]",nbin1, -300 , 0); 
  TH1D *hmomrefNeg    = new TH1D(" p_{Initial} "," p_{Initial} Distribution, p_{Fitting} < 0 MeV/c [MeV/c]",nbin1, 0 , 300);
  TH1D *hmomdiffNeg   = new TH1D(" p_{difference} ","  p_{difference} Distribution, p_{Fitting} < 0 MeV/c  [MeV/c]",nbin1, -300 , 0); 
  TH1D *hmomLow       = new TH1D(" p_{low} ","  p_{low} Distribution, p_{low} < 150 MeV/c  [MeV/c]",nbin1, 0 , 300); 
  TH2D *hmomshit = new TH2D("ALL", "mom vs # of hits", nbin1, 0, 300, nbin1/2, 0, 100);
  TH1D *hmomcheck = new TH1D(" p val = 0 "," Mom diff. p val =0  [MeV/c]",nbin1, -150 , 150);
  //2D his 
  TH2D *hChiRes = new TH2D("All","Chi2/Ndf vs Mom resolution", nbin1, 0 , 10 , nbin1, -10 , 10);
  TH2D *hChiRes200 = new TH2D("All","Chi2/Ndf vs Mom resolution mom 200", nbin1, 0 , 10 , nbin1, -10 , 10);
  
  double Eff[window];
  double EffO200 = 0.;
  double posiNum[window];
  double posiNumAll = 0.;
  double posiNumFit[window];
  double momSig[window];

  TH1D *pVal = new TH1D("pVal","p-value",nbin2,0.,1.);
  TH1D *chI2 = new TH1D("chI2","chi sqare",nbin2,0,100);
  TH1D *Ndf = new TH1D("NDF","NDF",50,0,100);
  TH1D *chI2N = new TH1D("chi2/NDF","chi2/NDF",nbin2,0,10);

  TH1D *covFx  = new TH1D("Final Covariance X ","Final Cov",nbin1,0,0.2);
  TH1D *covFy  = new TH1D("Final Covariance Y ","Final Cov",nbin1,0,0.2);
  TH1D *covFz  = new TH1D("Final Covariance Z","Final Cov",nbin1,0,0.2);
  TH1D *covFpos  = new TH1D("Final Covariance Pos","Final Pos Cov",nbin1,0,0.2);

  TH1D *covFpx  = new TH1D("Final Covariance Px","Final Cov",nbin1,0,1);
  TH1D *covFpy  = new TH1D("Final Covariance Py","Final Cov",nbin1,0,1);
  TH1D *covFpz  = new TH1D("Final Covariance Pz","Final Cov",nbin1,0,1);
  TH1D *covFmom  = new TH1D("Final Covariance Momentum","Final Mom Cov",nbin1,0,1);

  TH1D *weights = new TH1D("weights","Daf vs true weights",500,-1.01,1.01);

  TString fileName = "/home/wdlee/ssd1/work/build/g2wdPGID_build/g2wd.root";
  TFile* file = new TFile(fileName.Data());
  TFile* f = new TFile("FittingResult.root","recreate");
  TTree *tr= new TTree("data", "data");

  double momInit, momFitted , momDiff, momRes;

  tr -> Branch ("momInit",  &momInit , "momInit/D");
  tr -> Branch ("momFitted" ,  &momFitted  , "momFittied/D");
  tr -> Branch ("momDiff" ,  &momDiff  , "momDiff/D");
  tr -> Branch ("momRes" ,  &momRes  , "momRes/D");

  std::cout << "\033[1;32m [Notice] File " << fileName.Data() << " is open.\033[0m" << std::endl;
  TTree *treeHit = (TTree*)(file->Get("Hit"));

  int hitEventID;
  double hitTime;
  double hitPosX,hitPosY, hitPosZ;
  double hitPMag;
  double hitPX, hitPY, hitPZ;
  double hitRA;
  double hitR;
  double eDep;
  double hitAngle;
  double VolID;

  treeHit -> SetBranchAddress("hitEventID",&hitEventID);
  treeHit -> SetBranchAddress("hitTime",&hitTime);
  treeHit -> SetBranchAddress("hitPosX",&hitPosX);
  treeHit -> SetBranchAddress("hitPosY",&hitPosY);
  treeHit -> SetBranchAddress("hitPosZ",&hitPosZ);
  treeHit -> SetBranchAddress("hitPX",&hitPX);
  treeHit -> SetBranchAddress("hitPY",&hitPY);
  treeHit -> SetBranchAddress("hitPZ",&hitPZ);
  treeHit -> SetBranchAddress("hitPMag",&hitPMag);
  treeHit -> SetBranchAddress("hitAngle",&hitAngle);
  treeHit -> SetBranchAddress("eDep",&eDep);
  treeHit -> SetBranchAddress("VolID",&VolID);



  vector < vector <double> > hitT;
  vector < vector <double> > hitX;
  vector < vector <double> > hitY;
  vector < vector <double> > hitZ;

  vector < vector <double> > momX;
  vector < vector <double> > momY;
  vector < vector <double> > momZ;

  hitT.clear();
  hitX.clear();
  hitY.clear();
  hitZ.clear();

  momX.clear();
  momY.clear();
  momZ.clear();


  vector <double> hitTTemp;
  vector <double> hitXTemp;
  vector <double> hitYTemp;
  vector <double> hitZTemp;
  vector <double> momXTemp;
  vector <double> momYTemp;
  vector <double> momZTemp;



  int preEid ;
  int postEid;



  int eIDhit = treeHit -> GetEntries();
  std::cout << "\033[1;32m [Notice] " << eIDhit << " number of Entries. \033[0m" << std::endl;

  for (int i = 0; i <= eIDhit; i++)
  {
    std::cout << "\033[1;32m [Notice] event [ " << i << " ] is being fitted \033[0m" << std::endl;
    treeHit -> GetEntry(i);
    preEid = hitEventID;
    hitTTemp.push_back(hitTime);
    hitXTemp.push_back(hitPosX);
    hitYTemp.push_back(hitPosY);
    hitZTemp.push_back(hitPosZ);
    momXTemp.push_back(hitPX);
    momYTemp.push_back(hitPY);
    momZTemp.push_back(hitPZ);
    treeHit -> GetEntry(i+1);
    postEid = hitEventID;
    //cout << "preEid = " << preEid << endl;
    //cout << "postEid = " << postEid << endl;
    //cout << "hitTime = " << hitTime << endl;
    if(preEid == eIDhit)
    {
      hitT.push_back(hitTTemp);
      hitX.push_back(hitXTemp);
      hitY.push_back(hitYTemp);
      hitZ.push_back(hitZTemp);
      momX.push_back(momXTemp);
      momY.push_back(momYTemp);
      momZ.push_back(momZTemp);

      hitTTemp.clear();
      hitXTemp.clear();
      hitYTemp.clear();
      hitZTemp.clear();
      momXTemp.clear();
      momYTemp.clear();
      momZTemp.clear();
    }

    if(preEid != postEid)
    {
      hitT.push_back(hitTTemp);
      hitX.push_back(hitXTemp);
      hitY.push_back(hitYTemp);
      hitZ.push_back(hitZTemp);
      momX.push_back(momXTemp);
      momY.push_back(momYTemp);
      momZ.push_back(momZTemp);

      hitTTemp.clear();
      hitXTemp.clear();
      hitYTemp.clear();
      hitZTemp.clear();
      momXTemp.clear();
      momYTemp.clear();
      momZTemp.clear();
    }
  }


  int trackN = hitX.size();
  for (int i = 0; i < trackN; i++)
  {
    cout <<"hitX[" << i << " ] : " ;
    for (int j = 0; j < hitX[i].size(); j++)
    {
      cout <<  hitX[i][j] << " ";
    }
    cout << endl;
  }


  // Looping over entries
  for ( int i = 0; i < trackN; i++ )
  {
    //Array for momentums
    //chekc the real Hit number ;
    //set the initial state(for Kalman filter step 0) Units: cm / GeV/c
    TVector3 pos( hitX[i][0] / 10.  , hitY[i][0] / 10.    , hitZ[i][0] / 10.  );
    TVector3 mom( momX[i][0] / 1000. , momY[i][0] / 1000.  , momZ[i][0] / 1000. );
    //Calculate the Magnitude of Initial momentum - will be use in Momentum cutting
    double momMag = mom.Mag();
    //mom.Print();
    const int pdg = pdgWanted ;     // particle pdg
    //get the charge : why devide by 3? - Curious(from Dasilva) 
    //const double charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge()/(3.);
    //std::cout << "charge  = " << charge << std::endl;
    // smeared start values
    //const bool smearPosMom = true;     // init the Reps with smeared pos and mom
    const bool smearPosMom = false;     // init the Reps with smeared pos and mom
    const double posSmear = 0.019;     // cm
    //const double momSmear = 3. /180.*TMath::Pi();     // rad
    const double momSmear = 0.019;     // rad
    //const double momMagSmear = 0.019;   // relative
    //hmomshit -> Fill( momMag* 1000. , countHits);
    TVector3 posM = pos;
    TVector3 momM = mom;
    if (smearPosMom)
    {
      posM.SetX(gRandom->Gaus(posM.X(),posSmear));
      posM.SetY(gRandom->Gaus(posM.Y(),posSmear));
      posM.SetZ(gRandom->Gaus(posM.Z(),posSmear));
      momM.SetX(gRandom->Gaus(momM.X(),momSmear));
      momM.SetY(gRandom->Gaus(momM.Y(),momSmear));
      momM.SetZ(gRandom->Gaus(momM.Z(),momSmear));
    }
    // approximate covariance(Kalman filter step.1)

    int nHits = hitX[i].size();

    TMatrixDSym covM(6);
    for (int k = 0; k < 3; ++k)
    {
      covM(k,k) = resolution * resolution;
      //covM(k,k) = pow(posSmear,2);
    }
    for (int k = 3; k < 6; ++k)
    {
      covM(k,k) = pow(resolution / (sqrt(3) * nHits) , 2);
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
    //for (int l =0; l<nHits; l++)

    for (int j = 0; j < nHits; j++)
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
      hitCoords[1] = cPos.Perp()-18.425;
      genfit::PlanarMeasurement * measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, 0., j ,nullptr );
      measurement ->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(o,u,v) ), 0 );
      fitTrack.insertPoint(new genfit::TrackPoint(measurement,&fitTrack ));
    }
    hmomB -> Fill(momMag*1000.); 
    //check
    // assert(fitTrack.checkConsistency());
    // do the fit
    fitter->processTrack(&fitTrack);
    display->addEvent(&fitTrack);
    genfit::TrackPoint* tp = fitTrack.getPointWithMeasurementAndFitterInfo(0, rep);
    if (tp == nullptr) continue;
    //genfit::KalmanFittedStateOnPlane kfsop(*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getBackwardUpdate()));
    genfit::KalmanFittedStateOnPlane kfsop(*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getForwardUpdate()));
    const TVectorD& referenceState = stateRefOrig.getState();
    const TVectorD& state = kfsop.getState();
    const TMatrixDSym covf = kfsop.getCov();
    //genfit::FitStatus *fitstat = fitTrack.getFitStatus(rep);
    //bool fitconv = fitstat->isFitConverged();
    //if(!fitconv){
    //    cout << "Track could not be fitted successfully" << endl;
    //    continue;
    //}

    momInit = (1./referenceState[0])*1000.; 
    momFitted = (1./state[0])*1000.; 
    momDiff = momInit - momFitted; 
    momRes = (momInit - momFitted) / momInit;

    tr -> Fill();

    for (int j = 0; j < window; j ++)
    {
      if (j * windowRange < momMag &&  momMag < (j+1) * windowRange)
      {
        hmomRes[j]-> Fill( ((1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
        hmomDiff[j]-> Fill(  (1./state[0])*1000. - (1./referenceState[0] )*1000.);
        momSig[j] = hmomDiff[j] ->GetStdDev();
        double trueNTemp = 0.;
        double resTemp = 0.;
        resTemp = (( ( (1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0])));
        double resMTemp = abs(res25Temp);
        posiNum[j] = hmomRes[j] -> GetEntries();
        if(resMTemp < cutR )
        {
          hmomResCut[j]> Fill( ((1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
          hmomDiffCut[j]-> Fill(  (1./state[0])*1000. - (1./referenceState[0] )*1000.);
          double fittedNTemp = 0.;
          fittedNTemp = hmomResCut[j] ->GetEntries();
          Nf[j] = fittedNTemp;
          Eff[j] = (fittedNTemp / trueNTemp) * 100.;
        }
      }
      delete fitter;
  //f ->cd();
  //tr->Write(); 
  //f ->Close();
  ////Drawing Momentum Resolution
  //TCanvas* c6 = new TCanvas();
  //c6 -> Divide(2,1); 
  //c6 -> cd(1);
  //hmomDiffAll->Draw();
  ////hmomDiffAll->Fit("gaus");
  //c6 -> cd(2);
  //hmomResAll->Draw();
  ////Drawing pval, ndf, chi2 , ch2/ndf
  //TCanvas* c8 = new TCanvas();
  //c8 -> Divide(2,2);
  //c8 -> cd(1);
  //pVal -> Draw();
  //c8 -> cd(2);
  //chI2 -> Draw();
  //c8 -> cd(3);
  //Ndf -> Draw();
  //c8 -> cd(4);
  //chI2N ->Draw();
  //TCanvas* c12 = new TCanvas();
  //c12 -> SetGrid();
  //TEfficiency * fEff = new TEfficiency("eff", "Fitting Efficiency",12,0.,300.);
  //double y[12] = {Eff25 , Eff50 , Eff75 , Eff100 , Eff125 , Eff150 , Eff175 , Eff200 , Eff225 , Eff250 , Eff275 , Eff300};
  //double Nt[12] = {Nt25 , Nt50 , Nt75 , Nt100 , Nt125 , Nt150 , Nt175 , Nt200 , Nt225 , Nt250 , Nt275 , Nt300};
  //double Nf[12] = {Nf25 , Nf50 , Nf75 , Nf100 , Nf125 , Nf150 , Nf175 , Nf200 , Nf225 , Nf250 , Nf275 , Nf300};
  //double x[12]  ;
  //double ex[12] ;
  //double ey[12] ;
  //for (int i =0; i < 12 ; i++)
  //{
  //  x[i] = 25. * (i+1) - 12.5;
  //  ex[i] = 12.5;
  //  ey[i] = y[i]* TMath::Sqrt( ( 1 / Nf[i] ) + ( 1/ Nt[i]) );
  //}
  //auto Eff = new TGraphErrors(12,x,y,ex,ey);
  //Eff ->SetTitle("Efficiency");
  //Eff -> SetPointError( 25. , 0.);
  //Eff -> SetMarkerColor(2);
  //Eff ->GetXaxis() ->SetTitle("p_{initial} [MeV/c]");
  //Eff ->GetYaxis() ->SetTitle("Efficiency [%]");
  //Eff -> Draw("AP");
  //TCanvas* c14 = new TCanvas();
  //c14 -> SetGrid();
  //hmomResO200 ->GetXaxis() ->SetTitle("Momentum resolution");
  //hmomResO200 ->GetYaxis() ->SetTitle("# of Track(positron)");
  //hmomResO200 ->Draw();
  //TCanvas* c20 = new TCanvas();
  //c20 -> SetGrid();
  //double sig[12] = {momSig25 , momSig50 , momSig75 , momSig100 , momSig125 , momSig150 , momSig175 , momSig200 , momSig225 , momSig250 , momSig275 , momSig300};
  //double eSig[12] ;
  //double eMomSig[12];
  //for (int i =0; i < 12 ; i++)
  //{
  //  x[i] = 25. * (i+1) - 12.5;
  //  eMomSig[i] = sig[i] / x[i];
  //  eSig[i] = eMomSig[i]* TMath::Sqrt( 1/ Nt[i] );
  //}
  //auto mSig = new TGraphErrors(12 , x , eMomSig , ex , eSig);
  //mSig ->SetTitle("Performance");
  //mSig -> SetPointError( 25. , 0.);
  //mSig -> SetMarkerColor(2);
  //mSig ->GetXaxis() ->SetTitle("p_{initial} [MeV/c]");
  //mSig ->GetYaxis() ->SetTitle("#sigma_{p} / p_{initial}");
  //mSig -> Draw("AP");
  //// csv start 
  //{
  //  double * xE;
  //  double * yE;
  //  ofstream csv("Efficiency.csv");
  //  int nEff = Eff -> GetN();
  //  csv  << "Efficiency" << endl;
  //  for(Int_t i = 0; i < nEff ; i++)
  //  {
  //    xE = Eff -> GetX();
  //    yE = Eff -> GetY();
  //    csv << xE[i] << "," << yE[i] << endl;
  //  }
  //  csv << "Efficiency, 200MeV/c < p_{initial < 275MeV/c}" << endl;
  //  csv << EffO200 << endl;
  //  csv << "# of positrons" << endl;
  //  csv << "~25"<< "," << "~50" << "," << "50~75" << ","<< "75~100" << ","<< "100~125" << ","<< "125~150" << ","<< "150~175" << ","<< "175~200" << ","<< "200~225" << ","<< "225~250" << ","<< "250~275" << ","<< "275~300" << endl;
  //  csv << Nt25 << Nt50 << "," << Nt75 << "," << Nt100 << "," << Nt125 << "," << Nt150 << "," << Nt175 << "," << Nt200 << "," << Nt225 << "," << Nt250 << "," <<  Nt275 << "," << Nt300 << endl;  
  //  csv.close();
  //}
  //// ////////////////////////////////////////////////CSV is done/////////////////////////////////////////////////////////

    }
  display->setOptions("ABDEFHMPT"); // G show geometry
  //if (matFX) display->setOptions("ABDEFGHMPT"); // G show geometry
  display->open();


}
