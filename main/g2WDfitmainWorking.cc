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
  double Effa[12];
  double resolution =  pitch ;
  gRandom->SetSeed(14);
  // init MeasurementCreator
  genfit::MeasurementCreator measurementCreator;
  // init geometry and mag. field
  new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("/home/wdlee/ssd1/work/build/g2wdPGID_build/g2wdGeom.gdml");
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

  TH1D *hmomRes25 = new TH1D("P_{Initial} < 25 MeV/c","Momentum Resolution",nbin1, -mRr , mRr);
  TH1D *hmomResCut25 = new TH1D("P_{Initial} < 25 MeV/c","Momentum Resolution Cut",nbin1, -mRr , mRr);
  TH1D *hmomDiff25 = new TH1D("P_{Initial} < 25 MeV/c","Momentum Difference",nbin1, -mDr , mDr);
  TH1D *hmomDiffCut25 = new TH1D("P_{Initial} < 25 MeV/c","Momentum Difference Cut",nbin1, -mDr , mDr);
  TH1D *hmomSig25 = new TH1D("P_{Initial} < 25 MeV/c","Momentum Difference / true Momentum",nbin1, -mDr , mDr);

  TH1D *hmomRes50 = new TH1D("25 MeV/c < P_{Initial} < 50 MeV/c","Momentum Resolution",nbin1, -mRr , mRr);
  TH1D *hmomResCut50 = new TH1D("25 MeV/c < P_{Initial} < 50 MeV/c","Momentum Resolution Cut",nbin1, -mRr , mRr);
  TH1D *hmomDiff50 = new TH1D("25 MeV/c < P_{Initial} < 50 MeV/c","Momentum Difference",nbin1, -mDr , mDr);
  TH1D *hmomDiffCut50 = new TH1D("25 MeV/c < P_{Initial} < 50 MeV/c","Momentum Difference Cut",nbin1, -mDr , mDr);
  TH1D *hmomSig50 = new TH1D("25 MeV/c < P_{Initial} < 50 MeV/c","Momentum Difference / true Momentum",nbin1, -mDr , mDr);

  TH1D *hmomRes75 = new TH1D("50 MeV/c < P_{Initial} < 75 MeV/c","Momentum Resolution",nbin1, -mRr , mRr);
  TH1D *hmomResCut75 = new TH1D("50 MeV/c < P_{Initial} < 75 MeV/c","Momentum Resolution Cut",nbin1, -mRr , mRr);
  TH1D *hmomDiff75 = new TH1D("50 MeV/c < P_{Initial} < 75 MeV/c","Momentum Difference",nbin1, -mDr , mDr);
  TH1D *hmomDiffCut75 = new TH1D("50 MeV/c < P_{Initial} < 75 MeV/c","Momentum Difference Cut",nbin1, -mDr , mDr);
  TH1D *hmomSig75 = new TH1D("50 MeV/c < P_{Initial} < 75 MeV/c","Momentum Difference / true Momentum",nbin1, -mDr , mDr);

  TH1D *hmomRes100 = new TH1D("75 MeV/c < P_{Initial} < 100 MeV/c","Momentum Resolution",nbin1, -mRr , mRr);
  TH1D *hmomResCut100 = new TH1D("75 MeV/c < P_{Initial} < 100 MeV/c","Momentum Resolution Cut",nbin1, -mRr , mRr);
  TH1D *hmomDiff100 = new TH1D("75 MeV/c < P_{Initial} < 100 MeV/c","Momentum Difference",nbin1, -mDr , mDr);
  TH1D *hmomDiffCut100 = new TH1D("75 MeV/c < P_{Initial} < 100 MeV/c","Momentum Difference Cut",nbin1, -mDr , mDr);
  TH1D *hmomSig100 = new TH1D("75 MeV/c < P_{Initial} < 100 MeV/c","Momentum Difference / true Momentum",nbin1, -mDr , mDr);

  TH1D *hmomRes125 = new TH1D("100 MeV/c < P_{Initial} < 125 MeV/c","Momentum Resolution",nbin1, -mRr , mRr);
  TH1D *hmomResCut125 = new TH1D("100 MeV/c < P_{Initial} < 125 MeV/c","Momentum Resolution Cut",nbin1, -mRr , mRr);
  TH1D *hmomDiff125 = new TH1D("100 MeV/c < P_{Initial} < 125 MeV/c","Momentum Difference",nbin1, -mDr , mDr);
  TH1D *hmomDiffCut125 = new TH1D("100 MeV/c < P_{Initial} < 125 MeV/c","Momentum Difference Cut",nbin1, -mDr , mDr);
  TH1D *hmomSig125 = new TH1D("100 MeV/c < P_{Initial} < 125 MeV/c","Momentum Difference / true Momentum",nbin1, -mDr , mDr);

  TH1D *hmomRes150 = new TH1D("125 MeV/c < P_{Initial} < 150 MeV/c","Momentum Resolution",nbin1, -mRr , mRr);
  TH1D *hmomResCut150 = new TH1D("125 MeV/c < P_{Initial} < 150 MeV/c","Momentum Resolution Cut",nbin1, -mRr , mRr);
  TH1D *hmomDiff150 = new TH1D("125 MeV/c < P_{Initial} < 150 MeV/c","Momentum Difference",nbin1, -mDr , mDr);
  TH1D *hmomDiffCut150 = new TH1D("125 MeV/c < P_{Initial} < 150 MeV/c","Momentum Difference Cut",nbin1, -mDr , mDr);
  TH1D *hmomSig150 = new TH1D("125 MeV/c < P_{Initial} < 150 MeV/c","Momentum Difference / true Momentum",nbin1, -mDr , mDr);


  TH1D *hmomRes175 = new TH1D("150 MeV/c < P_{Initial} < 175 MeV/c","Momentum Resolution",nbin1, -mRr , mRr);
  TH1D *hmomResCut175 = new TH1D("150 MeV/c < P_{Initial} < 175 MeV/c","Momentum Resolution Cut",nbin1, -mRr , mRr);
  TH1D *hmomDiff175 = new TH1D("150 MeV/c < P_{Initial} < 175 MeV/c","Momentum Difference",nbin1, -mDr , mDr);
  TH1D *hmomDiffCut175 = new TH1D("150 MeV/c < P_{Initial} < 175 MeV/c","Momentum Difference Cut",nbin1, -mDr , mDr);
  TH1D *hmomSig175 = new TH1D("150 MeV/c < P_{Initial} < 175 MeV/c","Momentum Difference / true Momentum",nbin1, -mDr , mDr);


  TH1D *hmomRes200 = new TH1D("175 MeV/c < P_{Initial} < 200 MeV/c","Momentum Resolution",nbin1, -mRr , mRr);
  TH1D *hmomResCut200 = new TH1D("175 MeV/c < P_{Initial} < 200 MeV/c","Momentum Resolution Cut",nbin1, -mRr , mRr);
  TH1D *hmomDiff200 = new TH1D("175 MeV/c < P_{Initial} < 200 MeV/c","Momentum Difference",nbin1, -mDr , mDr);
  TH1D *hmomDiffCut200 = new TH1D("175 MeV/c < P_{Initial} < 200 MeV/c","Momentum Difference Cut",nbin1, -mDr , mDr);
  TH1D *hmomSig200 = new TH1D("175 MeV/c < P_{Initial} < 200 MeV/c","Momentum Difference / true Momentum",nbin1, -mDr , mDr);


  TH1D *hmomRes225 = new TH1D("200 MeV/c < P_{Initial} < 225 MeV/c","Momentum Resolution",nbin1, -mRr , mRr);
  TH1D *hmomResCut225 = new TH1D("200 MeV/c < P_{Initial} < 225 MeV/c","Momentum Resolution Cut",nbin1, -mRr , mRr);
  TH1D *hmomDiff225 = new TH1D("200 MeV/c < P_{Initial} < 225 MeV/c","Momentum Difference",nbin1, -mDr , mDr);
  TH1D *hmomDiffCut225 = new TH1D("200 MeV/c < P_{Initial} < 225 MeV/c","Momentum Difference Cut",nbin1, -mDr , mDr);
  TH1D *hmomSig225 = new TH1D("200 MeV/c < P_{Initial} < 225 MeV/c","Momentum Difference / true Momentum",nbin1, -mDr , mDr);


  TH1D *hmomRes250 = new TH1D("225 MeV/c < P_{Initial} < 250 MeV/c","Momentum Resolution",nbin1, -mRr , mRr);
  TH1D *hmomResCut250 = new TH1D("225 MeV/c < P_{Initial} < 250 MeV/c","Momentum Resolution Cut",nbin1, -mRr , mRr);
  TH1D *hmomDiff250 = new TH1D("225 MeV/c < P_{Initial} < 250 MeV/c","Momentum Difference",nbin1, -mDr , mDr);
  TH1D *hmomDiffCut250 = new TH1D("225 MeV/c < P_{Initial} < 250 MeV/c","Momentum Difference Cut",nbin1, -mDr , mDr);
  TH1D *hmomSig250 = new TH1D("225 MeV/c < P_{Initial} < 250 MeV/c","Momentum Difference / true Momentum",nbin1, -mDr , mDr);

  TH1D *hmomRes275 = new TH1D("250 MeV/c < P_{Initial} < 275 MeV/c","Momentum Resolution",nbin1, -mRr , mRr);
  TH1D *hmomResCut275 = new TH1D("250 MeV/c < P_{Initial} < 275 MeV/c","Momentum Resolution Cut",nbin1, -mRr , mRr);
  TH1D *hmomDiff275 = new TH1D("250 MeV/c < P_{Initial} < 275 MeV/c","Momentum Difference",nbin1, -mDr , mDr);
  TH1D *hmomDiffCut275 = new TH1D("250 MeV/c < P_{Initial} < 275 MeV/c","Momentum Difference Cut",nbin1, -mDr , mDr);
  TH1D *hmomSig275 = new TH1D("250 MeV/c < P_{Initial} < 275 MeV/c","Momentum Difference / true Momentum",nbin1, -mDr , mDr);

  TH1D *hmomRes300 = new TH1D("275 MeV/c < P_{Initial} < 300 MeV/c","Momentum Resolution",nbin1, -mRr , mRr);
  TH1D *hmomResCut300 = new TH1D("275 MeV/c < P_{Initial} < 300 MeV/c","Momentum Resolution Cut",nbin1, -mRr , mRr);
  TH1D *hmomDiff300 = new TH1D("275 MeV/c < P_{Initial} < 300 MeV/c","Momentum Difference",nbin1, -mDr , mDr);
  TH1D *hmomDiffCut300 = new TH1D("275 MeV/c < P_{Initial} < 300 MeV/c","Momentum Difference Cut",nbin1, -mDr , mDr);
  TH1D *hmomSig300 = new TH1D("275 MeV/c < P_{Initial} < 300 MeV/c","Momentum Difference / true Momentum",nbin1, -mDr , mDr);

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

  double Eff25 = 0.;
  double Eff50 = 0.;
  double Eff75 = 0.;
  double Eff100 = 0.;
  double Eff125 = 0.;
  double Eff150 = 0.; 
  double Eff175 = 0.;
  double Eff200 = 0.;
  double Eff225 = 0.;
  double Eff250 = 0.;
  double Eff275 = 0.;
  double Eff300 = 0.;

  double EffO200 = 0.;

  double Ntall = 0.;
  double Nt25 = 0.;
  double Nt50 = 0.;
  double Nt75 = 0.;
  double Nt100 = 0.;
  double Nt125 = 0.;
  double Nt150 = 0.;
  double Nt175 = 0.;
  double Nt200 = 0.;
  double Nt225 = 0.;
  double Nt250 = 0.;
  double Nt275 = 0.;
  double Nt300 = 0.;

  double Nfall = 0.;
  double Nf25 = 0.;
  double Nf50 = 0.;
  double Nf75 = 0.;
  double Nf100 = 0.;
  double Nf125 = 0.;
  double Nf150 = 0.;
  double Nf175 = 0.;
  double Nf200 = 0.;
  double Nf225 = 0.;
  double Nf250 = 0.;
  double Nf275 = 0.;
  double Nf300 = 0.;

  double momSig25 = 0.;
  double momSig50 = 0.;
  double momSig75 = 0.;
  double momSig100= 0.;
  double momSig125= 0.;
  double momSig150= 0.;
  double momSig175= 0.;
  double momSig200= 0.;
  double momSig225= 0.;
  double momSig250= 0.;
  double momSig275= 0.;
  double momSig300= 0.;


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


    if (0.2 < momMag && momMag <  0.275 )
    {
      hmomResO200-> Fill( ((1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));

      double resO200 = 0.;
      resO200 = (( ( (1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0])));
      double resO200A = abs(resO200);
      if(resO200A < 0.3 )
      {
        hmomResCutO200-> Fill( ((1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
      }
    }


    double cutR = 0.2;
    if (0.025 > momMag )
    {
      hmomRes25-> Fill( ((1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
      hmomDiff25-> Fill(  (1./state[0])*1000. - (1./referenceState[0] )*1000.);
      momSig25 = hmomDiff25 ->GetStdDev();
      double trueN25 = 0.;
      double res25 = 0.;
      res25 = (( ( (1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0])));
      double resM25 = abs(res25);
      trueN25 = hmomRes25 -> GetEntries();
      Nt25 = trueN25; 
      if(resM25 < cutR )
      {
        hmomResCut25-> Fill( ((1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
        hmomDiffCut25-> Fill(  (1./state[0])*1000. - (1./referenceState[0] )*1000.);
        double fittedN25 = 0.;
        fittedN25 = hmomResCut25 ->GetEntries();
        Nf25 = fittedN25;
        Eff25 = (fittedN25 / trueN25) * 100.;  
        Effa[0] = Eff25;
        Effv.push_back(Eff25);
      }
    }
    if ( 0.025 < momMag && momMag < 0.05)
    {
      hmomRes50-> Fill( ((1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
      hmomDiff50-> Fill(  (1./state[0])*1000. - (1./referenceState[0] )*1000.);
      hmomSig50-> Fill(  ((1./state[0])*1000. - (1./referenceState[0] )*1000.) / (1./referenceState[0] )*1000.);
      momSig50 = hmomDiff50 ->GetStdDev();
      tr -> Fill();
      double trueN50 = 0.;
      double res50 = 0.;
      res50 = (( ( (1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0])));
      double resM50 = abs(res50);
      trueN50 = hmomRes50 -> GetEntries();
      Nt50 = trueN50; 
      if(resM50 < cutR )
      {
        hmomResCut50-> Fill( ((1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
        hmomDiffCut50-> Fill(  (1./state[0])*1000. - (1./referenceState[0] )*1000.);
        double fittedN50 = 0.;
        fittedN50 = hmomResCut50 ->GetEntries();
        Nf50 = fittedN50;
        Eff50 = (fittedN50 / trueN50) * 100.;  
        Effa[1] = Eff50;
        Effv.push_back(Eff50);
      }
    }
    if (0.05 < momMag && momMag< 0.075)
    { 
      hmomRes75-> Fill( ( (1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
      hmomDiff75-> Fill(  (1./state[0])*1000. - (1./referenceState[0] )*1000.);
      hmomSig75-> Fill(  ((1./state[0])*1000. - (1./referenceState[0] )*1000.) / (1./referenceState[0] )*1000.);
      momSig75 = hmomDiff75 ->GetStdDev();
      double trueN75 = 0.;
      double res75 = 0.;
      res75 = (( ( (1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0])));
      double resM75 = abs(res75);
      trueN75 = hmomRes75 -> GetEntries();
      Nt75 = trueN75; 
      if(resM75 < cutR )
      {
        hmomResCut75-> Fill( ((1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
        hmomDiffCut75-> Fill(  (1./state[0])*1000. - (1./referenceState[0] )*1000.);
        double fittedN75 = 0.;
        fittedN75 = hmomResCut75 ->GetEntries();
        Nf75 = fittedN75;
        Eff75 = (fittedN75 / trueN75) * 100;  
        Effa[2] = Eff75;
        Effv.push_back(Eff75);
      }
    }
    if (0.075 < momMag && momMag< 0.1)
    { 
      hmomRes100-> Fill( ( (1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
      hmomDiff100-> Fill(  (1. / state[0])*1000. - (1. / referenceState[0] )*1000.);
      hmomSig100-> Fill(  ((1./state[0])*1000. - (1./referenceState[0] )*1000.) / (1./referenceState[0] )*1000.);
      momSig100 = hmomDiff100 ->GetStdDev();
      double trueN100 = 0.;
      double res100 = 0.;
      res100 = (( ( (1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0])));
      double resM100 = abs(res100);
      trueN100 = hmomRes100 -> GetEntries();
      Nt100 = trueN100; 
      if(resM100 < cutR )
      {
        hmomResCut100-> Fill( ((1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
        hmomDiffCut100-> Fill(  (1./state[0])*1000. - (1./referenceState[0] )*1000.);
        double fittedN100;
        fittedN100 = hmomResCut100 ->GetEntries();
        Nf100 = fittedN100;
        Eff100 = (fittedN100 / trueN100) * 100;  
        Effa[3] = Eff100;
        Effv.push_back(Eff100);
      } 
    }
    if (0.1 <  momMag && momMag< 0.125)
    { 
      hmomRes125-> Fill( ( (1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
      hmomDiff125-> Fill(  (1./state[0])*1000. - (1./referenceState[0] )*1000.);
      hmomSig125-> Fill(  ((1./state[0])*1000. - (1./referenceState[0] )*1000.) / (1./referenceState[0] )*1000.);
      momSig125 = hmomDiff125 ->GetStdDev();
      double trueN125 =0.;
      double res125 = 0.;
      res125 = (( ( (1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0])));
      double resM125 = abs(res125);
      trueN125 = hmomRes125 -> GetEntries();
      Nt125 = trueN125; 
      if(resM125 < cutR )
      {
        hmomResCut125-> Fill( ((1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
        hmomDiffCut125-> Fill(  (1./state[0])*1000. - (1./referenceState[0] )*1000.);
        double fittedN125 = 0.;
        fittedN125 = hmomResCut125 ->GetEntries();
        Nf125 = fittedN125;
        Eff125 = (fittedN125 / trueN125) * 100;  
        Effa[4] = Eff125;
        Effv.push_back(Eff125);
      } 
    }
    if (0.125 < momMag && momMag<0.15)
    { 
      hmomRes150-> Fill(  ((1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
      hmomDiff150-> Fill(  (1./state[0])*1000. - (1./referenceState[0] )*1000.);
      hmomSig150-> Fill(  ((1./state[0])*1000. - (1./referenceState[0] )*1000.) / (1./referenceState[0] )*1000.);
      momSig150 = hmomDiff150 ->GetStdDev();
      double trueN150 = 0.;
      double res150 = 0.;
      res150 = (( ( (1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0])));
      double resM150 = abs(res150);
      trueN150 = hmomRes150 -> GetEntries();
      Nt150 = trueN150; 
      if(resM150 < cutR )
      {
        hmomResCut150-> Fill( ((1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
        hmomDiffCut150-> Fill(  (1./state[0])*1000. - (1./referenceState[0] )*1000.);
        double fittedN150 = 0.;
        fittedN150 = hmomResCut150 ->GetEntries();
        Nf150 = fittedN150;
        Eff150 = (fittedN150 / trueN150) * 100;  
        Effa[5] = Eff150;
        Effv.push_back(Eff150);
      } 
    }
    if (0.15 < momMag && momMag<0.175)
    { 
      hmomRes175-> Fill(  ((1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
      hmomDiff175-> Fill(  (1./state[0])*1000. - (1./referenceState[0] )*1000.);
      hmomSig175-> Fill(  ((1./state[0])*1000. - (1./referenceState[0] )*1000.) / (1./referenceState[0] )*1000.);
      momSig175 = hmomDiff175 ->GetStdDev();
      double trueN175 = 0.;
      double res175 = 0.;
      res175 = (( ( (1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0])));
      double resM175 = abs(res175);
      trueN175 = hmomRes175 -> GetEntries();
      Nt175 = trueN175; 
      if(resM175 < cutR )
      {
        hmomResCut175-> Fill( ((1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
        hmomDiffCut175-> Fill(  (1./state[0])*1000. - (1./referenceState[0] )*1000.);
        double fittedN175 = 0.;
        fittedN175 = hmomResCut175 ->GetEntries();
        Nf175 = fittedN175;
        Eff175 = (fittedN175 / trueN175) * 100;  
        Effa[6] = Eff175;
        Effv.push_back(Eff175);
      } 
    }
    if (0.175 < momMag && momMag<0.2)
    { 
      hmomRes200-> Fill(  ((1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
      hmomDiff200-> Fill(  (1./state[0])*1000. - (1./referenceState[0] )*1000.);
      hmomSig200-> Fill(  ((1./state[0])*1000. - (1./referenceState[0] )*1000.) / (1./referenceState[0] )*1000.);
      momSig200 = hmomDiff200 ->GetStdDev();
      double trueN200 =0.;
      double res200 = 0.;
      res200 = (( ( (1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0])));
      double resM200 = abs(res200);
      trueN200 = hmomRes200 -> GetEntries();
      Nt200 = trueN200; 
      if(resM200 < cutR )
      {
        hmomResCut200-> Fill( ((1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
        hmomDiffCut200-> Fill(  (1./state[0])*1000. - (1./referenceState[0] )*1000.);
        double fittedN200 = 0.;
        fittedN200 = hmomResCut200 ->GetEntries();
        Nf200 = fittedN200;
        Eff200 = (fittedN200 / trueN200) * 100;  
        Effa[7] = Eff200;
        Effv.push_back(Eff200);
      } 
    }
    if (0.2 < momMag && momMag<0.225)
    { 
      hmomRes225-> Fill(  ((1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
      hmomDiff225-> Fill(  (1./state[0])*1000. - (1./referenceState[0] )*1000.);
      hmomSig225-> Fill(  ((1./state[0])*1000. - (1./referenceState[0] )*1000.) / (1./referenceState[0] )*1000.);
      momSig225 = hmomDiff225 ->GetStdDev();
      double trueN225 = 0.;
      double res225 = 0.;
      res225 = (( ( (1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0])));
      double resM225 = abs(res225);
      trueN225 = hmomRes225 -> GetEntries();
      Nt225 = trueN225; 
      if(resM225 < cutR )
      {
        hmomResCut225-> Fill( ((1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
        hmomDiffCut225-> Fill(  (1./state[0])*1000. - (1./referenceState[0] )*1000.);
        double fittedN225 = 0.;
        fittedN225 = hmomResCut225 ->GetEntries();
        Nf225 = fittedN225;
        Eff225 = (fittedN225 / trueN225) * 100;  
        Effa[8] = Eff225;
        Effv.push_back(Eff225);
      } 
    }
    if (0.225 < momMag && momMag<0.25)
    { 
      hmomRes250-> Fill(  ((1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
      hmomDiff250-> Fill(  (1./state[0])*1000. - (1./referenceState[0] )*1000.);
      hmomSig250-> Fill(  ((1./state[0])*1000. - (1./referenceState[0] )*1000.) / (1./referenceState[0] )*1000.);
      momSig250 = hmomDiff250 ->GetStdDev();
      double trueN250 = 0.;
      double res250 = 0.;
      res250 = (( ( (1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0])));
      double resM250 = abs(res250);
      trueN250 = hmomRes250 -> GetEntries();
      Nt250 = trueN250; 
      if(resM250 < cutR )
      {
        hmomResCut250-> Fill( ((1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
        hmomDiffCut250-> Fill(  (1./state[0])*1000. - (1./referenceState[0] )*1000.);
        double fittedN250;
        fittedN250 = hmomResCut250 ->GetEntries();
        Nf250 = fittedN250;
        Eff250 = (fittedN250 / trueN250) * 100;  
        Effa[9] = Eff250;
        Effv.push_back(Eff250);
      } 
    }
    if (0.25 < momMag && momMag<0.275)
    { 
      hmomRes275-> Fill(  ((1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
      hmomDiff275-> Fill(  (1./state[0])*1000. - (1./referenceState[0] )*1000.);
      hmomSig275-> Fill(  ((1./state[0])*1000. - (1./referenceState[0] )*1000.) / (1./referenceState[0] )*1000.);
      momSig275 = hmomDiff275 ->GetStdDev();
      double trueN275 = 0.;
      double res275 = 0.;
      res275= (( ( (1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0])));
      double resM275 = abs(res275);
      trueN275 = hmomRes275 -> GetEntries();
      Nt275 = trueN275; 
      if(resM275 < cutR )
      {
        hmomResCut275-> Fill( ((1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
        hmomDiffCut275-> Fill(  (1./state[0])*1000. - (1./referenceState[0] )*1000.);
        double fittedN275 = 0.;
        fittedN275 = hmomResCut275 ->GetEntries();
        Nf275 = fittedN275;
        Eff275 = (fittedN275 / trueN275) * 100;  
        Effa[10] = Eff275;
        Effv.push_back(Eff275);
      } 
    }
    if (0.275 < momMag && momMag<0.3)
    { 
      hmomRes300-> Fill(  ((1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
      hmomDiff300-> Fill(  (1./state[0])*1000. - (1./referenceState[0] )*1000.);
      hmomSig300-> Fill(  ((1./state[0])*1000. - (1./referenceState[0] )*1000.) / (1./referenceState[0] )*1000.);
      momSig300 = hmomDiff300 ->GetStdDev();
      double trueN300 = 0.;
      double res300 = 0.;
      res300= (( ( (1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0])));
      double resM300 = abs(res300);
      trueN300 = hmomRes300 -> GetEntries();
      Nt300 = trueN300;
      if(resM300 < cutR )
      {
        hmomResCut300-> Fill( ((1./state[0]) - (1./referenceState[0] )) / (1./referenceState[0]));
        hmomDiffCut300-> Fill(  (1./state[0])*1000. - (1./referenceState[0] )*1000.);
        double fittedN300 = 0.;
        fittedN300 = hmomResCut300 ->GetEntries();
        Nf300 = fittedN300;
        Eff300 = (fittedN300 / trueN300) * 100;  
        Effa[11] = Eff300;
        Effv.push_back(Eff300);
      } 
    }
    momFitting.push_back((1./ state[0])*1000.);
    momInitial.push_back((1./referenceState[0])*1000.);

  }
  delete fitter;

  f ->cd();
  tr->Write(); 
  f ->Close();

  //Drawing Momentum Resolution

  TCanvas* c6 = new TCanvas();
  c6 -> Divide(2,1); 
  c6 -> cd(1);
  hmomDiffAll->Draw();
  //hmomDiffAll->Fit("gaus");
  c6 -> cd(2);
  hmomResAll->Draw();

  //Drawing pval, ndf, chi2 , ch2/ndf
  TCanvas* c8 = new TCanvas();
  c8 -> Divide(2,2);
  c8 -> cd(1);
  pVal -> Draw();
  c8 -> cd(2);
  chI2 -> Draw();
  c8 -> cd(3);
  Ndf -> Draw();
  c8 -> cd(4);
  chI2N ->Draw();

  TCanvas* c12 = new TCanvas();
  c12 -> SetGrid();

  TEfficiency * fEff = new TEfficiency("eff", "Fitting Efficiency",12,0.,300.);

  double y[12] = {Eff25 , Eff50 , Eff75 , Eff100 , Eff125 , Eff150 , Eff175 , Eff200 , Eff225 , Eff250 , Eff275 , Eff300};
  double Nt[12] = {Nt25 , Nt50 , Nt75 , Nt100 , Nt125 , Nt150 , Nt175 , Nt200 , Nt225 , Nt250 , Nt275 , Nt300};
  double Nf[12] = {Nf25 , Nf50 , Nf75 , Nf100 , Nf125 , Nf150 , Nf175 , Nf200 , Nf225 , Nf250 , Nf275 , Nf300};
  double x[12]  ;
  double ex[12] ;
  double ey[12] ;
  for (int i =0; i < 12 ; i++)
  {
    x[i] = 25. * (i+1) - 12.5;
    ex[i] = 12.5;
    ey[i] = y[i]* TMath::Sqrt( ( 1 / Nf[i] ) + ( 1/ Nt[i]) );
  }

  auto Eff = new TGraphErrors(12,x,y,ex,ey);
  Eff ->SetTitle("Efficiency");
  Eff -> SetPointError( 25. , 0.);
  Eff -> SetMarkerColor(2);
  Eff ->GetXaxis() ->SetTitle("p_{initial} [MeV/c]");
  Eff ->GetYaxis() ->SetTitle("Efficiency [%]");

  Eff -> Draw("AP");

  TCanvas* c14 = new TCanvas();
  c14 -> SetGrid();
  hmomResO200 ->GetXaxis() ->SetTitle("Momentum resolution");
  hmomResO200 ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmomResO200 ->Draw();


  TCanvas* c15 = new TCanvas();
  c15 -> SetGrid();
  hmomResCutO200 ->GetXaxis() ->SetTitle("Momentum resolution");
  hmomResCutO200 ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmomResCutO200 ->Draw();



  TCanvas * c16 = new TCanvas();
  c16 -> Divide(2,2);
  c16 -> cd(1);
  hmomref ->GetXaxis() ->SetTitle("p_{initial}[MeV/c]");
  hmomref ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmomref -> Draw();
  c16 -> cd(2);
  hmoms ->GetXaxis() ->SetTitle("p_{Fitting}[MeV/c]");
  hmoms ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmoms -> Draw();
  c16 -> cd(3);
  hmomrefCut ->GetXaxis() ->SetTitle("p_{initial}[MeV/c]");
  hmomrefCut ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmomrefCut -> Draw();
  c16 -> cd(4);
  hmomsCut ->GetXaxis() ->SetTitle("p_{Fitting}[MeV/c]");
  hmomsCut ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmomsCut -> Draw();

  //TCanvas * c17 = new TCanvas();
  //hmomshit ->GetXaxis() ->SetTitle("p_{initial} [MeV/c]");
  //hmomshit ->GetYaxis() ->SetTitle("# of hits");
  //hmomshit -> Draw("colz");

  TCanvas * c18 = new TCanvas();
  c18 -> Divide(2,2);
  c18 -> cd(1); 
  hmomNeg ->GetXaxis() ->SetTitle("p_{Fitting}[MeV/c]");
  hmomNeg ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmomNeg -> Draw();     
  c18 -> cd(2); 
  hmomrefNeg ->GetXaxis() ->SetTitle("p_{initial}[MeV/c]");
  hmomrefNeg ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmomrefNeg -> Draw();
  c18 -> cd(3); 
  hmomdiffNeg ->GetXaxis() ->SetTitle("p_{difference}[MeV/c]");
  hmomdiffNeg ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmomdiffNeg -> Draw();   
  c18 -> cd(4);
  hmomLow ->GetXaxis() ->SetTitle("p_{initial}[MeV/c]");
  hmomLow ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmomLow -> Draw();

  TCanvas * c19 = new TCanvas();
  c19 -> Divide(3,4);
  c19 -> cd(1);
  hmomDiff25 ->GetXaxis() ->SetTitle("p_{difference} / p_{initial}[MeV/c]");
  hmomDiff25 ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmomDiff25 -> Draw();

  c19 -> cd(2);
  hmomDiff50 ->GetXaxis() ->SetTitle("p_{difference} / p_{initial}[MeV/c]");
  hmomDiff50 ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmomDiff50 -> Draw();

  c19 -> cd(3);
  hmomDiff75 ->GetXaxis() ->SetTitle("p_{difference} / p_{initial}[MeV/c]");
  hmomDiff75 ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmomDiff75 -> Draw();


  c19 -> cd(4);
  hmomDiff100 ->GetXaxis() ->SetTitle("p_{difference} / p_{initial}[MeV/c]");
  hmomDiff100 ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmomDiff100 -> Draw();

  c19 -> cd(5);
  hmomDiff125 ->GetXaxis() ->SetTitle("p_{difference} / p_{initial}[MeV/c]");
  hmomDiff125 ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmomDiff125 -> Draw();


  c19 -> cd(6);
  hmomDiff150 ->GetXaxis() ->SetTitle("p_{difference} / p_{initial}[MeV/c]");
  hmomDiff150 ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmomDiff150 -> Draw();


  c19 -> cd(7);
  hmomDiff175 ->GetXaxis() ->SetTitle("p_{difference} / p_{initial}[MeV/c]");
  hmomDiff175 ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmomDiff175 -> Draw();


  c19 -> cd(8);
  hmomDiff200 ->GetXaxis() ->SetTitle("p_{difference} / p_{initial}[MeV/c]");
  hmomDiff200 ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmomDiff200 -> Draw();

  c19 -> cd(9);
  hmomDiff225 ->GetXaxis() ->SetTitle("p_{difference} / p_{initial}[MeV/c]");
  hmomDiff225 ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmomDiff225 -> Draw();


  c19 -> cd(10);
  hmomDiff250 ->GetXaxis() ->SetTitle("p_{difference} / p_{initial}[MeV/c]");
  hmomDiff250 ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmomDiff250 -> Draw();

  c19 -> cd(11);
  hmomDiff275 ->GetXaxis() ->SetTitle("p_{difference} / p_{initial}[MeV/c]");
  hmomDiff275 ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmomDiff275 -> Draw();

  c19 -> cd(12);
  hmomDiff300 ->GetXaxis() ->SetTitle("p_{difference} / p_{initial}[MeV/c]");
  hmomDiff300 ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmomDiff300 -> Draw();


  TCanvas* c20 = new TCanvas();
  c20 -> SetGrid();
  double sig[12] = {momSig25 , momSig50 , momSig75 , momSig100 , momSig125 , momSig150 , momSig175 , momSig200 , momSig225 , momSig250 , momSig275 , momSig300};
  double eSig[12] ;
  double eMomSig[12];
  for (int i =0; i < 12 ; i++)
  {
    x[i] = 25. * (i+1) - 12.5;
    eMomSig[i] = sig[i] / x[i];
    eSig[i] = eMomSig[i]* TMath::Sqrt( 1/ Nt[i] );
  }
  auto mSig = new TGraphErrors(12 , x , eMomSig , ex , eSig);
  mSig ->SetTitle("Performance");


  mSig -> SetPointError( 25. , 0.);
  mSig -> SetMarkerColor(2);
  mSig ->GetXaxis() ->SetTitle("p_{initial} [MeV/c]");
  mSig ->GetYaxis() ->SetTitle("#sigma_{p} / p_{initial}");

  mSig -> Draw("AP");





  TCanvas * c21 = new TCanvas();
  c21 -> Divide(2,2);
  c21 -> cd(1);
  hmomDiff200 ->GetXaxis() ->SetTitle("p_{difference}");
  hmomDiff200 ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmomDiff200 -> Draw();
  hmomDiff200 -> Fit("gaus");



  c21 -> cd(2);
  hmomDiff225 ->GetXaxis() ->SetTitle("p_{difference}");
  hmomDiff225 ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmomDiff225 -> Draw();
  hmomDiff225 -> Fit("gaus");

  c21 -> cd(2);
  hmomDiff250 ->GetXaxis() ->SetTitle("p_{difference}");
  hmomDiff250 ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmomDiff250 -> Draw();
  hmomDiff250 -> Fit("gaus");


  c21 -> cd(1);
  hmomDiff275 ->GetXaxis() ->SetTitle("p_{difference}");
  hmomDiff275 ->GetYaxis() ->SetTitle("# of Track(positron)");
  hmomDiff275 -> Draw();
  hmomDiff275 -> Fit("gaus");



  // csv start 
  {
    double * xE;
    double * yE;
    ofstream csv("Efficiency.csv");
    int nEff = Eff -> GetN();
    csv  << "Efficiency" << endl;
    for(Int_t i = 0; i < nEff ; i++)
    {
      xE = Eff -> GetX();
      yE = Eff -> GetY();
      csv << xE[i] << "," << yE[i] << endl;
    }
    csv << "Efficiency, 200MeV/c < p_{initial < 275MeV/c}" << endl;
    csv << EffO200 << endl;
    csv << "# of positrons" << endl;

    csv << "~25"<< "," << "~50" << "," << "50~75" << ","<< "75~100" << ","<< "100~125" << ","<< "125~150" << ","<< "150~175" << ","<< "175~200" << ","<< "200~225" << ","<< "225~250" << ","<< "250~275" << ","<< "275~300" << endl;
    csv << Nt25 << Nt50 << "," << Nt75 << "," << Nt100 << "," << Nt125 << "," << Nt150 << "," << Nt175 << "," << Nt200 << "," << Nt225 << "," << Nt250 << "," <<  Nt275 << "," << Nt300 << endl;  
    csv.close();
  }
  // ////////////////////////////////////////////////CSV is done/////////////////////////////////////////////////////////

  display->setOptions("ABDEFHMPT"); // G show geometry
  //if (matFX) display->setOptions("ABDEFGHMPT"); // G show geometry
  display->open();


}
