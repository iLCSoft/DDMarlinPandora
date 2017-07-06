#ifndef DIGITIZER_DDCCALODIGI_H
#define DIGITIZER_DDCCALODIGI_H 1

#include "marlin/Processor.h"
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCFlagImpl.h>
#include "CalorimeterHitType.h"
#include "lcio.h"
#include <string>
#include <vector>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "DDScintillatorPpdDigi.h"
#include "CLHEP/Random/MTwistEngine.h"

using namespace lcio ;
using namespace marlin ;

const int MAX_LAYERS = 200;
const int MAX_STAVES =  16;

/** === DDCaloDigi Processor === <br>
 *  Simple calorimeter digitizer Processor. <br>
 *  Ported from ILDCaloDigi to use DD4hep
 *  Takes SimCalorimeterHit Collections and <br>
 *  produces CalorimeterHit Collections. <br>
 *  Simulated energy depositions in active <br>
 *  layers of calorimeters are <br>
 *  converted into physical energy. This is done <br>
 *  taking into account sampling fractions of <br>
 *  ECAL and HCAL. <br>
 *  User has to specify ECAL and HCAL SimCalorimeterHit <br>
 *  collections with processor parameters <br>
 *  HCALCollections and ECALCollections. <br>
 *  The names of the output CalorimeterHit Collections <br>
 *  are specified with processor parameters <br>
 *  ECALOutputCollection and HCALOutputCollection. <br>
 *  Conversion factors for ECAL and HCAL <br>
 *  are specified via processor parameters  <br>
 *  CalibrECAL and CalibrHCAL. <br>
 *  It should be noted that ECAL and HCAL may consist <br>
 *  of several sections with different sampling fractions. <br>
 *  To handle this situation, calibration coefficients for <br>
 *  ECAL and HCAL are passed as arrays of floats with each element <br>
 *  in this array corresponding to certain section with <br>
 *  a given sampling fraction. <br>
 *  List of layer numbers terminating each section are given through <br>
 *  processor parameters ECALLayers and HCALLayers <br>
 *  There is an option to perform digitization of <br> 
 *  both ECAL and HCAL in a digital mode. <br>
 *  Digital mode is activated by  <br>
 *  setting processor parameters <br>
 *  IfDigitalEcal / IfDigitalHcal to 1. <br>
 *  In this case CalibrECAL / CalibrHCAL will  <br>
 *  convert the number of hits into physical energy. <br>
 *  Thresholds on hit energies in ECAL and HCAL <br>
 *  are set with processor parameters <br>
 *  ECALThreshold and HCALThreshold.  <br>
 *  Relations between CalorimeterHits and SimCalorimeterHits <br>
 *  are held in the corresponding relation collection. <br>
 *  The name of this relation collection is specified <br>
 *  via processor parameter RelationOutputCollection. <br> 
 *  <h4>Input collections and prerequisites</h4>
 *  SimCalorimeterHit collections <br>
 *  <h4>Output</h4>
 *  CalorimeterHit collections for ECal and HCal. <br>
 *  Collection of relations <br>
 *  between CalorimeterHits and SimCalorimeterHits. <br> 
 *  For ECal Calorimeter hits the variable type is set to 0, <br>
 *  whereas for HCal Calorimeter hits the type is set to 1 <br>
 *  @author A. Raspereza (DESY) <br>
 *  @author M. Thomson (DESY) <br>
 *  @version $Id$ <br>
 */
class DDCaloDigi : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new DDCaloDigi ; }
  
  
  DDCaloDigi() ;
  DDCaloDigi(const DDCaloDigi&) = delete;
  DDCaloDigi& operator=(const DDCaloDigi&) = delete;
  
  virtual void init() ;
  
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  virtual void processEvent( LCEvent * evt ) ; 
   
  virtual void check( LCEvent * evt ) ; 
  
  virtual void end() ;

  virtual void fillECALGaps() ;
  
  float digitalHcalCalibCoeff(CHT::Layout,float energy );

  float analogueHcalCalibCoeff(CHT::Layout, int layer );

  float digitalEcalCalibCoeff(int layer );

  float analogueEcalCalibCoeff(int layer );

 protected:

  float ecalEnergyDigi(float energy, int id0, int id1);
  float ahcalEnergyDigi(float energy, int id0, int id1);

  float siliconDigi(float energy);
  float scintillatorDigi(float energy, bool isEcal);
  LCCollection* combineVirtualStripCells(LCCollection* col, bool isBarrel, int orientation );

  int getNumberOfVirtualCells();
  std::vector < std::pair <int, int> > & getLayerConfig();
  void checkConsistency(std::string colName, int layer);
  std::pair < int, int > getLayerProperties( std::string colName, int layer );
  int getStripOrientationFromColName( std::string colName );


  int _nRun = 0;
  int _nEvt = 0;
  
  LCFlagImpl _flag{};

  std::vector<std::string> _ecalCollections{};
  std::vector<std::string> _hcalCollections{};
  std::vector<std::string> _outputEcalCollections{};
  std::vector<std::string> _outputHcalCollections{};

  std::string _outputRelCollection = "";

  float _thresholdEcal = 5.0e-5;
  std::string _unitThresholdEcal = "GeV";
  std::vector<float> _thresholdHcal{};
  std::string _unitThresholdHcal = "GeV";

  int _digitalEcal = 0;
  int _mapsEcalCorrection = 0;
  int _digitalHcal = 0;

  //bool _ECAL_stripHits;

  std::vector<float> _calibrCoeffEcal{};
  std::vector<float> _calibrCoeffHcalBarrel{};
  std::vector<float> _calibrCoeffHcalEndCap{};
  std::vector<float> _calibrCoeffHcalOther{};

  std::vector<int> _ecalLayers{};
  std::vector<int> _hcalLayers{};

  int _ecalGapCorrection = 1;
  float _ecalGapCorrectionFactor = 1;
  float _ecalModuleGapCorrectionFactor = 0.5;
  float _ecalEndcapCorrectionFactor = 1.025;
  float _hcalEndcapCorrectionFactor = 1.025;
  int   _hcalGapCorrection = 1;
  float _hcalModuleGapCorrectionFactor = 0.5;

  std::vector<CalorimeterHitImpl*> _calHitsByStaveLayer[MAX_STAVES][MAX_LAYERS];
  std::vector<int> _calHitsByStaveLayerModule[MAX_STAVES][MAX_LAYERS];

  float _zOfEcalEndcap = 0.0;
  float _barrelPixelSizeT[MAX_LAYERS];
  float _barrelPixelSizeZ[MAX_LAYERS];
  float _endcapPixelSizeX[MAX_LAYERS];
  float _endcapPixelSizeY[MAX_LAYERS];
  float _barrelStaveDir[MAX_STAVES][2];
  
  int   _histograms = 0;

  // timing
  int   _useEcalTiming = 0;
  int   _ecalCorrectTimesForPropagation = 0;
  float _ecalTimeWindowMin = -10.0;
  float _ecalBarrelTimeWindowMax = 100.0;
  float _ecalEndcapTimeWindowMax = 100.0;
  float _ecalDeltaTimeHitResolution = 10.0;
  float _ecalTimeResolution = 10.0;
  bool  _ecalSimpleTimingCut = true;

  int   _useHcalTiming = 1;
  int   _hcalCorrectTimesForPropagation = 0;
  float _hcalTimeWindowMin = -10.0;
  float _hcalBarrelTimeWindowMax = 100.0;
  float _hcalEndcapTimeWindowMax = 100.0;
  float _hcalDeltaTimeHitResolution = 10.0;
  float _hcalTimeResolution = 10.0;
  bool  _hcalSimpleTimingCut = true;
  
  std::unique_ptr<DDScintillatorPpdDigi> _scEcalDigi{};
  std::unique_ptr<DDScintillatorPpdDigi> _scHcalDigi{};


  // parameters for extra ECAL digitization effects
  float _calibEcalMip = 1.0e-4;       // MIP calibration factor
  int   _applyEcalDigi = 0;           // which realistic calib to apply
  float _ecal_PPD_pe_per_mip = 7;     // # photoelectrons/MIP for MPPC
  int   _ecal_PPD_n_pixels = 10000;   // # pixels in MPPC
  float _ehEnergy = 3.6;              // energy to create e-h pair in silicon
  float _ecal_misCalibNpix = 0.05;    // miscalibration of # MPPC pixels

  float _misCalibEcal_uncorrel = 0.0; // general ECAL miscalibration (uncorrelated between channels)
  bool  _misCalibEcal_uncorrel_keep = false;// if true, use the same ECAL cell miscalibs in each event (requires more memory)
  float _misCalibEcal_correl = 0.0;     // general ECAL miscalibration (100% uncorrelated between channels)

  float _deadCellFractionEcal = 0.0;  // fraction of random dead channels
  bool  _deadCellEcal_keep = false;   // keep same cells dead between events?

  float _strip_abs_length = 1000000;  // absorption length along strip for non-uniformity modeling
  float _ecal_pixSpread = 0.05;       // relative spread of MPPC pixel signal
  float _ecal_elec_noise = 0;         // electronics noise (as fraction of MIP)
  float _ecalMaxDynMip = 2500;        // electronics dynamic range (in terms of MIPs)
  int _ecalStrip_default_nVirt = 9;   // # virtual cells used in Mokka simulation of strips (if available, this is taken from gear file)
  std::string _ecal_deafult_layer_config ="000000000000000";// ECAL layer configuration (if available, this is taken from gear file)

  // parameters for extra AHCAL digitization effects
  float _calibHcalMip = 1.0e-4;       // MIP calibration factor
  int   _applyHcalDigi = 0;           // which realistic calib to apply
  float _hcal_PPD_pe_per_mip = 10;    // # photoelectrons/MIP for MPPC
  int   _hcal_PPD_n_pixels= 400;      // # pixels in MPPC
  float _hcal_misCalibNpix = 0.05;    // miscalibration of # MPPC pixels

  float _misCalibHcal_uncorrel = 0.0; // general ECAL miscalibration (uncorrelated between channels)
  bool  _misCalibHcal_uncorrel_keep = false; // if true, use the same AHCAL cell miscalibs in each event (requires more memory)
  float _misCalibHcal_correl = 0.0;   // general ECAL miscalibration (100% uncorrelated between channels)

  float _deadCellFractionHcal = 0.0;  // fraction of random dead channels
  bool  _deadCellHcal_keep = false;   // keep same cells dead between events?
  float _hcal_pixSpread = 0.0;        // relative spread of MPPC pixel signal
  float _hcal_elec_noise = 0.0;       // electronics noise (as fraction of MIP)
  float _hcalMaxDynMip = 200;         // electronics dynamic range (in terms of MIPs)



  // internal variables
  std::vector < std::pair <int, int> > _layerTypes {};
  int _strip_virt_cells = 999;
  int _countWarnings = 0;
  std::string _ecalLayout = "";

  float _event_correl_miscalib_ecal = 0.0;
  float _event_correl_miscalib_hcal = 0.0;
  
  CLHEP::MTwistEngine *_randomEngineDeadCellEcal = NULL;
  CLHEP::MTwistEngine *_randomEngineDeadCellHcal = NULL;

  std::map < std::pair <int, int> , float > _ECAL_cell_miscalibs{};
  std::map < std::pair <int, int> , bool > _ECAL_cell_dead{};
  std::map < std::pair <int, int> , float > _HCAL_cell_miscalibs{};
  std::map < std::pair <int, int> , bool > _HCAL_cell_dead{};

  enum {
    SQUARE,
    STRIP_ALIGN_ALONG_SLAB,
    STRIP_ALIGN_ACROSS_SLAB,
    SIECAL=0,
    SCECAL
  };

  
  TH1F* fEcal = NULL;
  TH1F* fHcal = NULL;
  TH1F* fEcalC = NULL;
  TH1F* fHcalC = NULL;
  TH1F* fEcalC1 = NULL;
  TH1F* fHcalC1 = NULL;
  TH1F* fEcalC2 = NULL;
  TH1F* fHcalC2 = NULL;
  TH2F* fHcalCvsE = NULL;
  TH2F* fHcalLayer1 = NULL;
  TH2F* fHcalLayer11 = NULL;
  TH2F* fHcalLayer21 = NULL;
  TH2F* fHcalLayer31 = NULL;
  TH2F* fHcalLayer41 = NULL;
  TH2F* fHcalLayer51 = NULL;
  TH2F* fHcalLayer61 = NULL;
  TH2F* fHcalLayer71 = NULL;
  TH1F* fHcalRLayer1 = NULL;
  TH1F* fHcalRLayer11 = NULL;
  TH1F* fHcalRLayer21 = NULL;
  TH1F* fHcalRLayer31 = NULL;
  TH1F* fHcalRLayer41 = NULL;
  TH1F* fHcalRLayer51 = NULL;
  TH1F* fHcalRLayer61 = NULL;
  TH1F* fHcalRLayer71 = NULL;
  TH1F* fHcalRLayerNorm = NULL;

  TH1F* fEcalRLayerNorm = NULL;
  TH2F* fEcalLayer1 = NULL;
  TH2F* fEcalLayer11 = NULL;
  TH2F* fEcalLayer21 = NULL;
  TH1F* fEcalRLayer1 = NULL;
  TH1F* fEcalRLayer11 = NULL;
  TH1F* fEcalRLayer21 = NULL;

} ;

#endif



