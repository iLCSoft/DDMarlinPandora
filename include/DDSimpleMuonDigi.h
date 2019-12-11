#ifndef DDSimpleMuonDigi_H
#define DDSimpleMuonDigi_H 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>

#include "CalorimeterHitType.h"


using namespace lcio ;
using namespace marlin ;

namespace EVENT {
  class SimCalorimeterHit ;
}


/** === DDSimpleMuonDigi Processor === <br>
 *  Simple calorimeter digitizer for the muon detectors.
 *  Converts SimCalorimeterHit collections to one 
 *  CalorimeterHit collection applying a threshold and an calibration constant...
 * 
 *  @version $Id$
 */
class DDSimpleMuonDigi : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new DDSimpleMuonDigi ; }
  
  
  DDSimpleMuonDigi() ;
  
  virtual void init() ;
  
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  virtual void end() ;

  bool useLayer(CHT::Layout caloLayout, unsigned int layer) ;
  float computeHitTime( const EVENT::SimCalorimeterHit *h ) const ;
  
 protected:

  int _nRun = 0;
  int _nEvt = 0;

  IntVec _layersToKeepBarrelVec{}, _layersToKeepEndcapVec{};
  std::vector<bool>  _useLayersBarrelVec{}, _useLayersEndcapVec{};

  std::vector<std::string> _muonCollections{};

  std::string _outputMuonCollection = "";
  std::string _outputRelCollection = "";

  std::string _cellIDLayerString = "layer";

  float _thresholdMuon = 0.025;
  float _timeThresholdMuon = _thresholdMuon ;
  float _calibrCoeffMuon = 120000;
  float _maxHitEnergyMuon = 2.0;

  std::string _detectorNameBarrel = "YokeBarrel";
  std::string _detectorNameEndcap = "YokeEndcap";
  
  
} ;

#endif



