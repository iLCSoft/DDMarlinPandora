#ifndef DDSimpleMuonDigi_H
#define DDSimpleMuonDigi_H 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>

#include "CalorimeterHitType.h"


using namespace lcio ;
using namespace marlin ;


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
  
 protected:

  int _nRun ;
  int _nEvt ;

  IntVec _layersToKeepBarrelVec, _layersToKeepEndcapVec;
  std::vector<bool>  _useLayersBarrelVec, _useLayersEndcapVec;

  std::vector<std::string> _muonCollections;

  std::string _outputMuonCollection;
  std::string _outputRelCollection;

  std::string _cellIDLayerString ;

  float _thresholdMuon;
  float _calibrCoeffMuon;
  float _maxHitEnergyMuon;

  std::string _detectorNameBarrel;
  std::string _detectorNameEndcap;
  
  
} ;

#endif



