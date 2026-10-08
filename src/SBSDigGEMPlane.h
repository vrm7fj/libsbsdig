#ifndef SBSDIGGEMPLANE_H
#define SBSDIGGEMPLANE_H

#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <TROOT.h>
//#include "TH1D.h"
#include "TRandom3.h"
#include "TH2D.h"
#include "TString.h"
#include "TF1.h"

class TH2D;

struct Pedestal {
  double mean{};
  double rms{};
};

struct CommonMode {
  double mean{};
  double sigma{};
};

struct CMFit {
  double mean{}, sigma{}, alpha{}, n{}, lo{}, hi{};
};

struct ApvTS {
    int APV_val;
    int ts_val;

    bool operator<(const ApvTS& other) const {
        if (APV_val != other.APV_val)
            return APV_val < other.APV_val;
        return ts_val < other.ts_val;
    }
};

using pedmap = std::map<int, Pedestal>; // Map ped mean and rms value to each strip number.
using cmmap = std::map<int, CommonMode>; // Map CM mean and sigma values to each APV number.
using cmmeanmap = std::map<ApvTS, double>; // Map to hold true simulated CM mean value for a all the strips and TS.
using cmfitmap = std::map<int, CMFit>;

constexpr int fNChanAPV {128};

class SBSDigGEMPlane {
 public:
  SBSDigGEMPlane();
  SBSDigGEMPlane(short trackeruniqueid, short layer, short mod, int nstrips, int nsamples = 6, double thr = 100, double offset = 0, double roangle = 0);
  virtual ~SBSDigGEMPlane();
  void Clear();
  
  //void SetStripThreshold(double thr){Striphr = ;};

  Short_t Layer(){return fLayer;};
  Short_t Module(){return fModule;};
  Double_t dX(){return fdX;};
  Double_t Xoffset(){return fXoffset;};
  Double_t ROangle(){return fROangle;};
  Int_t GetNStrips(){return fNStrips;};
  Short_t GetADC(int strip, int samp){return  fStripADC[strip*fNSamples+samp];};
  Int_t GetADCSum(int strip){return fStripADCsum[strip];};
  Int_t GetGoodADCSum(int strip){return fStripMipADCsum[strip];};
  void SetADC(int strip, int samp, int adc){
    if(strip<fNStrips){
      fStripADCsum[strip]+= adc-fStripADC[strip*fNSamples+samp];
      fStripADC[strip*fNSamples+samp] = adc;
    }
  };
  void AddADC(int strip, int samp, int adc){
    if(strip<fNStrips){
      fStripADC[strip*fNSamples+samp]+=adc;
      fStripADCsum[strip]+=adc;
    }//else{
    //printf("strip = %d / %d, sample %d /%d \n", strip, fNStrips, samp, fNSamples);
    //}
  };
  void SetGoodADC(int strip, int samp, int adc){
    if(strip<fNStrips){
      fStripMipADCsum[strip]+= adc-fStripMipADC[strip*fNSamples+samp];
      fStripMipADC[strip*fNSamples+samp] = adc;
    }
  };
  void AddGoodADC(int strip, int samp, int adc){
    if(strip<fNStrips){
      fStripMipADC[strip*fNSamples+samp]+=adc;
      fStripMipADCsum[strip]+=adc;    
    }
  };
  Short_t GetGoodADC(int strip, int samp){ return fStripMipADC[strip*fNSamples+samp]; };

  void SetPlaneStripPed(const pedmap& PedestalMap){ // This function will set the pedestal mean and RMS for each strip in the GEM plane.
    fPedestalMap = PedestalMap;
  };

  void SetPlaneAPVCM(const cmmap& CommonmodeMap){  // This function will set CM mean and std.dev values for each APV in the GEM plane.
    fCommonmodeMap = CommonmodeMap;    
  };

  void SetPlaneAPVCMFit(const cmfitmap& fits);
  double SampleTrueCM(int apv, TRandom3* R);

  Double_t GetStripPedMean(int strip) {return fPedestalMap[strip].mean;};
  Double_t GetStripPedRMS(int strip) {return fPedestalMap[strip].rms;};
  Double_t GetAPVCMMean(int apv) {return fCommonmodeMap[apv].mean;};
  Double_t GetAPVCMSigma(int apv) {return fCommonmodeMap[apv].sigma;};

  void DoPedSub();
  double GetOnlineCommonMode(const std::array <int,fNChanAPV>&, int); // Danning method online version CM calculation.
  void ApplyOnlineCMCorr( const bool do_onlineCM_histos = false ); // Implement CM correction to digitized data.
  void WriteCMHistos(); // Write online CM histos.
  void ApplyOnlineZS(const double zs_thr_nsigma); // Implement zero suppression. Only to be done after ApplyOnlineCMCorr().
  //void InsertTrueCM( const int iAPV, const int isamp, const double true_cm_mean ) { fTrueCMmap.emplace(ApvTS{iAPV,isamp}, true_cm_mean); };
  void InsertTrueCM( const int iAPV, const int isamp, const double true_cm_mean ) {
    fTrueCMmap.emplace(ApvTS{iAPV,isamp}, true_cm_mean);
    fTrueCMbyAPV->Fill(iAPV, true_cm_mean);
  };
  void ClearTrueCMmap(){ fTrueCMmap.clear(); }; // Has to be cleared after each loop over the 'j' strips in the SBSDigGEMSimDig::Checkout() method.

 private:
  // ADC sampled value of strip array of each axis
  Short_t fLayer;
  Short_t fModule;
  Int_t fNStrips;
  Int_t fNSamples;
  Int_t fNAPVs; // Number of APV cards in the GEM plane.
  Double_t fStripThr;//threshold for ADC sum
  Int_t* fStripADCsum;
  Short_t* fStripADC;
  Short_t* fStripPedSubADC; // Pedestal subtracted ADC of the strip -> For online CM calcualtions.
  Short_t* fStripCMCorrADC; // Online Danning method CM corrected ADC of the strip -> For online CM calculations.
  Int_t* fStripMipADCsum;
  Short_t* fStripMipADC;    // ADC contribution from the MIP avalanche -> To be used as MC truth infomation.
  
  double fdX;
  double fXoffset;
  double fROangle;
  short fTrackerID;
  //ClassDef(SBSDigGEMPlane, 1)

  // Ped and CM values.
  //TRandom3* R = new TRandom3(0); //Random number generation for ped and CM testing purposes.
  // Double_t* fStripPedOffset;
  // Double_t* fStripPedRMS;
  // Double_t* fAPVCMMean;
  // Double_t* fAPVCMSigma;
  pedmap fPedestalMap;
  cmmap fCommonmodeMap;

  std::map <std::pair<int,int>, int> fOnlineCMbyAPVbySamp;
  double fCommonModeRange_nsigma = 5.0;
  double fCommonModeDanningMethod_NsigmaCut = 3.0; // As set in SBSGEMModule.cxx as the default.
  int fCommonModeMinStripsInRange = 25;

  // Define CM histogram(s) to validate 'online' CM calculation accuracy.
  TH2D* fCommonModeDiffOnlineMean;
  cmmeanmap fTrueCMmap; // std::map to hold the true simulated CM values.
  TH2D* fCommonModeDiffOnlineTrue;

  TH2D* fTrueCMbyAPV;    // generated (true) CM per APV
  TH2D* fOnlineCMbyAPV;  // online (emulated firmware) CM per APV

  // Crystal Ball fit per APV, for sampling the true CM.
  std::map<int, TF1*> fCMFitFunc;   
  
};
#endif
