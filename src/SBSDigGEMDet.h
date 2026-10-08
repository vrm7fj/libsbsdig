#ifndef SBSDIGGEMDET_H
#define SBSDIGGEMDET_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include "TObjString.h"
#include "SBSDigGEMPlane.h"

using trackerpedmap = std::map<int, pedmap>;
using trackercmmap = std::map<int, cmmap>;
using trackercmfitmap = std::map<int, cmfitmap>;  

//________________________________
class SBSDigGEMDet {
 public:
  SBSDigGEMDet();
  SBSDigGEMDet(UShort_t uinqueid, UInt_t nplanes, int* layer, int* nstrips, double* offset, double* roangle, int nsamp, double zsup_thr, bool do_ped_cm = false, std::string pedfile = "ped.dat", std::string cmfile = "cm.dat", std::string cmfitfile = "");
  virtual ~SBSDigGEMDet();
  void Clear();

  struct gemhit{
    int source;
    int module;
    double edep;
    //double tmin;
    //double tmax;
    double t;
    double xin;
    double yin;
    double zin;
    double xout;
    double yout;
    double zout;
    int mid;
  };

  std::vector<gemhit> fGEMhits;
  
  //private:
  UShort_t fUniqueID;
  UInt_t fNPlanes;
  Double_t fGateWidth;
  std::vector<Double_t> fZLayer;
  //std::map<int, SBSDigGEMPlane> GEMPlanes;
  std::vector<SBSDigGEMPlane> GEMPlanes;

  //containers to store ped and cm values for the tracker.
  trackerpedmap fTrackerPedMap;
  trackercmmap fTrackerCommonModeMap;
  trackercmfitmap fTrackerCMFitMap;                   
};

#endif
