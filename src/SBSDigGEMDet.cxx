#include "SBSDigGEMDet.h"
#include "TMath.h"
#define DEBUG 0

using namespace std;

SBSDigGEMDet::SBSDigGEMDet()
{
}

SBSDigGEMDet::SBSDigGEMDet(UShort_t uniqueid, UInt_t nplanes, int* layer, int* nstrips, double* offset, double* roangle, int nsamp, double zsup_thr, bool do_ped_cm, std::string pedfile, std::string cmfile, std::string cmfitfile):
  fUniqueID(uniqueid), fNPlanes(nplanes)
{

  if (do_ped_cm){

     const char* env_libsbsdig = std::getenv("LIBSBSDIG");

    //Let's read in pedestals first.
    ifstream in_ped(std::string(env_libsbsdig)+"/db/"+pedfile);

    if(!in_ped.is_open()){
      cout << "GEM pedestal file " << pedfile.c_str() << " does not exist!!!" << endl;
      exit(-1);
    }

    TString currentline;
    int plane_num = -1; 

    while ( currentline.ReadLine(in_ped) ){
      if ( !currentline.BeginsWith("#") ) {
        
        if ( currentline.BeginsWith("GEM_plane#") ) {
          TObjArray *currentline_tokens = currentline.Tokenize(" ");
          if( currentline_tokens->GetEntries() >= 2 ){
           TString sval = ( (TObjString*) (*currentline_tokens)[1] )->GetString();
           plane_num = sval.Atoi();
         }

         delete currentline_tokens;
        }
        else {
          TObjArray *currentline_tokens = currentline.Tokenize(" ");
          if( currentline_tokens->GetEntries() >= 3 ){
            int this_stripnum = (((TObjString*)(*currentline_tokens)[0])->GetString()).Atoi();
            double this_strippedmean = (((TObjString*)(*currentline_tokens)[1])->GetString()).Atof();
            double this_strippedrms = (((TObjString*)(*currentline_tokens)[2])->GetString()).Atof();
            fTrackerPedMap[plane_num][this_stripnum] = Pedestal{this_strippedmean, this_strippedrms};
          }

          delete currentline_tokens;
        }
      }
    }

    if (fTrackerPedMap.size() != nplanes){
      cout << "The number of GEM planes defined in the GEM pedestal file " << pedfile.c_str() << " is inaccurate!!!" << endl;
      exit(-1);
    }

    for ( const auto& [plane_num, plane_pedmap] : fTrackerPedMap ){

      if ( plane_pedmap.size() != nstrips[plane_num] ){
        cout << "The number of strips defined in the GEM pedestal file " << pedfile.c_str() << " for plane number " << plane_num << " is incorrect!!!" << endl;
        exit(-1);
      }
    }

    //Now read-in cm.
    ifstream in_cm(std::string(env_libsbsdig)+"/db/"+cmfile);

    if(!in_cm.is_open()){
      cout << "GEM common-mode file: " << cmfile.c_str() << " does not exist!!!" << endl;
      exit(-1);
    }

    plane_num = -1;

    while ( currentline.ReadLine(in_cm) ){
      if ( !currentline.BeginsWith("#") ) {
      
        if ( currentline.BeginsWith("GEM_plane#") ) {
          TObjArray *currentline_tokens = currentline.Tokenize(" ");
          if( currentline_tokens->GetEntries() >= 2 ){
           TString sval = ( (TObjString*) (*currentline_tokens)[1] )->GetString();
           plane_num = sval.Atoi();
         }

          delete currentline_tokens;
        }
        else {
          TObjArray *currentline_tokens = currentline.Tokenize(" ");
          if( currentline_tokens->GetEntries() >= 3 ){
            int this_apvnum = (((TObjString*)(*currentline_tokens)[0])->GetString()).Atoi();
            double this_apvcmmean = (((TObjString*)(*currentline_tokens)[1])->GetString()).Atof();
            double this_apvcmsigma = (((TObjString*)(*currentline_tokens)[2])->GetString()).Atof();            
            fTrackerCommonModeMap[plane_num][this_apvnum] = CommonMode{this_apvcmmean, this_apvcmsigma};
          }

          delete currentline_tokens;
        }
      }
    }
    // ---------------------------------------------------------------------------------------------
    // Optional: Crystal Ball fits of the beam CM per APV (from fitCMCrystalBall.C).
    // APVs missing here, or with status != 0, keep the Gaussian from cmfile.
    if ( !cmfitfile.empty() ) {
      ifstream in_cmfit(std::string(env_libsbsdig)+"/db/"+cmfitfile);
      if ( !in_cmfit.is_open() ) {
        cout << "GEM CM fit file: " << cmfitfile.c_str() << " does not exist!!!" << endl;
        exit(-1);
      }
      while ( currentline.ReadLine(in_cmfit) ) {
        if ( currentline.BeginsWith("#") ) continue;
        TObjArray *tok = currentline.Tokenize(" ");
        if ( tok->GetEntries() >= 11 ) {
          auto val = [&](int k){ return ((TObjString*)(*tok)[k])->GetString(); };
          int plane  = val(1).Atoi();
          int apv    = val(2).Atoi();
          int status = val(10).Atoi();
          if ( status == 0 ) {
            fTrackerCMFitMap[plane][apv] = CMFit{ val(6).Atof(), val(7).Atof(),   // mean, sigma
                                                  val(8).Atof(), val(9).Atof(),   // alpha, n
                                                  val(3).Atof(), val(4).Atof() }; // fit_lo, fit_hi
          }
        }
        delete tok;
      }
      for ( const auto& [p, m] : fTrackerCMFitMap )
        cout << "GEM plane " << p << ": beam CM from Crystal Ball fits for " << m.size() << " APVs" << endl;
    }
    //-----------------------------------------------------------------------------------------------

    if (fTrackerCommonModeMap.size() != nplanes){
      cout << "The number of GEM planes defined in the GEM common-mode file " << cmfile.c_str() << " is inaccurate!!!" << endl;
      exit(-1);
    }

    for ( const auto& [plane_num, plane_cmmap] : fTrackerCommonModeMap ){

      if ( plane_cmmap.size() != (nstrips[plane_num]+128-1)/128 ){
        cout << "The number of APVs defined in the GEM common-mode file " << cmfile.c_str() << " for plane number " << plane_num << " is incorrect!!!" << endl;
        exit(-1);
      }
    }
  }// end do_ped_cm

  //for(uint i = 0; i<fNPlanes; i++)GEMPlanes[i] = SBSDigGEMPlane(nstrips[i], nsamp, zsup_thr);
  for(uint i = 0; i<fNPlanes; i++){
    cout << i << " " << layer[i] << " " << nstrips[i] << " " << offset[i] << " " << roangle[i] << endl; 
    GEMPlanes.push_back(SBSDigGEMPlane(fUniqueID, layer[i], i/*/2*/, nstrips[i], nsamp, zsup_thr, offset[i], roangle[i]));

    if (do_ped_cm){
      GEMPlanes.at(i).SetPlaneStripPed(fTrackerPedMap[i]);
      GEMPlanes.at(i).SetPlaneAPVCM(fTrackerCommonModeMap[i]);
      
    if ( fTrackerCMFitMap.count(i) ) GEMPlanes.at(i).SetPlaneAPVCMFit(fTrackerCMFitMap[i]);
    }
  }
}

SBSDigGEMDet::~SBSDigGEMDet()
{
  
}

void SBSDigGEMDet::Clear()
{
  for(uint i = 0; i<fNPlanes; i++)GEMPlanes[i].Clear();
  fGEMhits.clear();
}
