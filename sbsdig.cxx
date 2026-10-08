//includes: standard
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>

//includes: root
#include <TROOT.h>
#include "TString.h"
#include "TObjString.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TCut.h"
#include "TEventList.h"
#include "TMath.h"
#include "TRandom3.h"

//includes: specific
#include "G4SBSRunData.hh"
#include "g4sbs_types.h"
#include "g4sbs_tree.h"
#include "SBSDigAuxi.h"
#include "SBSDigPMTDet.h"
#include "SBSDigGEMDet.h"
#include "SBSDigGEMSimDig.h"
#include "SBSDigBkgdGen.h"

#ifdef __APPLE__
#include "unistd.h"
#endif

#include <chrono>
/*
//Defining here the parameters for the new detectors.
//TODO: write a list of parameters that are not "frozen" (e.g. gain, pedestal parameters, etc...) and switch them into databases...
#define NPlanes_bbgem 32 // modules...
#define NChan_bbps 52
#define NChan_bbsh 189
#define NChan_bbhodo 180 
#define NChan_grinch 510 
#define NChan_hcal 288

#define TriggerJitter 3.0 //ns
#define ADCbits 12 
#define gatewidth_PMT 100 //ns
#define gatewidth_GEM 400 //ns

#define FADC_sampsize 4.0 //ns

//DB???
#define sigmapulse_bbpsSH 3.0 //ns / 1.2 of 
#define sigmapulse_bbhodo 1.6 //ns
#define sigmapulse_grinch 3.75 //ns

#define gain_bbps 2.e6
#define ped_bbps 600.0 // ADC channel
#define pedsigma_bbps 3.0 // ADC channel
#define trigoffset_bbps 18.2 //ns
#define ADCconv_bbps 50 //fC/ch

#define gain_bbsh 7.5e5
#define ped_bbsh 500.0 // ADC channel
#define pedsigma_bbsh 4.5 // ADC channel
#define trigoffset_bbsh 18.5 //ns
#define ADCconv_bbsh 50 //fC/ch

#define gain_grinch 7.0e6
#define ped_grinch 0.0 
#define pedsigma_grinch 0.0
#define trigoffset_grinch 15.3 //ns
#define threshold_grinch 3.e-3 //V
#define ADCconv_grinch 100 //fC/ch
#define TDCconv_grinch 1.0 //ns/channel
#define TDCbits_grinch 16 //ns/channel

#define gain_bbhodo 1.0e5
#define ped_bbhodo 0.0 
#define pedsigma_bbhodo 0.0
#define trigoffset_bbhodo 18.6 //ns
#define threshold_bbhodo 3.e-3 //V
#define ADCconv_bbhodo 100 //fC/ch
#define TDCconv_bbhodo 0.1 //ns/channel
#define TDCbits_bbhodo 19 //ns/channel

#define gain_hcal 1.0e6
#define ped_hcal 0.0 
#define pedsigma_hcal 0.0
#define trigoffset_hcal 81.0 //ns
#define threshold_hcal 3.e-3 //V
#define ADCconv_hcal 1.0 //fC/ch //??
#define TDCconv_hcal 0.12 //ns/channel
#define TDCbits_hcal 16 //ns/channel
*/

using namespace std;
//____________________________________________________
int main(int argc, char** argv){
  // Step 0: read out arguments
  string db_file, inputsigfile, inputbkgdfile = "";//sources of files
  ULong64_t Nentries = -1;//number of events to process
  //UShort_t Nbkgd = 0;//number of background files to add to each event
  double BkgdTimeWindow = 0, LumiFrac = 0;
  bool pmtbkgddig = false;
      
  if(argc<3 || argc>4){
    cout << "*** Inadequate number of arguments! ***" << endl
	 << " Arguments: database (mandatory); " << endl
	 << "           list_of_sig_input_files (str, mandatory); " << endl
	 << "          nb_of_sig_evts_to_process (int, def=-1); " << endl;
      //<< "         bkgd_histo_input_file (str, def=''); " << endl
      // << "        bkgd_lumi_frac (double, def=0); " << endl;
    return(-1);
  }
  
  db_file = argv[1];
  cout << " database file " << db_file << endl;
  inputsigfile = argv[2];
  cout << " Signal input files from: " << inputsigfile << endl;
  if(argc>3)Nentries = atoi(argv[3]);
  cout << " Number of (signal) events to process = " << Nentries << endl;
  /*
  if(argc>5){
    inputbkgdfile = argv[4];
    cout << " Background histgrams from: " << inputbkgdfile << endl;
    LumiFrac = max(0., atof(argv[5]));
    cout << " Fraction of background to superimpose to signal = " << LumiFrac << endl;
  }
  */
  
  // ------------------- // dev notes // ------------------- //
  // First, we want to extend the input tree (for signal only!!!)
  // I guess in order to avoid adding extra layers of code, 
  // the tree extension might have to be coded in the custom tree class

  
  std::vector<SBSDigPMTDet*> PMTdetectors;
  std::vector<int> detmap;
  std::vector<SBSDigGEMDet*> GEMdetectors;
  std::vector<SBSDigGEMSimDig*> GEMsimDig;
  std::vector<int> gemdetmap;
  
  // Variable parameters. 
  // Can be configured with the database, but are provided with defaults.
  Int_t Rseed = 0;
  Double_t TriggerJitter = 3.0;
  
  std::vector<TString> detectors_list;
  
  const int nparam_pmtdet_adc = 13;
  const int nparam_pmtdet_fadc = 12;
  const int nparam_pmtdet_fadc_hitposresp = 13;
  const int nparam_pmtdet_fadc_leadglass = 13;
  const int nparam_gemdet = 12;
  
  int nparam_bbps_read = 0;
  Double_t refindex_bbps = 1.68;
  Int_t NChan_bbps = 52;
  Double_t gatewidth_bbps = 100.;
  std::vector<Double_t> gain_bbps;
  std::vector<Int_t> nstages_bbps;
  std::vector<Double_t> ped_bbps;//
  std::vector<Double_t> pedsigma_bbps;
  Double_t trigoffset_bbps = 18.2;//
  Double_t threshold_bbps = 3.e-3;
  Double_t ADCconv_bbps = 50.;
  Int_t ADCbits_bbps = 12;
  Double_t TDCconv_bbps = 0.0625;
  Int_t TDCbits_bbps = 15;
  Double_t sigmapulse_bbps = 3.0;
    
  int nparam_bbsh_read = 0;
  Double_t refindex_bbsh = 1.68;
  Int_t NChan_bbsh = 189;
  Double_t gatewidth_bbsh = 100.;
  std::vector<Double_t> gain_bbsh;
  std::vector<Int_t> nstages_bbsh;
  std::vector<Double_t> ped_bbsh;
  std::vector<Double_t> pedsigma_bbsh;
  Double_t trigoffset_bbsh = 18.5;
  Double_t threshold_bbsh = 3.e-3;
  Double_t ADCconv_bbsh = 50.;
  Int_t ADCbits_bbsh = 12;
  Double_t TDCconv_bbsh = 0.0625;
  Int_t TDCbits_bbsh = 15;
  Double_t sigmapulse_bbsh = 3.0;
  
  int nparam_grinch_read = 0;
  Int_t NChan_grinch = 510;
  Double_t gatewidth_grinch = 100.;
  std::vector<Double_t> gain_grinch;
  std::vector<Int_t> nstages_grinch;
  std::vector<Double_t> ped_grinch;
  std::vector<Double_t> pedsigma_grinch;
  Double_t trigoffset_grinch = 15.3;
  Double_t threshold_grinch = 3.e-3;
  Double_t ADCconv_grinch = 100;
  Int_t ADCbits_grinch = 12;
  Double_t TDCconv_grinch = 1.;
  Int_t TDCbits_grinch = 16;
  Double_t sigmapulse_grinch = 3.75;
 
  int nparam_bbhodo_read = 0;
  Int_t NChan_bbhodo = 180;
  Double_t gatewidth_bbhodo = 100.;
  std::vector<Double_t> gain_bbhodo;
  std::vector<Int_t> nstages_bbhodo;
  std::vector<Double_t> ped_bbhodo;
  std::vector<Double_t> pedsigma_bbhodo;
  Double_t trigoffset_bbhodo = 18.6;
  Double_t threshold_bbhodo = 3.e3;
  Double_t ADCconv_bbhodo = 100.;
  Int_t ADCbits_bbhodo = 12;
  Double_t TDCconv_bbhodo = 0.1;
  Int_t TDCbits_bbhodo = 19;
  Double_t sigmapulse_bbhodo = 1.6;

  int nparam_hcal_read = 0;
  Int_t NChan_hcal = 288;
  Double_t gatewidth_hcal = 80;
  std::vector<Double_t> gain_hcal;
  std::vector<Int_t> nstages_hcal;
  std::vector<Double_t> ped_hcal;
  std::vector<Double_t> pedsigma_hcal;
  Double_t trigoffset_hcal = 81.;
  Double_t threshold_hcal = 3.e-3;
  Double_t ADCconv_hcal = 1.;
  Double_t TDCconv_hcal = 0.12;
  Int_t TDCbits_hcal = 16;
  Int_t FADC_ADCbits = 12;
  Double_t FADC_sampsize = 4.0;
  Double_t sigmapulse_hcal = 20.0;
  Double_t hcalzhitresp_p[3] = {5.242, 11.39, 10.41};
  
  int nparam_bbgem_read = 0;
  Int_t NPlanes_bbgem = 32;// number of planes/modules/readout
  Double_t gatewidth_bbgem = 400.;
  Double_t ZsupThr_bbgem = 240.;
  Int_t Nlayers_bbgem = 5;
  std::vector<Double_t> bbgem_layer_z;
  Int_t* layer_bbgem;
  Int_t* nstrips_bbgem;
  Double_t* offset_bbgem;
  Double_t* RO_angle_bbgem;
  Double_t* triggeroffset_bbgem;
  Double_t* gain_bbgem;//one gain per module
  Double_t* commonmode_array_bbgem;
  UShort_t nAPV_bbgem = 0;
  
  int nparam_sbsgem_read = 0;
  Int_t NPlanes_sbsgem = 32;// number of planes/modules/readout
  Double_t gatewidth_sbsgem = 400.;
  Double_t ZsupThr_sbsgem = 240.;
  Int_t Nlayers_sbsgem = 5;
  std::vector<Double_t> sbsgem_layer_z;
  Int_t* layer_sbsgem;
  Int_t* nstrips_sbsgem;
  Double_t* offset_sbsgem;
  Double_t* RO_angle_sbsgem;
  Double_t* triggeroffset_sbsgem;
  Double_t* gain_sbsgem;//one gain per module
  Double_t* commonmode_array_sbsgem;
  UShort_t nAPV_sbsgem = 0;
  
  //GEP detectors: parameters with dummy values...
  int nparam_ecal_read = 0;
  Double_t refindex_ecal = 1.68;
  Int_t NChan_ecal = 189;
  Double_t gatewidth_ecal = 100.;
  std::vector<Double_t> gain_ecal;
  std::vector<Int_t> nstages_ecal;
  std::vector<Double_t> ped_ecal;
  std::vector<Double_t> pedsigma_ecal;
  Double_t trigoffset_ecal = 18.5;
  Double_t threshold_ecal = 3.e-3;
  Double_t ADCconv_ecal = 50.;
  Int_t ADCbits_ecal = 12;
  Double_t TDCconv_ecal = 0.0625;
  Int_t TDCbits_ecal = 15;
  Double_t sigmapulse_ecal = 3.0;

  int nparam_cdet_read = 0;
  Int_t NChan_cdet = 189;
  Double_t gatewidth_cdet = 100.;
  std::vector<Double_t> gain_cdet;
  std::vector<Int_t> nstages_cdet;
  std::vector<Double_t> ped_cdet;
  std::vector<Double_t> pedsigma_cdet;
  Double_t trigoffset_cdet = 18.5;
  Double_t threshold_cdet = 3.e-3;
  Double_t ADCconv_cdet = 50.;
  Int_t ADCbits_cdet = 12;
  Double_t TDCconv_cdet = 0.0625;
  Int_t TDCbits_cdet = 15;
  Double_t sigmapulse_cdet = 3.0;
  
  int nparam_ft_read = 0;
  Int_t NPlanes_ft = 36;// number of planes/modules/readout
  Double_t gatewidth_ft = 400.;
  Double_t ZsupThr_ft = 240.;
  Int_t Nlayers_ft = 6;
  std::vector<Double_t> ft_layer_z;
  Int_t* layer_ft;
  Int_t* nstrips_ft;
  Double_t* offset_ft;
  Double_t* RO_angle_ft;
  Double_t* triggeroffset_ft;
  Double_t* gain_ft;//one gain per module
  Double_t* commonmode_array_ft;
  UShort_t nAPV_ft = 0;
  bool do_pedcm_ft = false;
  std::string pedfile_ft = "";
  std::string cmfile_ft = "";
  std::string cmfitfile_ft = "";
  bool do_online_cm_ft = false;
  bool do_online_zs_ft = false;
  Double_t online_zs_thr_nsigma_ft = 3.0;
  bool make_online_cmplots_ft = false;

  int nparam_fpp1_read = 0;
  Int_t NPlanes_fpp1 = 40;// number of planes/modules/readout
  Double_t gatewidth_fpp1 = 400.;
  Double_t ZsupThr_fpp1 = 240.;
  Int_t Nlayers_fpp1 = 5;
  std::vector<Double_t> fpp1_layer_z;
  Int_t* layer_fpp1;
  Int_t* nstrips_fpp1;
  Double_t* offset_fpp1;
  Double_t* RO_angle_fpp1;
  Double_t* triggeroffset_fpp1;
  Double_t* gain_fpp1;//one gain per module
  Double_t* commonmode_array_fpp1;
  UShort_t nAPV_fpp1 = 0;
  
  int nparam_fpp2_read = 0;
  Int_t NPlanes_fpp2 = 40;// number of planes/modules/readout
  Double_t gatewidth_fpp2 = 400.;
  Double_t ZsupThr_fpp2 = 240.;
  Int_t Nlayers_fpp2 = 5;
  std::vector<Double_t> fpp2_layer_z;
  Int_t* layer_fpp2;
  Int_t* nstrips_fpp2;
  Double_t* offset_fpp2;
  Double_t* RO_angle_fpp2;
  Double_t* triggeroffset_fpp2;
  Double_t* gain_fpp2;//one gain per module
  Double_t* commonmode_array_fpp2;
  UShort_t nAPV_fpp2 = 0;

  int nparam_cepol_front_read = 0;
  Int_t NPlanes_cepol_front = 28;// number of planes/modules/readout
  Double_t gatewidth_cepol_front = 400.;
  Double_t ZsupThr_cepol_front = 240.;
  Int_t Nlayers_cepol_front = 4;
  std::vector<Double_t> cepol_front_layer_z;
  Int_t* layer_cepol_front;
  Int_t* nstrips_cepol_front;
  Double_t* offset_cepol_front;
  Double_t* RO_angle_cepol_front;
  Double_t* triggeroffset_cepol_front;
  Double_t* gain_cepol_front;//one gain per module
  Double_t* commonmode_array_cepol_front;
  UShort_t nAPV_cepol_front = 0;

  int nparam_cepol_rear_read = 0;
  Int_t NPlanes_cepol_rear = 32;// number of planes/modules/readout
  Double_t gatewidth_cepol_rear = 400.;
  Double_t ZsupThr_cepol_rear = 240.;
  Int_t Nlayers_cepol_rear = 4;
  std::vector<Double_t> cepol_rear_layer_z;
  Int_t* layer_cepol_rear;
  Int_t* nstrips_cepol_rear;
  Double_t* offset_cepol_rear;
  Double_t* RO_angle_cepol_rear;
  Double_t* triggeroffset_cepol_rear;
  Double_t* gain_cepol_rear;//one gain per module
  Double_t* commonmode_array_cepol_rear;
  UShort_t nAPV_cepol_rear = 0;
  
  int nparam_prpolbs_gem_read = 0;
  Int_t NPlanes_prpolbs_gem = 16;// number of planes/modules/readout
  Double_t gatewidth_prpolbs_gem = 400.;
  Double_t ZsupThr_prpolbs_gem = 240.;
  Int_t Nlayers_prpolbs_gem = 2;
  std::vector<Double_t> prpolbs_gem_layer_z;
  Int_t* layer_prpolbs_gem;
  Int_t* nstrips_prpolbs_gem;
  Double_t* offset_prpolbs_gem;
  Double_t* RO_angle_prpolbs_gem;
  Double_t* triggeroffset_prpolbs_gem;
  Double_t* gain_prpolbs_gem;//one gain per module
  Double_t* commonmode_array_prpolbs_gem;
  UShort_t nAPV_prpolbs_gem = 0;
  
  int nparam_prpolfs_gem_read = 0;
  Int_t NPlanes_prpolfs_gem = 16;// number of planes/modules/readout
  Double_t gatewidth_prpolfs_gem = 400.;
  Double_t ZsupThr_prpolfs_gem = 240.;
  Int_t Nlayers_prpolfs_gem = 2;
  std::vector<Double_t> prpolfs_gem_layer_z;
  Int_t* layer_prpolfs_gem;
  Int_t* nstrips_prpolfs_gem;
  Double_t* offset_prpolfs_gem;
  Double_t* RO_angle_prpolfs_gem;
  Double_t* triggeroffset_prpolfs_gem;
  Double_t* gain_prpolfs_gem;//one gain per module
  Double_t* commonmode_array_prpolfs_gem;
  UShort_t nAPV_prpolfs_gem = 0;
  

  // ** How to add a new subsystem **
  // Add param for new detectors there...
  //polscint_bs
  int nparam_prpolscint_bs_read = 0;
  Int_t NChan_polscint_bs = 48;
  Double_t gatewidth_polscint_bs = 30.;
  std::vector<Double_t> gain_polscint_bs;
  std::vector<Int_t> nstages_polscint_bs;
  std::vector<Double_t> ped_polscint_bs;
  std::vector<Double_t> pedsigma_polscint_bs;
  Double_t trigoffset_polscint_bs = 37.6;
  Double_t threshold_polscint_bs = 3.e3;
  Double_t ADCconv_polscint_bs = 100.;
  Int_t ADCbits_polscint_bs = 12;
  Double_t TDCconv_polscint_bs = 0.1;
  Int_t TDCbits_polscint_bs = 19;
  Double_t sigmapulse_polscint_bs = 1.6;
  
  //polscint_fs
  int nparam_prpolscint_fs_read = 0;
  Int_t NChan_polscint_fs = 48;
  Double_t gatewidth_polscint_fs = 30.;
  std::vector<Double_t> gain_polscint_fs;
  std::vector<Int_t> nstages_polscint_fs;
  std::vector<Double_t> ped_polscint_fs;
  std::vector<Double_t> pedsigma_polscint_fs;
  Double_t trigoffset_polscint_fs = 37.6;
  Double_t threshold_polscint_fs = 3.e3;
  Double_t ADCconv_polscint_fs = 100.;
  Int_t ADCbits_polscint_fs = 12;
  Double_t TDCconv_polscint_fs = 0.1;
  Int_t TDCbits_polscint_fs = 19;
  Double_t sigmapulse_polscint_fs = 1.6;
  
  //activeana
  int nparam_activeana_read = 0;
  Int_t NChan_activeana = 32;
  Double_t gatewidth_activeana = 30.;
  std::vector<Double_t> gain_activeana;
  std::vector<Int_t> nstages_activeana;
  std::vector<Double_t> ped_activeana;
  std::vector<Double_t> pedsigma_activeana;
  Double_t trigoffset_activeana = 37.6;
  Double_t threshold_activeana = 3.e3;
  Double_t ADCconv_activeana = 100.;
  Int_t ADCbits_activeana = 19;
  Double_t TDCconv_activeana = 0.1;
  Int_t TDCbits_activeana = 19;
  Double_t sigmapulse_activeana = 1.6;

  //-----------------------------
  //  Read database
  //-----------------------------
  cout << "read database: " << db_file.c_str() << endl;
  ifstream in_db(db_file.c_str());
  if(!in_db.is_open()){
    cout << "database " << db_file.c_str() << " does not exist!!!" << endl;
    exit(-1);
  }
  
  TString currentline;
  while( currentline.ReadLine(in_db) && !currentline.BeginsWith("endconfig")){
    if( !currentline.BeginsWith("#") ){
      Int_t ntokens = 0;
      std::unique_ptr<TObjArray> tokens( currentline.Tokenize(", \t") );
      if( !tokens->IsEmpty() ) {
	ntokens = tokens->GetLast()+1;
      }
      //TObjArray *tokens = currentline.Tokenize(" ");//vg: def lost => versions prior to 6.06; should be fixed! ??? 
      //int ntokens = tokens->GetEntries();
      
      if( ntokens >= 2 ){
	TString skey = ( (TObjString*) (*tokens)[0] )->GetString();
	
	if(skey=="Rseed"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  Rseed = stemp.Atoi();
	}
	
	if(skey=="TriggerJitter"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TriggerJitter = stemp.Atof();
	}
	
	if(skey=="detectors_list"){
	  for(int k = 1; k<ntokens; k++){
	    TString sdet = ( (TObjString*) (*tokens)[k] )->GetString();
	    detectors_list.push_back(sdet);
	  }
	}
	
	//BBPS
	if(skey=="refindex_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  refindex_bbps = stemp.Atof();
	  nparam_bbps_read++;
	}
	
	if(skey=="NChan_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_bbps = stemp.Atoi();
	  nparam_bbps_read++;
	}
	
	if(skey=="gatewidth_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_bbps = stemp.Atof();
	  nparam_bbps_read++;
	}
	
	if(skey=="gain_bbps"){
	  gain_bbps.resize(NChan_bbps);
	  if(ntokens==NChan_bbps+1){
	    for(int k = 0; k<NChan_bbps; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      gain_bbps[k] = stemp.Atof()*qe;
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_bbps << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_bbps; k++){
	      gain_bbps[k] = stemp.Atof()*qe;
	    }
	  }
	  nparam_bbps_read++;
	}

	if(skey=="nstages_bbps"){
	  nstages_bbps.resize(NChan_bbps);
	  if(ntokens==NChan_bbps+1){
	    for(int k = 0; k<NChan_bbps; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      nstages_bbps[k] = stemp.Atoi();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_bbps << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_bbps; k++){
	      nstages_bbps[k] = stemp.Atoi();
	    }
	  }
	  nparam_bbps_read++;
	}

	if(skey=="ped_bbps"){
	  // TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  // ped_bbps = stemp.Atof();
	  ped_bbps.resize(NChan_bbps);
	  if(ntokens==NChan_bbps+1){
	    for(int k = 0; k<NChan_bbps; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      ped_bbps[k] = stemp.Atof();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_bbps << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_bbps; k++){
	      ped_bbps[k] = stemp.Atof();
	    }
	  }
	  nparam_bbps_read++;
	}
	
	if(skey=="pedsigma_bbps"){
	  // TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  // pedsigma_bbps = stemp.Atof();
	  pedsigma_bbps.resize(NChan_bbps);
	  if(ntokens==NChan_bbps+1){
	    for(int k = 0; k<NChan_bbps; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      pedsigma_bbps[k] = stemp.Atof();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_bbps << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_bbps; k++){
	      pedsigma_bbps[k] = stemp.Atof();
	    }
	  }
	  nparam_bbps_read++;
	}
	
	if(skey=="threshold_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  threshold_bbps = stemp.Atof();
	  nparam_bbps_read++;
	}
	
	if(skey=="trigoffset_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_bbps = stemp.Atof();
	  nparam_bbps_read++;
	}
	
	if(skey=="ADCconv_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_bbps = stemp.Atof();
	  nparam_bbps_read++;
	}	

	if(skey=="ADCbits_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_bbps = stemp.Atoi();
	  nparam_bbps_read++;
	}	
	
	if(skey=="TDCconv_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCconv_bbps = stemp.Atof();
	  nparam_bbps_read++;
	}	

	if(skey=="TDCbits_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCbits_bbps = stemp.Atoi();
	  nparam_bbps_read++;
	}	
	
	if(skey=="sigmapulse_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_bbps = stemp.Atof();
	  nparam_bbps_read++;
	}
	
	//BBSH
	if(skey=="refindex_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  refindex_bbsh = stemp.Atof();
	  nparam_bbsh_read++;
	}
	
	if(skey=="NChan_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_bbsh = stemp.Atoi();
	  nparam_bbsh_read++;
	}
	
	if(skey=="gatewidth_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_bbsh = stemp.Atof();
	  nparam_bbsh_read++;
	}
	
	if(skey=="gain_bbsh"){
	  gain_bbsh.resize(NChan_bbsh);
	  if(ntokens==NChan_bbsh+1){
	    for(int k = 0; k<NChan_bbsh; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      gain_bbsh[k] = stemp.Atof()*qe;
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_bbsh << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_bbsh; k++){
	      gain_bbsh[k] = stemp.Atof()*qe;
	    }
	  }
	  nparam_bbsh_read++;
	}
	
	if(skey=="nstages_bbsh"){
	  nstages_bbsh.resize(NChan_bbsh);
	  if(ntokens==NChan_bbsh+1){
	    for(int k = 0; k<NChan_bbsh; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      nstages_bbsh[k] = stemp.Atoi();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_bbsh << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_bbsh; k++){
	      nstages_bbsh[k] = stemp.Atoi();
	    }
	  }
	  nparam_bbsh_read++;
	}
	
	if(skey=="ped_bbsh"){
	  // TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  // ped_bbsh = stemp.Atof();
	  ped_bbsh.resize(NChan_bbsh);
	  if(ntokens==NChan_bbsh+1){
	    for(int k = 0; k<NChan_bbsh; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      ped_bbsh[k] = stemp.Atof();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_bbsh << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_bbsh; k++){
	      ped_bbsh[k] = stemp.Atof();
	    }
	  }
	  nparam_bbsh_read++;
	}
	
	if(skey=="pedsigma_bbsh"){
	  // TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  // pedsigma_bbsh = stemp.Atof();
	  pedsigma_bbsh.resize(NChan_bbsh);
	  if(ntokens==NChan_bbsh+1){
	    for(int k = 0; k<NChan_bbsh; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      pedsigma_bbsh[k] = stemp.Atof();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_bbsh << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_bbsh; k++){
	      pedsigma_bbsh[k] = stemp.Atof();
	    }
	  }
	  nparam_bbsh_read++;
	}
	
	if(skey=="trigoffset_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_bbsh = stemp.Atof();
	  nparam_bbsh_read++;
	}
	
	if(skey=="threshold_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  threshold_bbsh = stemp.Atof();
	  nparam_bbsh_read++;
	}
	
	if(skey=="ADCconv_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_bbsh = stemp.Atof();
	  nparam_bbsh_read++;
	}	

	if(skey=="ADCbits_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_bbsh = stemp.Atoi();
	  nparam_bbsh_read++;
	}	
	
	if(skey=="TDCconv_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCconv_bbsh = stemp.Atof();
	  nparam_bbsh_read++;
	}	

	if(skey=="TDCbits_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCbits_bbsh = stemp.Atoi();
	  nparam_bbsh_read++;
	}
	
	if(skey=="sigmapulse_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_bbsh = stemp.Atof();
	  nparam_bbsh_read++;
	}

	//GRINCH
	if(skey=="NChan_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_grinch = stemp.Atoi();
	  nparam_grinch_read++;
	}
	
	if(skey=="gatewidth_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_grinch = stemp.Atof();
	  nparam_grinch_read++;
	}
	
	if(skey=="gain_grinch"){
	  gain_grinch.resize(NChan_grinch);
	  if(ntokens==NChan_grinch+1){
	    for(int k = 0; k<NChan_grinch; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      gain_grinch[k] = stemp.Atof()*qe;
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_grinch << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_grinch; k++){
	      gain_grinch[k] = stemp.Atof()*qe;
	    }
	  }
	  nparam_grinch_read++;
	}
	
	if(skey=="nstages_grinch"){
	  nstages_grinch.resize(NChan_grinch);
	  if(ntokens==NChan_grinch+1){
	    for(int k = 0; k<NChan_grinch; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      nstages_grinch[k] = stemp.Atoi();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_grinch << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_grinch; k++){
	      nstages_grinch[k] = stemp.Atoi();
	    }
	  }
	  nparam_grinch_read++;
	}

	if(skey=="ped_grinch"){
	  // TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  // ped_grinch = stemp.Atof();
	  ped_grinch.resize(NChan_grinch);
	  if(ntokens==NChan_grinch+1){
	    for(int k = 0; k<NChan_grinch; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      ped_grinch[k] = stemp.Atof();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_grinch << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_grinch; k++){
	      ped_grinch[k] = stemp.Atof();
	    }
	  }
	  nparam_grinch_read++;
	}
	
	if(skey=="pedsigma_grinch"){
	  // TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  // pedsigma_grinch = stemp.Atof();
	  pedsigma_grinch.resize(NChan_grinch);
	  if(ntokens==NChan_grinch+1){
	    for(int k = 0; k<NChan_grinch; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      pedsigma_grinch[k] = stemp.Atof();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_grinch << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_grinch; k++){
	      pedsigma_grinch[k] = stemp.Atof();
	    }
	  }
	  nparam_grinch_read++;
	}
	
	if(skey=="trigoffset_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_grinch = stemp.Atof();
	  nparam_grinch_read++;
	}
	
	if(skey=="threshold_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  threshold_grinch = stemp.Atof();
	  nparam_grinch_read++;
	}
	
	if(skey=="ADCconv_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_grinch = stemp.Atof();
	  nparam_grinch_read++;
	}	

	if(skey=="ADCbits_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_grinch = stemp.Atoi();
	  nparam_grinch_read++;
	}	
	
	if(skey=="TDCconv_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCconv_grinch = stemp.Atof();
	  nparam_grinch_read++;
	}	
	
	if(skey=="TDCbits_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCbits_grinch = stemp.Atoi();
	  nparam_grinch_read++;
	}	
	
	if(skey=="sigmapulse_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_grinch = stemp.Atof();
	  nparam_grinch_read++;
	}
	
	//BBHODO
	if(skey=="NChan_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_bbhodo = stemp.Atoi();
	  nparam_bbhodo_read++;
	}
	
	if(skey=="gatewidth_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_bbhodo = stemp.Atof();
	  nparam_bbhodo_read++;
	}
	
	if(skey=="gain_bbhodo"){
	  gain_bbhodo.resize(NChan_bbhodo);
	  if(ntokens==NChan_bbhodo+1){
	    for(int k = 0; k<NChan_bbhodo; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      gain_bbhodo[k] = stemp.Atof()*qe;
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_bbhodo << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_bbhodo; k++){
	      gain_bbhodo[k] = stemp.Atof()*qe;
	    }
	  }
	  nparam_bbhodo_read++;
	}
	
	if(skey=="nstages_bbhodo"){
	  nstages_bbhodo.resize(NChan_bbhodo);
	  if(ntokens==NChan_bbhodo+1){
	    for(int k = 0; k<NChan_bbhodo; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      nstages_bbhodo[k] = stemp.Atoi();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_bbhodo << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_bbhodo; k++){
	      nstages_bbhodo[k] = stemp.Atoi();
	    }
	  }
	  nparam_bbhodo_read++;
	}

	if(skey=="ped_bbhodo"){
	  // TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  // ped_bbhodo = stemp.Atof();
	  ped_bbhodo.resize(NChan_bbhodo);
	  if(ntokens==NChan_bbhodo+1){
	    for(int k = 0; k<NChan_bbhodo; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      ped_bbhodo[k] = stemp.Atof();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_bbhodo << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_bbhodo; k++){
	      ped_bbhodo[k] = stemp.Atof();
	    }
	  }
	  nparam_bbhodo_read++;
	}
	
	if(skey=="pedsigma_bbhodo"){
	  // TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  // pedsigma_bbhodo = stemp.Atof();
	  pedsigma_bbhodo.resize(NChan_bbhodo);
	  if(ntokens==NChan_bbhodo+1){
	    for(int k = 0; k<NChan_bbhodo; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      pedsigma_bbhodo[k] = stemp.Atof();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_bbhodo << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_bbhodo; k++){
	      pedsigma_bbhodo[k] = stemp.Atof();
	    }
	  }
	  nparam_bbhodo_read++;
	}
	
	if(skey=="trigoffset_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_bbhodo = stemp.Atof();
	  nparam_bbhodo_read++;
	}
	
	if(skey=="threshold_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  threshold_bbhodo = stemp.Atof();
	  nparam_bbhodo_read++;
	}
	
	if(skey=="ADCconv_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_bbhodo = stemp.Atof();
	  nparam_bbhodo_read++;
	}	

	if(skey=="ADCbits_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_bbhodo = stemp.Atoi();
	  nparam_bbhodo_read++;
	}	
	
	if(skey=="TDCconv_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCconv_bbhodo = stemp.Atof();
	  nparam_bbhodo_read++;
	}	
	
	if(skey=="TDCbits_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCbits_bbhodo = stemp.Atoi();
	  nparam_bbhodo_read++;
	}	
	
	if(skey=="sigmapulse_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_bbhodo = stemp.Atof();
	  nparam_bbhodo_read++;
	}
	
	//HCal
	if(skey=="NChan_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_hcal = stemp.Atoi();
	  nparam_hcal_read++;
	}
	
	if(skey=="gatewidth_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_hcal = stemp.Atof();
	  nparam_hcal_read++;
	}
	
	if(skey=="gain_hcal"){
	  gain_hcal.resize(NChan_hcal);
	  if(ntokens==NChan_hcal+1){
	    for(int k = 0; k<NChan_hcal; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      gain_hcal[k] = stemp.Atof()*qe;
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_hcal << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_hcal; k++){
	      gain_hcal[k] = stemp.Atof()*qe;
	    }
	  }
	  nparam_hcal_read++;
	}

	if(skey=="nstages_hcal"){
	  nstages_hcal.resize(NChan_hcal);
	  if(ntokens==NChan_hcal+1){
	    for(int k = 0; k<NChan_hcal; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      nstages_hcal[k] = stemp.Atoi();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_hcal << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_hcal; k++){
	      nstages_hcal[k] = stemp.Atoi();
	    }
	  }
	  nparam_hcal_read++;
	}

	if(skey=="ped_hcal"){
	  // TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  // ped_hcal = stemp.Atof();
	  ped_hcal.resize(NChan_hcal);
	  if(ntokens==NChan_hcal+1){
	    for(int k = 0; k<NChan_hcal; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      ped_hcal[k] = stemp.Atof();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_hcal << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_hcal; k++){
	      ped_hcal[k] = stemp.Atof();
	    }
	  }
	  nparam_hcal_read++;
	}
	
	if(skey=="pedsigma_hcal"){
	  // TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  // pedsigma_hcal = stemp.Atof();
	  pedsigma_hcal.resize(NChan_hcal);
	  if(ntokens==NChan_hcal+1){
	    for(int k = 0; k<NChan_hcal; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      pedsigma_hcal[k] = stemp.Atof();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_hcal << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_hcal; k++){
	      pedsigma_hcal[k] = stemp.Atof();
	    }
	  }
	  nparam_hcal_read++;
	}
	
	if(skey=="trigoffset_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_hcal = stemp.Atof();
	  nparam_hcal_read++;
	}
	
	if(skey=="threshold_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  threshold_hcal = stemp.Atof();
	  nparam_hcal_read++;
	}
	
	if(skey=="ADCconv_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_hcal = stemp.Atof();
	  nparam_hcal_read++;
	}	
	
	if(skey=="TDCconv_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCconv_hcal = stemp.Atof();
	  nparam_hcal_read++;
	}	
	
	if(skey=="TDCbits_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCbits_hcal = stemp.Atoi();
	  nparam_hcal_read++;
	}	
	
	if(skey=="sigmapulse_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_hcal = stemp.Atof();
	  nparam_hcal_read++;
	}
	
	if(skey=="zhitresp_hcal"){
	  if(ntokens==3+1){
	    for(int k = 0; k<3; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      hcalzhitresp_p[k] = stemp.Atof(); 
	    }
	    nparam_hcal_read++;
	  }
	}
	
	if(skey=="FADC_ADCbits"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  FADC_ADCbits = stemp.Atoi();
	}
	
	if(skey=="FADC_sampsize"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  FADC_sampsize = stemp.Atof();
	}

	//Ecal
	if(skey=="refindex_ecal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  refindex_ecal = stemp.Atof();
	  nparam_ecal_read++;
	}
	
	if(skey=="NChan_ecal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_ecal = stemp.Atoi();
	  nparam_ecal_read++;
	}
	
	if(skey=="gatewidth_ecal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_ecal = stemp.Atof();
	  nparam_ecal_read++;
	}
	
	if(skey=="gain_ecal"){
	  gain_ecal.resize(NChan_ecal);
	  if(ntokens==NChan_ecal+1){
	    for(int k = 0; k<NChan_ecal; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      gain_ecal[k] = stemp.Atof()*qe;
	    }
	  }else{
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_ecal; k++){
	      gain_ecal[k] = stemp.Atof()*qe;
	    }
	  }
	  nparam_ecal_read++;
	}

	if(skey=="nstages_ecal"){
	  nstages_ecal.resize(NChan_ecal);
	  if(ntokens==NChan_ecal+1){
	    for(int k = 0; k<NChan_ecal; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      nstages_ecal[k] = stemp.Atoi();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_ecal << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_ecal; k++){
	      nstages_ecal[k] = stemp.Atoi();
	    }
	  }
	  nparam_ecal_read++;
	}

	if(skey=="ped_ecal"){
	  // TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  // ped_ecal = stemp.Atof();
	  ped_ecal.resize(NChan_ecal);
	  if(ntokens==NChan_ecal+1){
	    for(int k = 0; k<NChan_ecal; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      ped_ecal[k] = stemp.Atof();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_ecal << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_ecal; k++){
	      ped_ecal[k] = stemp.Atof();
	    }
	  }
	  nparam_ecal_read++;
	}
	
	if(skey=="pedsigma_ecal"){
	  // TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  // pedsigma_ecal = stemp.Atof();
	  pedsigma_ecal.resize(NChan_ecal);
	  if(ntokens==NChan_ecal+1){
	    for(int k = 0; k<NChan_ecal; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      pedsigma_ecal[k] = stemp.Atof();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_ecal << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_ecal; k++){
	      pedsigma_ecal[k] = stemp.Atof();
	    }
	  }
	  nparam_ecal_read++;
	}
	
	if(skey=="trigoffset_ecal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_ecal = stemp.Atof();
	  nparam_ecal_read++;
	}
	
	if(skey=="threshold_ecal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  threshold_ecal = stemp.Atof();
	  nparam_ecal_read++;
	}
	
	if(skey=="ADCconv_ecal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_ecal = stemp.Atof();
	  nparam_ecal_read++;
	}	
	
	if(skey=="TDCconv_ecal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCconv_ecal = stemp.Atof();
	  nparam_ecal_read++;
	}	
	
	if(skey=="TDCbits_ecal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCbits_ecal = stemp.Atoi();
	  nparam_ecal_read++;
	}	
	
	if(skey=="sigmapulse_ecal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_ecal = stemp.Atof();
	  nparam_ecal_read++;
	}
	
	//CDET
	if(skey=="NChan_cdet"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_cdet = stemp.Atoi();
	  nparam_cdet_read++;
	}
	
	if(skey=="gatewidth_cdet"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_cdet = stemp.Atof();
	  nparam_cdet_read++;
	}
	
	if(skey=="gain_cdet"){
	  gain_cdet.resize(NChan_cdet);
	  if(ntokens==NChan_cdet+1){
	    for(int k = 0; k<NChan_cdet; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      gain_cdet[k] = stemp.Atof()*qe;
	    }
	  }else{
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_cdet; k++){
	      gain_cdet[k] = stemp.Atof()*qe;
	    }
	  }
	  nparam_cdet_read++;
	}
	
	if(skey=="nstages_cdet"){
	  nstages_cdet.resize(NChan_cdet);
	  if(ntokens==NChan_cdet+1){
	    for(int k = 0; k<NChan_cdet; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      nstages_cdet[k] = stemp.Atoi();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_cdet << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_cdet; k++){
	      nstages_cdet[k] = stemp.Atoi();
	    }
	  }
	  nparam_cdet_read++;
	}
	
	if(skey=="ped_cdet"){
	  // TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  // ped_cdet = stemp.Atof();
	  ped_cdet.resize(NChan_cdet);
	  if(ntokens==NChan_cdet+1){
	    for(int k = 0; k<NChan_cdet; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      ped_cdet[k] = stemp.Atof();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_cdet << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_cdet; k++){
	      ped_cdet[k] = stemp.Atof();
	    }
	  }
	  nparam_cdet_read++;
	}
	
	if(skey=="pedsigma_cdet"){
	  // TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  // pedsigma_cdet = stemp.Atof();
	  pedsigma_cdet.resize(NChan_cdet);
	  if(ntokens==NChan_cdet+1){
	    for(int k = 0; k<NChan_cdet; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      pedsigma_cdet[k] = stemp.Atof();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_cdet << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_cdet; k++){
	      pedsigma_cdet[k] = stemp.Atof();
	    }
	  }
	  nparam_cdet_read++;
	}
	
	if(skey=="trigoffset_cdet"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_cdet = stemp.Atof();
	  nparam_cdet_read++;
	}
	
	if(skey=="threshold_cdet"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  threshold_cdet = stemp.Atof();
	  nparam_cdet_read++;
	}
	
	if(skey=="ADCconv_cdet"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_cdet = stemp.Atof();
	  nparam_cdet_read++;
	}	

	if(skey=="ADCbits_cdet"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_cdet = stemp.Atoi();
	  nparam_cdet_read++;
	}	
	
	if(skey=="TDCconv_cdet"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCconv_cdet = stemp.Atof();
	  nparam_cdet_read++;
	}	
	
	if(skey=="TDCbits_cdet"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCbits_cdet = stemp.Atoi();
	  nparam_cdet_read++;
	}	
	
	if(skey=="sigmapulse_cdet"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_cdet = stemp.Atof();
	  nparam_cdet_read++;
	}
	
	
	// ** How to add a new subsystem **
	// Add reading of param from other detectors there...
	//GEn-RP Hodoscopes bs
	if(skey=="NChan_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_polscint_bs = stemp.Atoi();
	  nparam_prpolscint_bs_read++;
	}
	
	if(skey=="gatewidth_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_polscint_bs = stemp.Atof();
	  nparam_prpolscint_bs_read++;
	}
	
	if(skey=="gain_polscint_bs"){
	  gain_polscint_bs.resize(NChan_polscint_bs);
	  if(ntokens==NChan_polscint_bs+1){
	    for(int k = 0; k<NChan_polscint_bs; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      gain_polscint_bs[k] = stemp.Atof()*qe;
	    }
	  }else{
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_polscint_bs; k++){
	      gain_polscint_bs[k] = stemp.Atof()*qe;
	    }
	  }
	  nparam_prpolscint_bs_read++;
	}
	
	if(skey=="nstages_polscint_bs"){
	  nstages_polscint_bs.resize(NChan_polscint_bs);
	  if(ntokens==NChan_polscint_bs+1){
	    for(int k = 0; k<NChan_polscint_bs; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      nstages_polscint_bs[k] = stemp.Atoi();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_polscint_bs << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_polscint_bs; k++){
	      nstages_polscint_bs[k] = stemp.Atoi();
	    }
	  }
	  nparam_prpolscint_bs_read++;
	}
	
	if(skey=="ped_polscint_bs"){
	  // TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  // ped_polscint_bs = stemp.Atof();
	  ped_polscint_bs.resize(NChan_polscint_bs);
	  if(ntokens==NChan_polscint_bs+1){
	    for(int k = 0; k<NChan_polscint_bs; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      ped_polscint_bs[k] = stemp.Atof();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_polscint_bs << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_polscint_bs; k++){
	      ped_polscint_bs[k] = stemp.Atof();
	    }
	  }
	  nparam_prpolscint_bs_read++;
	}
	
	if(skey=="pedsigma_polscint_bs"){
	  // TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  // pedsigma_polscint_bs = stemp.Atof();
	  pedsigma_polscint_bs.resize(NChan_polscint_bs);
	  if(ntokens==NChan_polscint_bs+1){
	    for(int k = 0; k<NChan_polscint_bs; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      pedsigma_polscint_bs[k] = stemp.Atof();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_polscint_bs << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_polscint_bs; k++){
	      pedsigma_polscint_bs[k] = stemp.Atof();
	    }
	  }
	  nparam_prpolscint_bs_read++;
	}
	
	if(skey=="trigoffset_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_polscint_bs = stemp.Atof();
	  nparam_prpolscint_bs_read++;
	}
	
	if(skey=="threshold_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  threshold_polscint_bs = stemp.Atof();
	  nparam_prpolscint_bs_read++;
	}
	
	if(skey=="ADCconv_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_polscint_bs = stemp.Atof();
	  nparam_prpolscint_bs_read++;
	}	

	if(skey=="ADCbits_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_polscint_bs = stemp.Atoi();
	  nparam_prpolscint_bs_read++;
	}	
	
	if(skey=="TDCconv_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCconv_polscint_bs = stemp.Atof();
	  nparam_prpolscint_bs_read++;
	}	
	
	if(skey=="TDCbits_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCbits_polscint_bs = stemp.Atoi();
	  nparam_prpolscint_bs_read++;
	}	
	
	if(skey=="sigmapulse_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_polscint_bs = stemp.Atof();
	  nparam_prpolscint_bs_read++;
	}
	
	//GEn-RP Hodoscopes fs
	if(skey=="NChan_polscint_fs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_polscint_fs = stemp.Atoi();
	  nparam_prpolscint_fs_read++;
	}
	
	if(skey=="gatewidth_polscint_fs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_polscint_fs = stemp.Atof();
	  nparam_prpolscint_fs_read++;
	}
	
	if(skey=="gain_polscint_fs"){
	  gain_polscint_fs.resize(NChan_polscint_fs);
	  if(ntokens==NChan_polscint_fs+1){
	    for(int k = 0; k<NChan_polscint_fs; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      gain_polscint_fs[k] = stemp.Atof()*qe;
	    }
	  }else{
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_polscint_fs; k++){
	      gain_polscint_fs[k] = stemp.Atof()*qe;
	    }
	  }
	  nparam_prpolscint_fs_read++;
	}
	
	if(skey=="nstages_polscint_fs"){
	  nstages_polscint_fs.resize(NChan_polscint_fs);
	  if(ntokens==NChan_polscint_fs+1){
	    for(int k = 0; k<NChan_polscint_fs; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      nstages_polscint_fs[k] = stemp.Atoi();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_polscint_fs << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_polscint_fs; k++){
	      nstages_polscint_fs[k] = stemp.Atoi();
	    }
	  }
	  nparam_prpolscint_fs_read++;
	}

	if(skey=="ped_polscint_fs"){
	  // TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  // ped_polscint_fs = stemp.Atof();
	  ped_polscint_fs.resize(NChan_polscint_fs);
	  if(ntokens==NChan_polscint_fs+1){
	    for(int k = 0; k<NChan_polscint_fs; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      ped_polscint_fs[k] = stemp.Atof();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_polscint_fs << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_polscint_fs; k++){
	      ped_polscint_fs[k] = stemp.Atof();
	    }
	  }
	  nparam_prpolscint_fs_read++;
	}
	
	if(skey=="pedsigma_polscint_fs"){
	  // TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  // pedsigma_polscint_fs = stemp.Atof();
	  pedsigma_polscint_fs.resize(NChan_polscint_fs);
	  if(ntokens==NChan_polscint_fs+1){
	    for(int k = 0; k<NChan_polscint_fs; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      pedsigma_polscint_fs[k] = stemp.Atof();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_polscint_fs << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_polscint_fs; k++){
	      pedsigma_polscint_fs[k] = stemp.Atof();
	    }
	  }
	  nparam_prpolscint_fs_read++;
	}
	
	if(skey=="trigoffset_polscint_fs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_polscint_fs = stemp.Atof();
	  nparam_prpolscint_fs_read++;
	}
	
	if(skey=="threshold_polscint_fs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  threshold_polscint_fs = stemp.Atof();
	  nparam_prpolscint_fs_read++;
	}
	
	if(skey=="ADCconv_polscint_fs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_polscint_fs = stemp.Atof();
	  nparam_prpolscint_fs_read++;
	}	

	if(skey=="ADCbits_polscint_fs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_polscint_fs = stemp.Atoi();
	  nparam_prpolscint_fs_read++;
	}	
	
	if(skey=="TDCconv_polscint_fs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCconv_polscint_fs = stemp.Atof();
	  nparam_prpolscint_fs_read++;
	}	
	
	if(skey=="TDCbits_polscint_fs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCbits_polscint_fs = stemp.Atoi();
	  nparam_prpolscint_fs_read++;
	}	
	
	if(skey=="sigmapulse_polscint_fs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_polscint_fs = stemp.Atof();
	  nparam_prpolscint_fs_read++;
	}

	//GEn-RP activeana
	if(skey=="NChan_activeana"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_activeana = stemp.Atoi();
	  nparam_activeana_read++;
	}
	
	if(skey=="gatewidth_activeana"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_activeana = stemp.Atof();
	  nparam_activeana_read++;
	}
	
	if(skey=="gain_activeana"){
	  gain_activeana.resize(NChan_activeana);
	  if(ntokens==NChan_activeana+1){
	    for(int k = 0; k<NChan_activeana; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      gain_activeana[k] = stemp.Atof();
	    }
	  }else{
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_activeana; k++){
	      gain_activeana[k] = stemp.Atof();
	    }
	  }
	  nparam_activeana_read++;
	}
	
	if(skey=="nstages_activeana"){
	  nstages_activeana.resize(NChan_activeana);
	  if(ntokens==NChan_activeana+1){
	    for(int k = 0; k<NChan_activeana; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      nstages_activeana[k] = stemp.Atoi();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_activeana << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_activeana; k++){
	      nstages_activeana[k] = stemp.Atoi();
	    }
	  }
	  nparam_activeana_read++;
	}
	
	if(skey=="ped_activeana"){
	  // TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  // ped_activeana = stemp.Atof();
	  ped_activeana.resize(NChan_activeana);
	  if(ntokens==NChan_activeana+1){
	    for(int k = 0; k<NChan_activeana; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      ped_activeana[k] = stemp.Atof();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_activeana << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_activeana; k++){
	      ped_activeana[k] = stemp.Atof();
	    }
	  }
	  nparam_activeana_read++;
	}
	
	if(skey=="pedsigma_activeana"){
	  // TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  // pedsigma_activeana = stemp.Atof();
	  pedsigma_activeana.resize(NChan_activeana);
	  if(ntokens==NChan_activeana+1){
	    for(int k = 0; k<NChan_activeana; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k+1] )->GetString();
	      pedsigma_activeana[k] = stemp.Atof();
	    }
	  }else{
	    cout << ntokens-1 << " entries for " << skey << " dont match number of channels = "  
		 << NChan_activeana << endl <<  " applying first value on all planes " << endl;
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    for(int k = 0; k<NChan_activeana; k++){
	      pedsigma_activeana[k] = stemp.Atof();
	    }
	  }
	  nparam_activeana_read++;
	}
	
	if(skey=="trigoffset_activeana"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_activeana = stemp.Atof();
	  nparam_activeana_read++;
	}
	
	if(skey=="threshold_activeana"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  threshold_activeana = stemp.Atof();
	  nparam_activeana_read++;
	}
	
	if(skey=="ADCconv_activeana"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_activeana = stemp.Atof();
	  nparam_activeana_read++;
	}	

	if(skey=="ADCbits_activeana"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_activeana = stemp.Atoi();
	  nparam_activeana_read++;
	}	
	
	if(skey=="TDCconv_activeana"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCconv_activeana = stemp.Atof();
	  nparam_activeana_read++;
	}	
	
	if(skey=="TDCbits_activeana"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCbits_activeana = stemp.Atoi();
	  nparam_activeana_read++;
	}	
	
	if(skey=="sigmapulse_activeana"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_activeana = stemp.Atof();
	  nparam_activeana_read++;
	}
	


	//GEMs
	if(skey=="NPlanes_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NPlanes_bbgem = stemp.Atoi();
	  
	  layer_bbgem = new Int_t[NPlanes_bbgem];
	  nstrips_bbgem = new Int_t[NPlanes_bbgem];
	  offset_bbgem = new Double_t[NPlanes_bbgem];
	  RO_angle_bbgem = new Double_t[NPlanes_bbgem];
	  triggeroffset_bbgem = new Double_t[NPlanes_bbgem/2];
	  gain_bbgem = new Double_t[NPlanes_bbgem/2];
	  nparam_bbgem_read++;
	}
	
	if(skey=="gatewidth_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_bbgem = stemp.Atof();
	  nparam_bbgem_read++;
	}
		
	if(skey=="ZsupThr_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ZsupThr_bbgem = stemp.Atof();
	  nparam_bbgem_read++;
	}

	if(skey=="nlayers_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  Nlayers_bbgem = stemp.Atof();
	  nparam_bbgem_read++;
	}
	
	if(skey=="bbgem_layer_z"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  if(ntokens==Nlayers_bbgem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      bbgem_layer_z.push_back(stemp.Atof());
	    }
	  }else{
	    cout << "number of entries for bbgem_layer_z = " << ntokens-1 
		 << " don't match nlayers = " << Nlayers_bbgem << endl;
	    cout << "fix your db " << endl;
	  }
	  nparam_bbgem_read++;
	}
	
	if(skey=="layer_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_bbgem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      layer_bbgem[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for layer_bbgem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_bbgem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_bbgem_read++;
	}
	
	if(skey=="nstrips_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_bbgem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      nstrips_bbgem[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for nstrips_bbgem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_bbgem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_bbgem_read++;
	}
	
	if(skey=="offset_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_bbgem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      offset_bbgem[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for offset_bbgem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_bbgem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_bbgem_read++;
	}
	
	if(skey=="RO_angle_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_bbgem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      RO_angle_bbgem[k-1] = stemp.Atof()*TMath::DegToRad();
	    }
	  }else{
	    cout << "number of entries for RO_angle_bbgem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_bbgem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_bbgem_read++;
	}
	
	if(skey=="triggeroffset_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_bbgem/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      triggeroffset_bbgem[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for triggeroffset_bbgem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_bbgem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_bbgem_read++;
	}
	
	if(skey=="gain_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_bbgem/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      gain_bbgem[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for gain_bbgem = " << ntokens-1 
		 << " don't match Nplanes/2 = " << NPlanes_bbgem/2 << endl;
	    if(ntokens>=2){
	      cout << "applying first value on all planes " << endl;
	      TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	      for(int k = 0; k<NPlanes_bbgem/2; k++){
		gain_bbgem[k] = stemp.Atof();
	      }
	    }else{
	      cout << "fix your db " << endl;
	      exit(-1);
	    }
	  }
	  nparam_bbgem_read++;
	}
	
	if(skey=="commonmode_array_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  commonmode_array_bbgem = new Double_t[ntokens-1];
	  for(int k = 1; k<ntokens; k++){
	    nAPV_bbgem++;
	    TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	    commonmode_array_bbgem[k-1] = stemp.Atof();
	  }
	  nparam_bbgem_read++;
	}
	
	//SBSGEM
	if(skey=="NPlanes_sbsgem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NPlanes_sbsgem = stemp.Atoi();
	  
	  layer_sbsgem = new Int_t[NPlanes_sbsgem];
	  nstrips_sbsgem = new Int_t[NPlanes_sbsgem];
	  offset_sbsgem = new Double_t[NPlanes_sbsgem];
	  RO_angle_sbsgem = new Double_t[NPlanes_sbsgem];
	  triggeroffset_sbsgem = new Double_t[NPlanes_sbsgem/2];
	  gain_sbsgem = new Double_t[NPlanes_sbsgem/2];
	  nparam_sbsgem_read++;
	}
	
	if(skey=="gatewidth_sbsgem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_sbsgem = stemp.Atof();
	  nparam_sbsgem_read++;
	}
		
	if(skey=="ZsupThr_sbsgem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ZsupThr_sbsgem = stemp.Atof();
	  nparam_sbsgem_read++;
	}

	if(skey=="nlayers_sbsgem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  Nlayers_sbsgem = stemp.Atof();
	  nparam_sbsgem_read++;
	}
	
	if(skey=="sbsgem_layer_z"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  if(ntokens==Nlayers_sbsgem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      sbsgem_layer_z.push_back(stemp.Atof());
	    }
	  }else{
	    cout << "number of entries for sbsgem_layer_z = " << ntokens-1 
		 << " don't match nlayers = " << Nlayers_sbsgem << endl;
	    cout << "fix your db " << endl;
	  }
	  nparam_sbsgem_read++;
	}
	
	if(skey=="layer_sbsgem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_sbsgem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      layer_sbsgem[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for layer_sbsgem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_sbsgem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_sbsgem_read++;
	}
	
	if(skey=="nstrips_sbsgem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_sbsgem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      nstrips_sbsgem[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for nstrips_sbsgem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_sbsgem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_sbsgem_read++;
	}
	
	if(skey=="offset_sbsgem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_sbsgem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      offset_sbsgem[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for offset_sbsgem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_sbsgem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_sbsgem_read++;
	}
	
	if(skey=="RO_angle_sbsgem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_sbsgem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      RO_angle_sbsgem[k-1] = stemp.Atof()*TMath::DegToRad();
	    }
	  }else{
	    cout << "number of entries for RO_angle_sbsgem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_sbsgem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_sbsgem_read++;
	}
	
	if(skey=="triggeroffset_sbsgem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_sbsgem/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      triggeroffset_sbsgem[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for triggeroffset_sbsgem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_sbsgem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_sbsgem_read++;
	}
	
	if(skey=="gain_sbsgem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_sbsgem/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      gain_sbsgem[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for gain_sbsgem = " << ntokens-1 
		 << " don't match Nplanes/2 = " << NPlanes_sbsgem/2 << endl;
	    if(ntokens>=2){
	      cout << "applying first value on all planes " << endl;
	      TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	      for(int k = 0; k<NPlanes_sbsgem/2; k++){
		gain_sbsgem[k] = stemp.Atof();
	      }
	    }else{
	      cout << "fix your db " << endl;
	      exit(-1);
	    }
	  }
	  nparam_sbsgem_read++;
	}
	
	if(skey=="commonmode_array_sbsgem"){
	  cout << "reading " << skey.Data() << endl;
	  commonmode_array_sbsgem = new Double_t[ntokens-1];
	  for(int k = 1; k<ntokens; k++){
	    nAPV_sbsgem++;
	    TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	    commonmode_array_sbsgem[k-1] = stemp.Atof();
	  }
	  nparam_sbsgem_read++;
	}
	
	// GEMs GEn-RP front
	if(skey=="NPlanes_cepol_front"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NPlanes_cepol_front = stemp.Atoi();
	  
	  layer_cepol_front = new Int_t[NPlanes_cepol_front];
	  nstrips_cepol_front = new Int_t[NPlanes_cepol_front];
	  offset_cepol_front = new Double_t[NPlanes_cepol_front];
	  RO_angle_cepol_front = new Double_t[NPlanes_cepol_front];
	  triggeroffset_cepol_front = new Double_t[NPlanes_cepol_front/2];
	  gain_cepol_front = new Double_t[NPlanes_cepol_front/2];
	  nparam_cepol_front_read++;
	}
	
	if(skey=="gatewidth_cepol_front"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_cepol_front = stemp.Atof();
	  nparam_cepol_front_read++;
	}
		
	if(skey=="ZsupThr_cepol_front"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ZsupThr_cepol_front = stemp.Atof();
	  nparam_cepol_front_read++;
	}

	if(skey=="nlayers_cepol_front"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  Nlayers_cepol_front = stemp.Atof();
	  nparam_cepol_front_read++;
	}
	
	if(skey=="cepol_front_layer_z"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  if(ntokens==Nlayers_cepol_front+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      cepol_front_layer_z.push_back(stemp.Atof());
	    }
	  }else{
	    cout << "number of entries for cepol_front_layer_z = " << ntokens-1 
		 << " don't match nlayers = " << Nlayers_cepol_front << endl;
	    cout << "fix your db " << endl;
	  }
	  nparam_cepol_front_read++;
	}
	
	if(skey=="layer_cepol_front"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_cepol_front+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      layer_cepol_front[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for layer_cepol_front = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_cepol_front << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_cepol_front_read++;
	}
	
	if(skey=="nstrips_cepol_front"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_cepol_front+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      nstrips_cepol_front[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for nstrips_cepol_front = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_cepol_front << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_cepol_front_read++;
	}
	
	if(skey=="offset_cepol_front"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_cepol_front+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      offset_cepol_front[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for offset_cepol_front = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_cepol_front << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_cepol_front_read++;
	}
	
	if(skey=="RO_angle_cepol_front"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_cepol_front+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      RO_angle_cepol_front[k-1] = stemp.Atof()*TMath::DegToRad();
	    }
	  }else{
	    cout << "number of entries for RO_angle_cepol_front = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_cepol_front << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_cepol_front_read++;
	}
	
	if(skey=="triggeroffset_cepol_front"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_cepol_front/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      triggeroffset_cepol_front[k-1] = stemp.Atof();
	      // if this is affected at k-1, program crashes... it should not and that's confusing...
	    }
	  }else{
	    cout << "number of entries for triggeroffset_cepol_front = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_cepol_front << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_cepol_front_read++;
	}
	
	if(skey=="gain_cepol_front"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_cepol_front/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      gain_cepol_front[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for gain_cepol_front = " << ntokens-1 
		 << " don't match Nplanes/2 = " << NPlanes_cepol_front/2 << endl;
	    if(ntokens>=2){
	      cout << "applying first value on all planes " << endl;
	      TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	      for(int k = 0; k<NPlanes_cepol_front/2; k++){
		gain_cepol_front[k] = stemp.Atof();
	      }
	    }else{
	      cout << "fix your db " << endl;
	      exit(-1);
	    }
	  }
	  nparam_cepol_front_read++;
	}
		
	if(skey=="commonmode_array_cepol_front"){
	  cout << "reading " << skey.Data() << endl;
	  commonmode_array_cepol_front = new Double_t[ntokens-1];
	  for(int k = 1; k<ntokens; k++){
	    nAPV_cepol_front++;
	    TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	    commonmode_array_cepol_front[k-1] = stemp.Atof();
	  }
	  nparam_cepol_front_read++;
	}

	// GEMs GEn-RP rear
	if(skey=="NPlanes_cepol_rear"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NPlanes_cepol_rear = stemp.Atoi();
	  
	  layer_cepol_rear = new Int_t[NPlanes_cepol_rear];
	  nstrips_cepol_rear = new Int_t[NPlanes_cepol_rear];
	  offset_cepol_rear = new Double_t[NPlanes_cepol_rear];
	  RO_angle_cepol_rear = new Double_t[NPlanes_cepol_rear];
	  triggeroffset_cepol_rear = new Double_t[NPlanes_cepol_rear/2];
	  gain_cepol_rear = new Double_t[NPlanes_cepol_rear/2];
	  nparam_cepol_rear_read++;
	}
	
	if(skey=="gatewidth_cepol_rear"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_cepol_rear = stemp.Atof();
	  nparam_cepol_rear_read++;
	}
		
	if(skey=="ZsupThr_cepol_rear"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ZsupThr_cepol_rear = stemp.Atof();
	  nparam_cepol_rear_read++;
	}

	if(skey=="nlayers_cepol_rear"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  Nlayers_cepol_rear = stemp.Atof();
	  nparam_cepol_rear_read++;
	}
	
	if(skey=="cepol_rear_layer_z"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  if(ntokens==Nlayers_cepol_rear+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      cepol_rear_layer_z.push_back(stemp.Atof());
	    }
	  }else{
	    cout << "number of entries for cepol_rear_layer_z = " << ntokens-1 
		 << " don't match nlayers = " << Nlayers_cepol_rear << endl;
	    cout << "fix your db " << endl;
	  }
	  nparam_cepol_rear_read++;
	}
	
	if(skey=="layer_cepol_rear"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_cepol_rear+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      layer_cepol_rear[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for layer_cepol_rear = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_cepol_rear << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_cepol_rear_read++;
	}
	
	if(skey=="nstrips_cepol_rear"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_cepol_rear+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      nstrips_cepol_rear[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for nstrips_cepol_rear = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_cepol_rear << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_cepol_rear_read++;
	}
	
	if(skey=="offset_cepol_rear"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_cepol_rear+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      offset_cepol_rear[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for offset_cepol_rear = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_cepol_rear << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_cepol_rear_read++;
	}
	
	if(skey=="RO_angle_cepol_rear"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_cepol_rear+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      RO_angle_cepol_rear[k-1] = stemp.Atof()*TMath::DegToRad();
	    }
	  }else{
	    cout << "number of entries for RO_angle_cepol_rear = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_cepol_rear << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_cepol_rear_read++;
	}
	
	if(skey=="triggeroffset_cepol_rear"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_cepol_rear/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      triggeroffset_cepol_rear[k-1] = stemp.Atof();
	      // if this is affected at k-1, program crashes... it should not and that's confusing...
	    }
	  }else{
	    cout << "number of entries for triggeroffset_cepol_rear = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_cepol_rear << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_cepol_rear_read++;
	}
	
	if(skey=="gain_cepol_rear"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_cepol_rear/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      gain_cepol_rear[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for gain_cepol_rear = " << ntokens-1 
		 << " don't match Nplanes/2 = " << NPlanes_cepol_rear/2 << endl;
	    if(ntokens>=2){
	      cout << "applying first value on all planes " << endl;
	      TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	      for(int k = 0; k<NPlanes_cepol_rear/2; k++){
		gain_cepol_rear[k] = stemp.Atof();
	      }
	    }else{
	      cout << "fix your db " << endl;
	      exit(-1);
	    }
	  }
	  nparam_cepol_rear_read++;
	}
		
	if(skey=="commonmode_array_cepol_rear"){
	  cout << "reading " << skey.Data() << endl;
	  commonmode_array_cepol_rear = new Double_t[ntokens-1];
	  for(int k = 1; k<ntokens; k++){
	    nAPV_cepol_rear++;
	    TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	    commonmode_array_cepol_rear[k-1] = stemp.Atof();
	  }
	  nparam_cepol_rear_read++;
	}

	// GEMs prpolbs_gem
	if(skey=="NPlanes_prpolbs_gem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NPlanes_prpolbs_gem = stemp.Atoi();
	  
	  layer_prpolbs_gem = new Int_t[NPlanes_prpolbs_gem];
	  nstrips_prpolbs_gem = new Int_t[NPlanes_prpolbs_gem];
	  offset_prpolbs_gem = new Double_t[NPlanes_prpolbs_gem];
	  RO_angle_prpolbs_gem = new Double_t[NPlanes_prpolbs_gem];
	  triggeroffset_prpolbs_gem = new Double_t[NPlanes_prpolbs_gem/2];
	  gain_prpolbs_gem = new Double_t[NPlanes_prpolbs_gem/2];
	  nparam_prpolbs_gem_read++;
	}
	
	if(skey=="gatewidth_prpolbs_gem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_prpolbs_gem = stemp.Atof();
	  nparam_prpolbs_gem_read++;
	}
		
	if(skey=="ZsupThr_prpolbs_gem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ZsupThr_prpolbs_gem = stemp.Atof();
	  nparam_prpolbs_gem_read++;
	}

	if(skey=="nlayers_prpolbs_gem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  Nlayers_prpolbs_gem = stemp.Atof();
	  nparam_prpolbs_gem_read++;
	}
	
	if(skey=="prpolbs_gem_layer_z"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  if(ntokens==Nlayers_prpolbs_gem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      prpolbs_gem_layer_z.push_back(stemp.Atof());
	    }
	  }else{
	    cout << "number of entries for prpolbs_gem_layer_z = " << ntokens-1 
		 << " don't match nlayers = " << Nlayers_prpolbs_gem << endl;
	    cout << "fix your db " << endl;
	  }
	  nparam_prpolbs_gem_read++;
	}
	
	if(skey=="layer_prpolbs_gem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_prpolbs_gem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      layer_prpolbs_gem[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for layer_prpolbs_gem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_prpolbs_gem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_prpolbs_gem_read++;
	}
	
	if(skey=="nstrips_prpolbs_gem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_prpolbs_gem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      nstrips_prpolbs_gem[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for nstrips_prpolbs_gem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_prpolbs_gem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_prpolbs_gem_read++;
	}
	
	if(skey=="offset_prpolbs_gem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_prpolbs_gem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      offset_prpolbs_gem[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for offset_prpolbs_gem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_prpolbs_gem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_prpolbs_gem_read++;
	}
	
	if(skey=="RO_angle_prpolbs_gem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_prpolbs_gem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      RO_angle_prpolbs_gem[k-1] = stemp.Atof()*TMath::DegToRad();
	    }
	  }else{
	    cout << "number of entries for RO_angle_prpolbs_gem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_prpolbs_gem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_prpolbs_gem_read++;
	}
	
	if(skey=="triggeroffset_prpolbs_gem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_prpolbs_gem/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      triggeroffset_prpolbs_gem[k-1] = stemp.Atof();
	      // if this is affected at k-1, program crashes... it should not and that's confusing...
	    }
	  }else{
	    cout << "number of entries for triggeroffset_prpolbs_gem = " << ntokens-1 
		 << " don't match Nplanes/2 = " << NPlanes_prpolbs_gem/2 << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_prpolbs_gem_read++;
	}
	
	if(skey=="gain_prpolbs_gem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_prpolbs_gem/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      gain_prpolbs_gem[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for gain_prpolbs_gem = " << ntokens-1 
		 << " don't match Nplanes/2 = " << NPlanes_prpolbs_gem/2 << endl;
	    if(ntokens>=2){
	      cout << "applying first value on all planes " << endl;
	      TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	      for(int k = 0; k<NPlanes_prpolbs_gem/2; k++){
		gain_prpolbs_gem[k] = stemp.Atof();
	      }
	    }else{
	      cout << "fix your db " << endl;
	      exit(-1);
	    }
	  }
	  nparam_prpolbs_gem_read++;
	}
	
	if(skey=="commonmode_array_prpolbs_gem"){
	  cout << "reading " << skey.Data() << endl;
	  commonmode_array_prpolbs_gem = new Double_t[ntokens-1];
	  for(int k = 1; k<ntokens; k++){
	    nAPV_prpolbs_gem++;
	    TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	    commonmode_array_prpolbs_gem[k-1] = stemp.Atof();
	  }
	  nparam_prpolbs_gem_read++;
	}
	// GEMs prpolfs_gem
	if(skey=="NPlanes_prpolfs_gem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NPlanes_prpolfs_gem = stemp.Atoi();
	  
	  layer_prpolfs_gem = new Int_t[NPlanes_prpolfs_gem];
	  nstrips_prpolfs_gem = new Int_t[NPlanes_prpolfs_gem];
	  offset_prpolfs_gem = new Double_t[NPlanes_prpolfs_gem];
	  RO_angle_prpolfs_gem = new Double_t[NPlanes_prpolfs_gem];
	  triggeroffset_prpolfs_gem = new Double_t[NPlanes_prpolfs_gem/2];
	  gain_prpolfs_gem = new Double_t[NPlanes_prpolfs_gem/2];
	  nparam_prpolfs_gem_read++;
	}
	
	if(skey=="gatewidth_prpolfs_gem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_prpolfs_gem = stemp.Atof();
	  nparam_prpolfs_gem_read++;
	}
		
	if(skey=="ZsupThr_prpolfs_gem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ZsupThr_prpolfs_gem = stemp.Atof();
	  nparam_prpolfs_gem_read++;
	}

	if(skey=="nlayers_prpolfs_gem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  Nlayers_prpolfs_gem = stemp.Atof();
	  nparam_prpolfs_gem_read++;
	}
	
	if(skey=="prpolfs_gem_layer_z"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  if(ntokens==Nlayers_prpolfs_gem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      prpolfs_gem_layer_z.push_back(stemp.Atof());
	    }
	  }else{
	    cout << "number of entries for prpolfs_gem_layer_z = " << ntokens-1 
		 << " don't match nlayers = " << Nlayers_prpolfs_gem << endl;
	    cout << "fix your db " << endl;
	  }
	  nparam_prpolfs_gem_read++;
	}
	
	if(skey=="layer_prpolfs_gem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_prpolfs_gem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      layer_prpolfs_gem[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for layer_prpolfs_gem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_prpolfs_gem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_prpolfs_gem_read++;
	}
	
	if(skey=="nstrips_prpolfs_gem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_prpolfs_gem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      nstrips_prpolfs_gem[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for nstrips_prpolfs_gem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_prpolfs_gem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_prpolfs_gem_read++;
	}
	
	if(skey=="offset_prpolfs_gem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_prpolfs_gem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      offset_prpolfs_gem[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for offset_prpolfs_gem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_prpolfs_gem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_prpolfs_gem_read++;
	}
	
	if(skey=="RO_angle_prpolfs_gem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_prpolfs_gem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      RO_angle_prpolfs_gem[k-1] = stemp.Atof()*TMath::DegToRad();
	    }
	  }else{
	    cout << "number of entries for RO_angle_prpolfs_gem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_prpolfs_gem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_prpolfs_gem_read++;
	}
	
	if(skey=="triggeroffset_prpolfs_gem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_prpolfs_gem/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      triggeroffset_prpolfs_gem[k-1] = stemp.Atof();
	      // if this is affected at k-1, program crashes... it should not and that's confusing...
	    }
	  }else{
	    cout << "number of entries for triggeroffset_prpolfs_gem = " << ntokens-1 
		 << " don't match Nplanes/2 = " << NPlanes_prpolfs_gem/2 << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_prpolfs_gem_read++;
	}
	
	if(skey=="gain_prpolfs_gem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_prpolfs_gem/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      gain_prpolfs_gem[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for gain_prpolfs_gem = " << ntokens-1 
		 << " don't match Nplanes/2 = " << NPlanes_prpolfs_gem/2 << endl;
	    if(ntokens>=2){
	      cout << "applying first value on all planes " << endl;
	      TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	      for(int k = 0; k<NPlanes_prpolfs_gem/2; k++){
		gain_prpolfs_gem[k] = stemp.Atof();
	      }
	    }else{
	      cout << "fix your db " << endl;
	      exit(-1);
	    }
	  }
	  nparam_prpolfs_gem_read++;
	}
	
	if(skey=="commonmode_array_prpolfs_gem"){
	  cout << "reading " << skey.Data() << endl;
	  commonmode_array_prpolfs_gem = new Double_t[ntokens-1];
	  for(int k = 1; k<ntokens; k++){
	    nAPV_prpolfs_gem++;
	    TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	    commonmode_array_prpolfs_gem[k-1] = stemp.Atof();
	  }
	  nparam_prpolfs_gem_read++;
	}

	//FT
	if(skey=="NPlanes_ft"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NPlanes_ft = stemp.Atoi();
	  
	  layer_ft = new Int_t[NPlanes_ft];
	  nstrips_ft = new Int_t[NPlanes_ft];
	  offset_ft = new Double_t[NPlanes_ft];
	  RO_angle_ft = new Double_t[NPlanes_ft];
	  triggeroffset_ft = new Double_t[NPlanes_ft/2];
	  gain_ft = new Double_t[NPlanes_ft/2];
	  nparam_ft_read++;
	}
	
	if(skey=="gatewidth_ft"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_ft = stemp.Atof();
	  nparam_ft_read++;
	}
	
	if(skey=="ZsupThr_ft"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ZsupThr_ft = stemp.Atof();
	  nparam_ft_read++;
	}

	if(skey=="nlayers_ft"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  Nlayers_ft = stemp.Atof();
	  nparam_ft_read++;
	}
	
	if(skey=="ft_layer_z"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  if(ntokens==Nlayers_ft+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      ft_layer_z.push_back(stemp.Atof());
	    }
	  }else{
	    cout << "number of entries for ft_layer_z = " << ntokens-1 
		 << " don't match nlayers = " << Nlayers_ft << endl;
	    cout << "fix your db " << endl;
	  }
	  nparam_ft_read++;
	}
		
	if(skey=="layer_ft"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_ft+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      layer_ft[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for layer_ft = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_ft << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_ft_read++;
	}
	
	if(skey=="nstrips_ft"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_ft+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      nstrips_ft[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for nstrips_ft = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_ft << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_ft_read++;
	}
	
	if(skey=="offset_ft"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_ft+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      offset_ft[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for offset_ft = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_ft << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_ft_read++;
	}
	
	if(skey=="RO_angle_ft"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_ft+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      RO_angle_ft[k-1] = stemp.Atof()*TMath::DegToRad();
	    }
	  }else{
	    cout << "number of entries for RO_angle_ft = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_ft << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_ft_read++;
	}
	
	if(skey=="triggeroffset_ft"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_ft/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      triggeroffset_ft[k-1] = stemp.Atof();
	      // if this is affected at k-1, program crashes... it should not and that's confusing...
	    }
	  }else{
	    cout << "number of entries for triggeroffset_ft = " << ntokens-1 
		 << " don't match Nplanes/2 = " << NPlanes_ft/2 << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_ft_read++;
	}
	
	if(skey=="gain_ft"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_ft/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      gain_ft[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for gain_ft = " << ntokens-1 
		 << " don't match Nplanes/2 = " << NPlanes_ft/2 << endl;
	    if(ntokens>=2){
	      cout << "applying first value on all planes " << endl;
	      TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	      for(int k = 0; k<NPlanes_ft/2; k++){
		gain_ft[k] = stemp.Atof();
	      }
	    }else{
	      cout << "fix your db " << endl;
	      exit(-1);
	    }
	  }
	  nparam_ft_read++;
	}

	if(skey=="commonmode_array_ft"){
	  cout << "reading " << skey.Data() << endl;
	  commonmode_array_ft = new Double_t[ntokens-1];
	  for(int k = 1; k<ntokens; k++){
	    nAPV_ft++;
	    TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	    commonmode_array_ft[k-1] = stemp.Atof();
	  }
	  nparam_ft_read++;
	}

  if(skey=="do_pedcm_ft"){
    cout << "reading " << skey.Data() << endl;
    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
    do_pedcm_ft = stemp.Atoi();
    // nparam_ft_read++; // Intentionally not incrementing 'nparam_ft_read' to ensure backward compatibility with old database files.
  }

  if(skey=="pedfile_ft"){
    cout << "reading " << skey.Data() << endl;
    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
    pedfile_ft = stemp.Data();
    // nparam_ft_read++; // Intentionally not incrementing 'nparam_ft_read' to ensure backward compatibility with old database files.
  }

  if(skey=="cmfile_ft"){
    cout << "reading " << skey.Data() << endl;
    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
    cmfile_ft = stemp.Data();
    // nparam_ft_read++; // Intentionally not incrementing 'nparam_ft_read' to ensure backward compatibility with old database files.
  }
  
  if(skey=="cmfitfile_ft"){
    cout << "reading " << skey.Data() << endl;
    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
    cmfitfile_ft = stemp.Data();
  }

  if(skey=="do_online_cm_ft"){
    cout << "reading " << skey.Data() << endl;
    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
    do_online_cm_ft = stemp.Atoi();
  }

  if(skey=="do_online_zs_ft"){
    cout << "reading " << skey.Data() << endl;
    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
    do_online_zs_ft = stemp.Atoi();
  }

  if (skey=="online_zs_thr_nsigma_ft"){
    cout << "reading " << skey.Data() << endl;
    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
    online_zs_thr_nsigma_ft = stemp.Atoi();
  }

  if (skey=="make_online_cmplots_ft"){
    cout << "reading " << skey.Data() << endl;
    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
    make_online_cmplots_ft = stemp.Atoi();
  }

	
	//FPP1
	if(skey=="NPlanes_fpp1"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NPlanes_fpp1 = stemp.Atoi();
	  
	  layer_fpp1 = new Int_t[NPlanes_fpp1];
	  nstrips_fpp1 = new Int_t[NPlanes_fpp1];
	  offset_fpp1 = new Double_t[NPlanes_fpp1];
	  RO_angle_fpp1 = new Double_t[NPlanes_fpp1];
	  triggeroffset_fpp1 = new Double_t[NPlanes_fpp1/2];
	  gain_fpp1 = new Double_t[NPlanes_fpp1/2];
	  nparam_fpp1_read++;
	}
	
	if(skey=="gatewidth_fpp1"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_fpp1 = stemp.Atof();
	  nparam_fpp1_read++;
	}
	
	if(skey=="ZsupThr_fpp1"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ZsupThr_fpp1 = stemp.Atof();
	  nparam_fpp1_read++;
	}
	
	if(skey=="nlayers_fpp1"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  Nlayers_fpp1 = stemp.Atof();
	  nparam_fpp1_read++;
	}
	
	if(skey=="fpp1_layer_z"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  if(ntokens==Nlayers_fpp1+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      fpp1_layer_z.push_back(stemp.Atof());
	    }
	  }else{
	    cout << "number of entries for fpp1_layer_z = " << ntokens-1 
		 << " don't match nlayers = " << Nlayers_fpp1 << endl;
	    cout << "fix your db " << endl;
	  }
	  nparam_fpp1_read++;
	}
	
	if(skey=="layer_fpp1"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_fpp1+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      layer_fpp1[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for layer_fpp1 = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_fpp1 << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_fpp1_read++;
	}
	
	if(skey=="nstrips_fpp1"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_fpp1+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      nstrips_fpp1[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for nstrips_fpp1 = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_fpp1 << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_fpp1_read++;
	}
	
	if(skey=="offset_fpp1"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_fpp1+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      offset_fpp1[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for offset_fpp1 = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_fpp1 << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_fpp1_read++;
	}
	
	if(skey=="RO_angle_fpp1"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_fpp1+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      RO_angle_fpp1[k-1] = stemp.Atof()*TMath::DegToRad();
	    }
	  }else{
	    cout << "number of entries for RO_angle_fpp1 = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_fpp1 << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_fpp1_read++;
	}
	
	if(skey=="triggeroffset_fpp1"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_fpp1/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      triggeroffset_fpp1[k-1] = stemp.Atof();
	      // if this is affected at k-1, program crashes... it should not and that's confusing...
	    }
	  }else{
	    cout << "number of entries for triggeroffset_fpp1 = " << ntokens-1 
		 << " don't match Nplanes/2 = " << NPlanes_fpp1/2 << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_fpp1_read++;
	}
	
	if(skey=="gain_fpp1"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_fpp1/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      gain_fpp1[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for gain_fpp1 = " << ntokens-1 
		 << " don't match Nplanes/2 = " << NPlanes_fpp1/2 << endl;
	    if(ntokens>=2){
	      cout << "applying first value on all planes " << endl;
	      TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	      for(int k = 0; k<NPlanes_fpp1/2; k++){
		gain_fpp1[k] = stemp.Atof();
	      }
	    }else{
	      cout << "fix your db " << endl;
	      exit(-1);
	    }
	  }
	  nparam_fpp1_read++;
	}

	if(skey=="commonmode_array_fpp1"){
	  cout << "reading " << skey.Data() << endl;
	  commonmode_array_fpp1 = new Double_t[ntokens-1];
	  for(int k = 1; k<ntokens; k++){
	    nAPV_fpp1++;
	    TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	    commonmode_array_fpp1[k-1] = stemp.Atof();
	  }
	  nparam_fpp1_read++;
	}
	
	//FPP2
	if(skey=="NPlanes_fpp2"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NPlanes_fpp2 = stemp.Atoi();
	  
	  layer_fpp2 = new Int_t[NPlanes_fpp2];
	  nstrips_fpp2 = new Int_t[NPlanes_fpp2];
	  offset_fpp2 = new Double_t[NPlanes_fpp2];
	  RO_angle_fpp2 = new Double_t[NPlanes_fpp2];
	  triggeroffset_fpp2 = new Double_t[NPlanes_fpp2/2];
	  gain_fpp2 = new Double_t[NPlanes_fpp2/2];
	  nparam_fpp2_read++;
	}
	
	if(skey=="gatewidth_fpp2"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_fpp2 = stemp.Atof();
	  nparam_fpp2_read++;
	}
	
	if(skey=="ZsupThr_fpp2"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ZsupThr_fpp2 = stemp.Atof();
	  nparam_fpp2_read++;
	}
	
	if(skey=="nlayers_fpp2"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  Nlayers_fpp2 = stemp.Atof();
	  nparam_fpp2_read++;
	}
	
	if(skey=="fpp2_layer_z"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  if(ntokens==Nlayers_fpp2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      fpp2_layer_z.push_back(stemp.Atof());
	    }
	  }else{
	    cout << "number of entries for fpp2_layer_z = " << ntokens-1 
		 << " don't match nlayers = " << Nlayers_fpp2 << endl;
	    cout << "fix your db " << endl;
	  }
	  nparam_fpp2_read++;
	}
	
	if(skey=="layer_fpp2"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_fpp2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      layer_fpp2[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for layer_fpp2 = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_fpp2 << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_fpp2_read++;
	}
	
	if(skey=="nstrips_fpp2"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_fpp2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      nstrips_fpp2[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for nstrips_fpp2 = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_fpp2 << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_fpp2_read++;
	}
	
	if(skey=="offset_fpp2"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_fpp2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      offset_fpp2[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for offset_fpp2 = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_fpp2 << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_fpp2_read++;
	}
	
	if(skey=="RO_angle_fpp2"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_fpp2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      RO_angle_fpp2[k-1] = stemp.Atof()*TMath::DegToRad();
	    }
	  }else{
	    cout << "number of entries for RO_angle_fpp2 = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_fpp2 << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_fpp2_read++;
	}
	
	if(skey=="triggeroffset_fpp2"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_fpp2/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      triggeroffset_fpp2[k-1] = stemp.Atof();
	      // if this is affected at k-1, program crashes... it should not and that's confusing...
	    }
	  }else{
	    cout << "number of entries for triggeroffset_fpp2 = " << ntokens-1 
		 << " don't match Nplanes/2 = " << NPlanes_fpp2/2 << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_fpp2_read++;
	}
	
	if(skey=="gain_fpp2"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_fpp2/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      gain_fpp2[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for gain_fpp2 = " << ntokens-1 
		 << " don't match Nplanes/2 = " << NPlanes_fpp2/2 << endl;
	    if(ntokens>=2){
	      cout << "applying first value on all planes " << endl;
	      TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	      for(int k = 0; k<NPlanes_fpp2/2; k++){
		gain_fpp2[k] = stemp.Atof();
	      }
	    }else{
	      cout << "fix your db " << endl;
	      exit(-1);
	    }
	  }
	  nparam_fpp2_read++;
	}
	
	if(skey=="commonmode_array_fpp2"){
	  cout << "reading " << skey.Data() << endl;
	  commonmode_array_fpp2 = new Double_t[ntokens-1];
	  for(int k = 1; k<ntokens; k++){
	    nAPV_fpp2++;
	    TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	    commonmode_array_fpp2[k-1] = stemp.Atof();
	  }
	  nparam_fpp2_read++;
	}
	
      }//end if( ntokens >= 2 )
      tokens->~TObjArray();// ineffective... :(
    }//end if( !currentline.BeginsWith("#"))
  }//end while
  
  //-----------------------------
  //  Declare detectors
  //-----------------------------
  cout << " declaring detectors " << endl;
  for(int k = 0; k<detectors_list.size(); k++){
    cout << "detector: " << detectors_list[k].Data() << "... " << endl;
    if(detectors_list[k] == "bbgem"){
      if(nparam_bbgem_read!=nparam_gemdet){
	cout << detectors_list[k] <<  " does not have the right number of parameters!!! " << endl << " fix database and retry! " << endl;
	exit(-1);
      }

      SBSDigGEMDet* bbgem = new SBSDigGEMDet(BBGEM_UNIQUE_DETID, NPlanes_bbgem, layer_bbgem, nstrips_bbgem, offset_bbgem, RO_angle_bbgem, 6, ZsupThr_bbgem);
      SBSDigGEMSimDig* gemdig = new SBSDigGEMSimDig(NPlanes_bbgem/2, triggeroffset_bbgem, gain_bbgem, ZsupThr_bbgem, nAPV_bbgem, commonmode_array_bbgem);
      for(int m = 0; m<Nlayers_bbgem; m++){
	bbgem->fZLayer.push_back(bbgem_layer_z[m]);
      }
      bbgem->fGateWidth = gatewidth_bbgem;
      
      GEMdetectors.push_back(bbgem);
      gemdetmap.push_back(BBGEM_UNIQUE_DETID);
      GEMsimDig.push_back(gemdig);
      cout << " set up! " << endl;
    }
    if(detectors_list[k] == "sbsgem"){
      if(nparam_sbsgem_read!=nparam_gemdet){
	cout << detectors_list[k] <<  " does not have the right number of parameters!!! " << endl << " fix database and retry! " << endl;
	exit(-1);
      }
      
      SBSDigGEMDet* sbsgem = new SBSDigGEMDet(SBSGEM_UNIQUE_DETID, NPlanes_sbsgem, layer_sbsgem, nstrips_sbsgem, offset_sbsgem, RO_angle_sbsgem, 6, ZsupThr_sbsgem);
      SBSDigGEMSimDig* gemdig = new SBSDigGEMSimDig(NPlanes_sbsgem/2, triggeroffset_sbsgem, gain_sbsgem, ZsupThr_sbsgem, nAPV_sbsgem, commonmode_array_sbsgem);
      for(int m = 0; m<Nlayers_sbsgem; m++){
	sbsgem->fZLayer.push_back(sbsgem_layer_z[m]);
      }
      sbsgem->fGateWidth = gatewidth_sbsgem;
      
      GEMdetectors.push_back(sbsgem);
      gemdetmap.push_back(SBSGEM_UNIQUE_DETID);
      GEMsimDig.push_back(gemdig);
      cout << " set up! " << endl;
    }
    if(detectors_list[k] == "cepol_front"){
      if(nparam_cepol_front_read!=nparam_gemdet){
	cout << detectors_list[k] <<  " does not have the right number of parameters!!! " << endl << " fix database and retry! " << endl;
	exit(-1);
      }

      SBSDigGEMDet* cepol_front = new SBSDigGEMDet(CEPOL_GEMFRONT_UNIQUE_DETID, NPlanes_cepol_front, layer_cepol_front, nstrips_cepol_front, offset_cepol_front, RO_angle_cepol_front, 6, ZsupThr_cepol_front);
      SBSDigGEMSimDig* gemdig = new SBSDigGEMSimDig(NPlanes_cepol_front/2, triggeroffset_cepol_front, gain_cepol_front, ZsupThr_cepol_front, nAPV_cepol_front, commonmode_array_cepol_front);
      for(int m = 0; m<Nlayers_cepol_front; m++){
	cepol_front->fZLayer.push_back(cepol_front_layer_z[m]);
      }
      cepol_front->fGateWidth = gatewidth_cepol_front;
      
      GEMdetectors.push_back(cepol_front);
      gemdetmap.push_back(CEPOL_GEMFRONT_UNIQUE_DETID);
      GEMsimDig.push_back(gemdig);
      cout << " set up! " << endl;
    }
    if(detectors_list[k] == "cepol_rear"){
      if(nparam_cepol_rear_read!=nparam_gemdet){
	cout << detectors_list[k] <<  " does not have the right number of parameters!!! " << endl << " fix database and retry! " << endl;
	exit(-1);
      }
      
      SBSDigGEMDet* cepol_rear = new SBSDigGEMDet(CEPOL_GEMREAR_UNIQUE_DETID, NPlanes_cepol_rear, layer_cepol_rear, nstrips_cepol_rear, offset_cepol_rear, RO_angle_cepol_rear, 6, ZsupThr_cepol_rear);
      SBSDigGEMSimDig* gemdig = new SBSDigGEMSimDig(NPlanes_cepol_rear/2, triggeroffset_cepol_rear, gain_cepol_rear, ZsupThr_cepol_rear, nAPV_cepol_rear, commonmode_array_cepol_rear);
      for(int m = 0; m<Nlayers_cepol_rear; m++){
	cepol_rear->fZLayer.push_back(cepol_rear_layer_z[m]);
      }
      cepol_rear->fGateWidth = gatewidth_cepol_rear;
      
      GEMdetectors.push_back(cepol_rear);
      gemdetmap.push_back(CEPOL_GEMREAR_UNIQUE_DETID);
      GEMsimDig.push_back(gemdig);
      cout << " set up! " << endl;
    }
    if(detectors_list[k] == "prpolbs_gem"){
      if(nparam_prpolbs_gem_read!=nparam_gemdet){
	cout << detectors_list[k] <<  " does not have the right number of parameters!!! " << endl << " fix database and retry! " << endl;
	exit(-1);
      }

     SBSDigGEMDet* prpolbs_gem = new SBSDigGEMDet(PRPOLBS_GEM_UNIQUE_DETID, NPlanes_prpolbs_gem, layer_prpolbs_gem, nstrips_prpolbs_gem, offset_prpolbs_gem, RO_angle_prpolbs_gem, 6, ZsupThr_prpolbs_gem);
      SBSDigGEMSimDig* gemdig = new SBSDigGEMSimDig(NPlanes_prpolbs_gem/2, triggeroffset_prpolbs_gem, gain_prpolbs_gem, ZsupThr_prpolbs_gem, nAPV_prpolbs_gem, commonmode_array_prpolbs_gem);
      for(int m = 0; m<Nlayers_prpolbs_gem; m++){
	prpolbs_gem->fZLayer.push_back(prpolbs_gem_layer_z[m]);
      }
      prpolbs_gem->fGateWidth = gatewidth_prpolbs_gem;
      
      GEMdetectors.push_back(prpolbs_gem);
      gemdetmap.push_back(PRPOLBS_GEM_UNIQUE_DETID);
      GEMsimDig.push_back(gemdig);
    }
    
    if(detectors_list[k] == "prpolfs_gem"){
      if(nparam_prpolfs_gem_read!=nparam_gemdet){
	cout << detectors_list[k] <<  " does not have the right number of parameters!!! " << endl << " fix database and retry! " << endl;
	exit(-1);
      }
      
      SBSDigGEMDet* prpolfs_gem = new SBSDigGEMDet(PRPOLFS_GEM_UNIQUE_DETID, NPlanes_prpolfs_gem, layer_prpolfs_gem, nstrips_prpolfs_gem, offset_prpolfs_gem, RO_angle_prpolfs_gem, 6, ZsupThr_prpolfs_gem);
      SBSDigGEMSimDig* gemdig = new SBSDigGEMSimDig(NPlanes_prpolfs_gem/2, triggeroffset_prpolfs_gem, gain_prpolfs_gem, ZsupThr_prpolfs_gem, nAPV_prpolfs_gem, commonmode_array_prpolfs_gem);
      for(int m = 0; m<Nlayers_prpolfs_gem; m++){
	prpolfs_gem->fZLayer.push_back(prpolfs_gem_layer_z[m]);
      }
      prpolfs_gem->fGateWidth = gatewidth_prpolfs_gem;
      
      GEMdetectors.push_back(prpolfs_gem);
      gemdetmap.push_back(PRPOLFS_GEM_UNIQUE_DETID);
      GEMsimDig.push_back(gemdig);
    }
    
    if(detectors_list[k] == "ft"){
      if(nparam_ft_read!=nparam_gemdet){
	cout << detectors_list[k] <<  " does not have the right number of parameters!!! " << endl << " fix database and retry! " << endl;
	exit(-1);
      }
      
      //SBSDigGEMDet* ft = new SBSDigGEMDet(FT_UNIQUE_DETID, NPlanes_ft, layer_ft, nstrips_ft, offset_ft, RO_angle_ft, 6, ZsupThr_ft, do_pedcm_ft, pedfile_ft, cmfile_ft);
      
      SBSDigGEMDet* ft = new SBSDigGEMDet(FT_UNIQUE_DETID, NPlanes_ft, layer_ft, nstrips_ft, offset_ft, RO_angle_ft, 6, ZsupThr_ft, do_pedcm_ft, pedfile_ft, cmfile_ft, cmfitfile_ft);
      
      SBSDigGEMSimDig* gemdig = new SBSDigGEMSimDig(NPlanes_ft/2, triggeroffset_ft, gain_ft, ZsupThr_ft, nAPV_ft, commonmode_array_ft, do_pedcm_ft, do_online_cm_ft, do_online_zs_ft, online_zs_thr_nsigma_ft, make_online_cmplots_ft);
      for(int m = 0; m<Nlayers_ft; m++){
	     ft->fZLayer.push_back(ft_layer_z[m]);
      }
      ft->fGateWidth = gatewidth_ft;
      
      GEMdetectors.push_back(ft);
      gemdetmap.push_back(FT_UNIQUE_DETID);
      GEMsimDig.push_back(gemdig);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "fpp1"){
      if(nparam_fpp1_read!=nparam_gemdet){
	cout << detectors_list[k] <<  " does not have the right number of parameters!!! " << endl << " fix database and retry! " << endl;
	exit(-1);
      }
      
      SBSDigGEMDet* fpp1 = new SBSDigGEMDet(FPP1_UNIQUE_DETID, NPlanes_fpp1, layer_fpp1, nstrips_fpp1, offset_fpp1, RO_angle_fpp1, 6, ZsupThr_fpp1);
      SBSDigGEMSimDig* gemdig = new SBSDigGEMSimDig(NPlanes_fpp1/2, triggeroffset_fpp1, gain_fpp1, ZsupThr_fpp1, nAPV_fpp1, commonmode_array_fpp1);
      for(int m = 0; m<Nlayers_fpp1; m++){
	fpp1->fZLayer.push_back(fpp1_layer_z[m]);
      }
      fpp1->fGateWidth = gatewidth_fpp1;
      
      GEMdetectors.push_back(fpp1);
      gemdetmap.push_back(FPP1_UNIQUE_DETID);
      GEMsimDig.push_back(gemdig);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "fpp2"){
      if(nparam_fpp2_read!=nparam_gemdet){
	cout << detectors_list[k] <<  " does not have the right number of parameters!!! " << endl << " fix database and retry! " << endl;
	exit(-1);
      }
      
      SBSDigGEMDet* fpp2 = new SBSDigGEMDet(FPP2_UNIQUE_DETID, NPlanes_fpp2, layer_fpp2, nstrips_fpp2, offset_fpp2, RO_angle_fpp2, 6, ZsupThr_fpp2);
      SBSDigGEMSimDig* gemdig = new SBSDigGEMSimDig(NPlanes_fpp2/2, triggeroffset_fpp2, gain_fpp2, ZsupThr_fpp2, nAPV_fpp2, commonmode_array_fpp2);
      for(int m = 0; m<Nlayers_fpp2; m++){
	fpp2->fZLayer.push_back(fpp2_layer_z[m]);
      }
      fpp2->fGateWidth = gatewidth_fpp2;
      
      GEMdetectors.push_back(fpp2);
      gemdetmap.push_back(FPP2_UNIQUE_DETID);
      GEMsimDig.push_back(gemdig);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "bbps"){
      if(nparam_bbps_read!=nparam_pmtdet_fadc_leadglass){
	cout << detectors_list[k] <<  " does not have the right number of parameters!!! " << endl << " fix database and retry! " << endl;
	exit(-1);
      }
      //SBSDigPMTDet* bbps = new SBSDigPMTDet(BBPS_UNIQUE_DETID, NChan_bbps, gain_bbps*qe, sigmapulse_bbps, gatewidth_bbps);
      SBSDigPMTDet* bbps = new SBSDigPMTDet(BBPS_UNIQUE_DETID, NChan_bbps, gain_bbps, nstages_bbps);
      bbps->fRefIndex = refindex_bbps;
      bbps->fGain = gain_bbps;
      bbps->fNstages = nstages_bbps;
      bbps->fPedestal = ped_bbps;
      bbps->fPedSigma = pedsigma_bbps;
      bbps->fTrigOffset = trigoffset_bbps;
      bbps->fThreshold = threshold_bbps*spe_unit/ROimpedance;
      bbps->fGateWidth = gatewidth_bbps;
      bbps->fADCconv = ADCconv_bbps;
      //bbps->fADCbits = ADCbits_bbps;
      bbps->fADCbits = FADC_ADCbits;
      bbps->fTDCconv = TDCconv_bbps;
      bbps->fTDCbits = TDCbits_bbps;
      bbps->fSigmaPulse = sigmapulse_bbps;
      bbps->SetSamples(FADC_sampsize);
      
      PMTdetectors.push_back(bbps);
      detmap.push_back(BBPS_UNIQUE_DETID);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "bbsh"){
      if(nparam_bbsh_read!=nparam_pmtdet_fadc_leadglass){
	cout << detectors_list[k] <<  " does not have the right number of parameters!!! " << endl << " fix database and retry! " << endl;
	exit(-1);
      }
      //SBSDigPMTDet* bbsh = new SBSDigPMTDet(BBSH_UNIQUE_DETID, NChan_bbsh, gain_bbsh*qe, sigmapulse_bbsh, gatewidth_bbsh);
      SBSDigPMTDet* bbsh = new SBSDigPMTDet(BBSH_UNIQUE_DETID, NChan_bbsh, gain_bbsh, nstages_bbsh);
      
      bbsh->fRefIndex = refindex_bbsh;
      bbsh->fGain = gain_bbsh;
      bbsh->fNstages = nstages_bbsh;
      bbsh->fPedestal = ped_bbsh;
      bbsh->fPedSigma = pedsigma_bbsh;
      bbsh->fTrigOffset = trigoffset_bbsh;
      bbsh->fThreshold = threshold_bbsh*spe_unit/ROimpedance;
      bbsh->fGateWidth = gatewidth_bbsh;
      bbsh->fADCconv = ADCconv_bbsh;
      bbsh->fADCbits = ADCbits_bbsh;    
      bbsh->fTDCconv = TDCconv_bbsh;
      bbsh->fTDCbits = TDCbits_bbsh;
      bbsh->fSigmaPulse = sigmapulse_bbsh;
      bbsh->SetSamples(FADC_sampsize);
            
      PMTdetectors.push_back(bbsh);
      detmap.push_back(BBSH_UNIQUE_DETID);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "grinch"){
      if(nparam_grinch_read!=nparam_pmtdet_fadc){
	cout << detectors_list[k] <<  " does not have the right number of parameters!!! " << endl << " fix database and retry! " << endl;
	exit(-1);
      }
      SBSDigPMTDet* grinch = new SBSDigPMTDet(GRINCH_UNIQUE_DETID, NChan_grinch, gain_grinch, nstages_grinch, sigmapulse_grinch, gatewidth_grinch);
  
      grinch->fGain = gain_grinch;
      grinch->fNstages = nstages_grinch;
      grinch->fPedestal = ped_grinch;
      grinch->fPedSigma = pedsigma_grinch;
      grinch->fTrigOffset = trigoffset_grinch;
      grinch->fThreshold = threshold_grinch*spe_unit/ROimpedance;
      grinch->fGateWidth = gatewidth_grinch;
      grinch->fADCconv = ADCconv_grinch;
      grinch->fADCbits = ADCbits_grinch;
      grinch->fTDCconv = TDCconv_grinch;
      grinch->fTDCbits = TDCbits_grinch;
      
      PMTdetectors.push_back(grinch);
      detmap.push_back(GRINCH_UNIQUE_DETID);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "bbhodo"){
      if(nparam_bbhodo_read!=nparam_pmtdet_fadc){
	cout << detectors_list[k] <<  " does not have the right number of parameters!!! " << endl << " fix database and retry! " << endl;
	exit(-1);
      }
      SBSDigPMTDet* bbhodo = new SBSDigPMTDet(HODO_UNIQUE_DETID, NChan_bbhodo, gain_bbhodo, nstages_bbhodo, sigmapulse_bbhodo, gatewidth_bbhodo);
      
      bbhodo->fGain = gain_bbhodo;
      bbhodo->fNstages = nstages_bbhodo;
      bbhodo->fPedestal = ped_bbhodo;
      bbhodo->fPedSigma = pedsigma_bbhodo;
      bbhodo->fTrigOffset = trigoffset_bbhodo;
      bbhodo->fThreshold = threshold_bbhodo*spe_unit/ROimpedance;
      bbhodo->fGateWidth = gatewidth_bbhodo;
      bbhodo->fADCconv = ADCconv_bbhodo;
      bbhodo->fADCbits = ADCbits_bbhodo;
      bbhodo->fTDCconv = TDCconv_bbhodo;
      bbhodo->fTDCbits = TDCbits_bbhodo; 
      
      PMTdetectors.push_back(bbhodo);
      detmap.push_back(HODO_UNIQUE_DETID);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "hcal"){
      if(nparam_hcal_read!=nparam_pmtdet_fadc_hitposresp){
	cout << detectors_list[k] <<  " does not have the right number of parameters!!! " << endl << " fix database and retry! " << endl;
	exit(-1);
      }
      SBSDigPMTDet* hcal = new SBSDigPMTDet(HCAL_UNIQUE_DETID, NChan_hcal, gain_hcal, nstages_hcal);
      
      hcal->fGain = gain_hcal;
      hcal->fNstages = nstages_hcal;
      hcal->fPedestal = ped_hcal;
      hcal->fPedSigma = pedsigma_hcal;
      hcal->fTrigOffset = trigoffset_hcal;
      hcal->fThreshold = threshold_hcal*spe_unit/ROimpedance;
      hcal->fGateWidth = gatewidth_hcal;
      hcal->fADCconv = ADCconv_hcal;
      hcal->fADCbits = FADC_ADCbits;
      hcal->fTDCconv = TDCconv_hcal;
      hcal->fTDCbits = TDCbits_hcal; 
      hcal->fSigmaPulse = sigmapulse_hcal; 
      hcal->SetSamples(FADC_sampsize);
      hcal->SetHitPosRespParam(3, hcalzhitresp_p);
      
      //ordered by increasing uinque id
      PMTdetectors.push_back(hcal);
      detmap.push_back(HCAL_UNIQUE_DETID);
      cout << " set up! " << endl;
    }

    if(detectors_list[k] == "ecal"){
      if(nparam_ecal_read!=nparam_pmtdet_fadc_leadglass){
	cout << detectors_list[k] <<  " does not have the right number of parameters!!! " << endl << " fix database and retry! " << endl;
	exit(-1);
      }
      SBSDigPMTDet* ecal = new SBSDigPMTDet(ECAL_UNIQUE_DETID, NChan_ecal, gain_ecal, nstages_ecal);
      
      ecal->fRefIndex = refindex_ecal;
      ecal->fGain = gain_ecal;
      ecal->fNstages = nstages_ecal;
      ecal->fPedestal = ped_ecal;
      ecal->fPedSigma = pedsigma_ecal;
      ecal->fTrigOffset = trigoffset_ecal;
      ecal->fThreshold = threshold_ecal*spe_unit/ROimpedance;
      ecal->fGateWidth = gatewidth_ecal;
      ecal->fADCconv = ADCconv_ecal;
      ecal->fADCbits = FADC_ADCbits;
      ecal->fTDCconv = TDCconv_ecal;
      ecal->fTDCbits = TDCbits_ecal; 
      ecal->fSigmaPulse = sigmapulse_ecal; 
      ecal->SetSamples(FADC_sampsize);
      
      //ordered by increasing uinque id
      PMTdetectors.push_back(ecal);
      detmap.push_back(ECAL_UNIQUE_DETID);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "cdet"){
      if(nparam_cdet_read!=nparam_pmtdet_fadc){
	cout << detectors_list[k] <<  " does not have the right number of parameters!!! " << endl << " fix database and retry! " << endl;
	exit(-1);
      }
      SBSDigPMTDet* cdet = new SBSDigPMTDet(CDET_UNIQUE_DETID, NChan_cdet, gain_cdet, nstages_cdet, sigmapulse_cdet, gatewidth_cdet);
      
      cdet->fGain = gain_cdet;
      cdet->fNstages = nstages_cdet;
      cdet->fPedestal = ped_cdet;
      cdet->fPedSigma = pedsigma_cdet;
      cdet->fTrigOffset = trigoffset_cdet;
      cdet->fThreshold = threshold_cdet*spe_unit/ROimpedance;
      cdet->fGateWidth = gatewidth_cdet;
      cdet->fADCconv = ADCconv_cdet;
      cdet->fADCbits = FADC_ADCbits;
      cdet->fTDCconv = TDCconv_cdet;
      cdet->fTDCbits = TDCbits_cdet; 
      
      //ordered by increasing uinque id
      PMTdetectors.push_back(cdet);
      detmap.push_back(CDET_UNIQUE_DETID);
      cout << " set up! " << endl;
    }
    
    // ** How to add a new subsystem **
    // Add the new detector here!
    if(detectors_list[k] == "prpolscint_bs"){
      if(nparam_prpolscint_bs_read!=nparam_pmtdet_fadc){
	cout << detectors_list[k] <<  " does not have the right number of parameters!!! " << endl << " fix database and retry! " << endl;
	exit(-1);
      }
      SBSDigPMTDet* polscint_bs = new SBSDigPMTDet(PRPOLBS_SCINT_UNIQUE_DETID, NChan_polscint_bs, gain_polscint_bs, nstages_polscint_bs, sigmapulse_polscint_bs, gatewidth_polscint_bs);
      
      polscint_bs->fGain = gain_polscint_bs;
      polscint_bs->fNstages = nstages_polscint_bs;
      polscint_bs->fPedestal = ped_polscint_bs;
      polscint_bs->fPedSigma = pedsigma_polscint_bs;
      polscint_bs->fTrigOffset = trigoffset_polscint_bs;
      polscint_bs->fThreshold = threshold_polscint_bs*spe_unit/ROimpedance;
      polscint_bs->fGateWidth = gatewidth_polscint_bs;
      polscint_bs->fADCconv = ADCconv_polscint_bs;
      polscint_bs->fADCbits = ADCbits_polscint_bs;
      polscint_bs->fTDCconv = TDCconv_polscint_bs;
      polscint_bs->fTDCbits = TDCbits_polscint_bs; 
      
      PMTdetectors.push_back(polscint_bs);
      detmap.push_back(PRPOLBS_SCINT_UNIQUE_DETID);
      cout << " set up! " << endl;
    } 
    
    if(detectors_list[k] == "prpolscint_fs"){
      if(nparam_prpolscint_fs_read!=nparam_pmtdet_fadc){
	cout << detectors_list[k] <<  " does not have the right number of parameters!!! " << endl << " fix database and retry! " << endl;
	exit(-1);
      }
      SBSDigPMTDet* polscint_fs = new SBSDigPMTDet(PRPOLFS_SCINT_UNIQUE_DETID, NChan_polscint_fs, gain_polscint_fs, nstages_polscint_fs, sigmapulse_polscint_fs, gatewidth_polscint_fs);
      
      polscint_fs->fGain = gain_polscint_fs;
      polscint_fs->fNstages = nstages_polscint_fs;
      polscint_fs->fPedestal = ped_polscint_fs;
      polscint_fs->fPedSigma = pedsigma_polscint_fs;
      polscint_fs->fTrigOffset = trigoffset_polscint_fs;
      polscint_fs->fThreshold = threshold_polscint_fs*spe_unit/ROimpedance;
      polscint_fs->fGateWidth = gatewidth_polscint_fs;
      polscint_fs->fADCconv = ADCconv_polscint_fs;
      polscint_fs->fADCbits = ADCbits_polscint_fs;
      polscint_fs->fTDCconv = TDCconv_polscint_fs;
      polscint_fs->fTDCbits = TDCbits_polscint_fs; 
      
      PMTdetectors.push_back(polscint_fs);
      detmap.push_back(PRPOLFS_SCINT_UNIQUE_DETID);
      cout << " set up! " << endl;
    } 

    if(detectors_list[k] == "activeana"){
      if(nparam_activeana_read!=nparam_pmtdet_fadc){
	cout << detectors_list[k] <<  " does not have the right number of parameters!!! " << endl << " fix database and retry! " << endl;
	exit(-1);
      }
      //SBSDigPMTDet* activeana = new SBSDigPMTDet(ACTIVEANA_UNIQUE_DETID, NChan_activeana, gain_activeana*qe, sigmapulse_activeana, gatewidth_activeana);
      SBSDigPMTDet* activeana = new SBSDigPMTDet(ACTIVEANA_UNIQUE_DETID, NChan_activeana, gain_activeana, nstages_activeana);
      
      activeana->fGain = gain_activeana;
      activeana->fNstages = nstages_activeana;
      activeana->fPedestal = ped_activeana;
      activeana->fPedSigma = pedsigma_activeana;
      activeana->fTrigOffset = trigoffset_activeana;
      activeana->fThreshold = threshold_activeana*spe_unit/ROimpedance;
      activeana->fGateWidth = gatewidth_activeana;
      activeana->fADCconv = ADCconv_activeana;
      activeana->fADCbits = ADCbits_activeana;
      activeana->fTDCconv = TDCconv_activeana;
      activeana->fTDCbits = TDCbits_activeana; 
      activeana->fSigmaPulse = sigmapulse_activeana; 
      activeana->SetSamples(FADC_sampsize);
      
      PMTdetectors.push_back(activeana);
      detmap.push_back(ACTIVEANA_UNIQUE_DETID);
      cout << " set up! " << endl;
    } 
  }
  
  TRandom3* R = new TRandom3(Rseed);
  
  // Step 1: read input files build the input chains
  // build signal chain
  ifstream sig_inputfile(inputsigfile);
  TChain *C_s = new TChain("T");
  while( currentline.ReadLine(sig_inputfile) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C_s->Add(currentline.Data());
    }
  }
  TObjArray *fileElements_s=C_s->GetListOfFiles();
  TIter next_s(fileElements_s);
  TChainElement *chEl_s=0;

  //restart reading to get background info
  while( currentline.ReadLine(sig_inputfile) && !currentline.BeginsWith("end_bkgdinfo")){
    if( !currentline.BeginsWith("#") ){
      Int_t ntokens = 0;
      std::unique_ptr<TObjArray> tokens( currentline.Tokenize(", \t") );
      if( !tokens->IsEmpty() ) {
	ntokens = tokens->GetLast()+1;
      }
      if(ntokens>=3){
	inputbkgdfile = ( (TObjString*) (*tokens)[0] )->GetString();
	TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	// EPAF: this "bkgd time window" notion is confusing. 
	// I will replace it with more straightforward stuff:
	// number of events total and current: 
	// then I will convert this to a to make it transparent to the code
	double nevtotal = stemp.Atof();
	stemp = ( (TObjString*) (*tokens)[2] )->GetString();
	double Ibeam = stemp.Atof();//in uA
	BkgdTimeWindow = nevtotal*qe/Ibeam/spe_unit;// in ns!
	LumiFrac = 1.0;
	if(Ibeam<=0)LumiFrac=0.0;
	if(ntokens==4){
	  stemp = ( (TObjString*) (*tokens)[3] )->GetString();
	  pmtbkgddig = stemp.Atoi();
	}
      }
      
    }
  }
  
  TFile* f_bkgd;
  SBSDigBkgdGen* BkgdGenerator;
  if(LumiFrac>0){
    f_bkgd = TFile::Open(inputbkgdfile.c_str());
    if(f_bkgd->IsZombie()){
      LumiFrac = 0;
    }else{
      BkgdGenerator = new SBSDigBkgdGen(f_bkgd, detectors_list, BkgdTimeWindow, pmtbkgddig);
      cout << "Includes background from file: " << inputbkgdfile.c_str() 
	   << " (integrated on " << BkgdTimeWindow << " ns time window);" 
	   << endl << " assuming " << LumiFrac*100 << "% luminosity." << endl
	   << " Doing PMT background digitization? " << pmtbkgddig << endl;
    }
  }
  
  double Theta_SBS, D_HCal;
  
  ULong64_t Nev_fs;
  ULong64_t ev_s;
  
  ULong64_t NEventsTotal = 0;
  
  int i_fs = 0;
  bool has_data;
  
  double timeZero;
  
  std::chrono::time_point<std::chrono::steady_clock> start = 
    std::chrono::steady_clock::now();
  
  while (( chEl_s=(TChainElement*)next_s() )) {
    if(NEventsTotal>=Nentries){
      break;
    }
    TFile f_s(chEl_s->GetTitle(), "UPDATE");
    if(f_s.IsZombie())cout << "File " << chEl_s->GetTitle() << " cannot be found. Please check the path of your file." << endl; 
    //run_data = (G4SBSRunData*)f_s.Get("run_data");
    G4SBSRunData* run_data = (G4SBSRunData*)f_s.Get("run_data");
    if(run_data==0){
      cout << "File does not have run data available!!! skip!" << endl;
      continue;
    }
    Theta_SBS = run_data->fSBStheta;
    D_HCal = run_data->fHCALdist;
    C_s = (TChain*)f_s.Get("T");
    g4sbs_tree *T_s = new g4sbs_tree(C_s, detectors_list, bool(LumiFrac));
    
    Nev_fs = C_s->GetEntries();
    
    for(ev_s = 0; ev_s<Nev_fs; ev_s++, NEventsTotal++){
      if(NEventsTotal>=Nentries)break;
      if(NEventsTotal%100==0)
	cout << NEventsTotal << "/" << Nentries << endl;
      
      timeZero = R->Gaus(0.0, TriggerJitter);
      
      //Clear detectors
      for(int k = 0; k<PMTdetectors.size(); k++){
	if(detmap[k]==HCAL_UNIQUE_DETID || detmap[k]==ECAL_UNIQUE_DETID || detmap[k]==BBPS_UNIQUE_DETID || detmap[k]==BBSH_UNIQUE_DETID || detmap[k]==ACTIVEANA_UNIQUE_DETID){
	  PMTdetectors[k]->Clear(true);
	}else{
	  PMTdetectors[k]->Clear();
	}
      }
      for(int k = 0; k<GEMdetectors.size(); k++){
	GEMdetectors[k]->Clear();
      }
      
      has_data = false;
      
      T_s->ClearDigBranches();
      T_s->GetEntry(ev_s);
      
      // unfold the thing then... but where???
      has_data = UnfoldData(T_s, Theta_SBS, D_HCal, R, PMTdetectors, detmap, GEMdetectors, gemdetmap, timeZero, 0);
      if(!has_data)continue;
      
      // if we want to add background, add background
      if(LumiFrac>0){
	//first digitize signal only...
	//	for(int k = 0; k<GEMdetectors.size(); k++){
	//  GEMsimDig[k]->Digitize(GEMdetectors[k], R);
	//  GEMsimDig[k]->CheckOut(GEMdetectors[k], gemdetmap[k], R, T_s, bool(LumiFrac>0));
	//}
	BkgdGenerator->GenerateBkgd(R, PMTdetectors, detmap, GEMdetectors, gemdetmap, LumiFrac);
      }
      
      //Digitize: PMT detectors
      for(int k = 0; k<PMTdetectors.size(); k++){
	PMTdetectors[k]->Digitize(T_s,R);
      }
      //digitize: GEMs
      for(int k = 0; k<GEMdetectors.size(); k++){
	GEMsimDig[k]->Digitize(GEMdetectors[k], R, bool(LumiFrac>0));
	GEMsimDig[k]->CheckOut(GEMdetectors[k], gemdetmap[k], R, T_s);
      }
      T_s->FillDigBranches();
    }// end loop on signal events 
    // if there are debugging histos to write, write them...
    for(int k = 0; k<GEMdetectors.size(); k++){
      cout << "GEM det ID: " << GEMdetectors[k]->fUniqueID << endl;
      GEMsimDig[k]->write_histos(GEMdetectors[k]);
      GEMsimDig[k]->print_time_execution();
    }
    if(LumiFrac>0)BkgdGenerator->WriteXCHistos();
    //write expanded tree
    T_s->fChain->Write("", TObject::kOverwrite);
    f_s.Write();
    f_s.Close();
    i_fs++;
  }// end loop on signal files
  
  std::chrono::time_point<std::chrono::steady_clock> end = 
    std::chrono::steady_clock::now();
  
  std::chrono::duration<double> diff = end-start;
  cout << " Total time " << std::setprecision(9) << diff.count() << " s "<< endl;
  
  exit(0);
}
