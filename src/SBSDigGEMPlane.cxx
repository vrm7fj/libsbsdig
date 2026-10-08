#include "SBSDigGEMPlane.h"
#include "TMath.h"

using namespace std;

//
// Class SBSDigGEMPlane
//
SBSDigGEMPlane::SBSDigGEMPlane() :
  fLayer(0), fNStrips(3840), fNSamples(6), fStripThr(100), fXoffset(0), fROangle(0)
{  
  // SetPlaneStripPed();
  // SetPlaneAPVCM();
}

SBSDigGEMPlane::SBSDigGEMPlane(short trackeruniqueid, short layer, short mod, int nstrips, int nsamples, double thr, double offset, double roangle) :
  fTrackerID(trackeruniqueid), fLayer(layer), fModule(mod), fNStrips(nstrips), fNSamples(nsamples), fStripThr(thr), fXoffset(offset), fROangle(roangle)
{
  //fModule = mod;
  fdX = fNStrips*4.e-4;
  
  fStripADCsum = new Int_t[fNStrips];
  fStripADC = new Short_t[fNStrips*fNSamples];
  fStripMipADCsum = new Int_t[fNStrips];
  fStripMipADC = new Short_t[fNStrips*fNSamples];
  fStripPedSubADC = new Short_t[fNStrips*fNSamples];
  fStripCMCorrADC = new Short_t[fNStrips*fNSamples];

  fNAPVs = int(fNStrips/fNChanAPV);

  TString trackerName = "undefined";

  if ( fTrackerID == 42 ) trackerName = "ft";
  else if ( fTrackerID == 43 ) trackerName = "fpp1";

  fCommonModeDiffOnlineMean = new TH2D( TString::Format("hDiffOnlineCMandCMmean_%s_mod%i",trackerName.Data(),fModule), "Online Danning CM - CM mean; APV card; CM online - CM mean (user)", fNAPVs, -0.5, fNAPVs-0.5, 2500, -1000.0, 1000.0 );
  fCommonModeDiffOnlineTrue = new TH2D( TString::Format("hDiffOnlineCMandTrueCM_%s_mod%i",trackerName.Data(),fModule), "Online Danning CM - True CM; APV card; CM online - True CM (simulated)", fNAPVs, -0.5, fNAPVs-0.5, 2500, -1000.0, 1000.0 );

  fTrueCMbyAPV   = new TH2D( TString::Format("hTrueCM_by_APV_%s_mod%i",trackerName.Data(),fModule),   "Generated (true) CM; APV card; true CM (ADC)",   fNAPVs, -0.5, fNAPVs-0.5, 1024, -0.5, 4095.5 );
  fOnlineCMbyAPV = new TH2D( TString::Format("hOnlineCM_by_APV_%s_mod%i",trackerName.Data(),fModule), "Online Danning CM; APV card; online CM (ADC)", fNAPVs, -0.5, fNAPVs-0.5, 1024, -0.5, 4095.5 );

  Clear();
}

SBSDigGEMPlane::~SBSDigGEMPlane()
{
  Clear();
}


void SBSDigGEMPlane::Clear()
{
  memset(fStripADCsum, 0, fNStrips*sizeof(Int_t));
  memset(fStripADC, 0, fNStrips*fNSamples*sizeof(Short_t));
  memset(fStripPedSubADC, 0, fNStrips*fNSamples*sizeof(Short_t));
  memset(fStripCMCorrADC, 0, fNStrips*fNSamples*sizeof(Short_t));
  memset(fStripMipADCsum, 0, fNStrips*sizeof(Int_t));
  memset(fStripMipADC, 0, fNStrips*fNSamples*sizeof(Short_t));
}

//ClassImp(SBSDigGEMPlane);

// 'Online' CM and ZS.

void SBSDigGEMPlane::DoPedSub(){ // First subtract pedestals from each channel before calculating online CM.

  for ( int istrip=0; istrip < fNStrips; istrip++ ){
    for ( int isamp=0; isamp < fNSamples; isamp++ ){
      fStripPedSubADC[istrip*fNSamples+isamp] = fStripADC[istrip*fNSamples+isamp]; //- fPedestalMap[istrip].mean; // Removed pedestal offset addition in SBSDigGEMSimDig.cxx on Feb 17, 2026. ADR.
    }    
  }
}

double SBSDigGEMPlane::GetOnlineCommonMode(const std::array <int,fNChanAPV>& adcarray, int apvnum){ // We are replicating 'falg 4' in SBSGEMModule::GetCommonMode() method here.
  
  double cm_mean = GetAPVCMMean(apvnum);
  double cm_rms = GetAPVCMSigma(apvnum); // Here I'm changing CM sigma to RMS to go with the SBS_OFFLINE convension. But shouldn't it be sigma as true CM values is not 0 ?

  double cm_temp = 0.;

  int n_keep = 0;

  for (int iter=0; iter<3; iter++){

    double cm_min = cm_mean - fCommonModeRange_nsigma*cm_rms;

    if (iter==0) cm_min = 0.;

    double cm_max = cm_mean + fCommonModeRange_nsigma*cm_rms;
    double sumADCinrange = 0.;
    n_keep = 0;

    for (int ichan=0; ichan<fNChanAPV; ichan++){//Size of the array will be 128.

      double ADCtemp = (double)adcarray[ichan];
      double pedrmstemp = fPedestalMap[apvnum*fNChanAPV+ichan].rms;

      if (iter != 0){
        cm_min = cm_temp - fCommonModeDanningMethod_NsigmaCut*2.5*pedrmstemp;
        cm_max = cm_temp + fCommonModeDanningMethod_NsigmaCut*2.5*pedrmstemp;
      }

      if ( ADCtemp >= cm_min && ADCtemp <= cm_max ){
        n_keep++;
        sumADCinrange += ADCtemp;
      }
    }

    if (n_keep == 0) return cm_mean;

    cm_temp = sumADCinrange / n_keep;
  }

  if( n_keep < fCommonModeMinStripsInRange ) return cm_mean;

  return cm_temp;

}

void SBSDigGEMPlane::ApplyOnlineCMCorr( const bool do_onlineCM_histos ){ // Calculate CM for per 128 strip-set (per APV) for each TS.
// Then apply the CM correction to the fStripADC array. 

  //DoPedSub();// First have to subtract pedestal offsets. --> We are NOT adding pedestal offsets now. ADR 02/19/2026.

  for ( int istrip=0; istrip < fNStrips; istrip += fNChanAPV ){

    int iAPV = int(istrip/fNChanAPV);

    for ( int isamp=0; isamp < fNSamples; isamp++ ){

      std::array <int, fNChanAPV> thisAPVthisSampPedSubADC{-1500};

      for ( int ichan=0; ichan < fNChanAPV; ichan++ ){
        //thisAPVthisSampPedSubADC[ichan] = fStripPedSubADC[(istrip+ichan)*fNSamples+isamp];
        thisAPVthisSampPedSubADC[ichan] = fStripADC[(istrip+ichan)*fNSamples+isamp];
      }

      double thisAPVthisSampOnineCM = GetOnlineCommonMode(thisAPVthisSampPedSubADC, iAPV);

      // if ( do_onlineCM_histos ){
      //   fCommonModeDiffOnlineMean->Fill( iAPV, thisAPVthisSampOnineCM - GetAPVCMMean(iAPV) );
      //   fCommonModeDiffOnlineTrue->Fill( iAPV, thisAPVthisSampOnineCM - fTrueCMmap.at(ApvTS{iAPV, isamp}) );
      // }

      if ( do_onlineCM_histos ){
        fCommonModeDiffOnlineMean->Fill( iAPV, thisAPVthisSampOnineCM - GetAPVCMMean(iAPV) );
        fCommonModeDiffOnlineTrue->Fill( iAPV, thisAPVthisSampOnineCM - fTrueCMmap.at(ApvTS{iAPV, isamp}) );
        fOnlineCMbyAPV->Fill( iAPV, thisAPVthisSampOnineCM );   // new
      }

      for ( int ichan=0; ichan < fNChanAPV; ichan++ ){
        //fStripCMCorrADC[(istrip+ichan)*fNSamples+isamp] = TMath::Nint(fStripPedSubADC[(istrip+ichan)*fNSamples+isamp] - thisAPVthisSampOnineCM);
        // fStripADC[(istrip+ichan)*fNSamples+isamp] = fStripCMCorrADC[(istrip+ichan)*fNSamples+isamp];
        fStripCMCorrADC[(istrip+ichan)*fNSamples+isamp] = TMath::Nint(fStripADC[(istrip+ichan)*fNSamples+isamp] - thisAPVthisSampOnineCM);
        SetADC((istrip+ichan), isamp, fStripCMCorrADC[(istrip+ichan)*fNSamples+isamp]);
      }
    }
  }  
}

void SBSDigGEMPlane::ApplyOnlineZS(const double zs_thr_nsigma){

  for ( int istrip=0; istrip < fNStrips; istrip++ ){

    bool stripPassZS = GetADCSum(istrip)/double(fNSamples) > zs_thr_nsigma*fPedestalMap[istrip].rms;

    // zero-suppress the strips that did not pass the threshold.
    if (!stripPassZS){
      for ( int isamp=0; isamp < fNSamples; isamp++ ){
        SetADC(istrip, isamp, 0);
      }
    }
  }
}

// void SBSDigGEMPlane::WriteCMHistos(){

//   fCommonModeDiffOnlineMean->Write();
//   fCommonModeDiffOnlineTrue->Write();
// }

void SBSDigGEMPlane::WriteCMHistos(){
  fCommonModeDiffOnlineMean->Write();
  fCommonModeDiffOnlineTrue->Write();
  fTrueCMbyAPV->Write();     // new
  fOnlineCMbyAPV->Write();   // new
}

// Build one Crystal Ball TF1 per APV from the beam-CM fit table.
void SBSDigGEMPlane::SetPlaneAPVCMFit(const cmfitmap& fits){
  for ( const auto& [apv, f] : fits ){
    TF1* cb = new TF1(Form("cmcb_t%d_p%d_a%d", fTrackerID, fModule, apv), "crystalball", f.lo, f.hi);
    cb->SetParameters(1., f.mean, f.sigma, f.alpha, f.n); // constant irrelevant for sampling
    cb->SetNpx(2000);
    fCMFitFunc[apv] = cb;
  }
}

// Draw a true CM value: Crystal Ball if fitted, else Gaussian from the CM file.
double SBSDigGEMPlane::SampleTrueCM(int apv, TRandom3* R){
  auto it = fCMFitFunc.find(apv);
  if ( it != fCMFitFunc.end() ) return it->second->GetRandom(R);
  return R->Gaus(GetAPVCMMean(apv), GetAPVCMSigma(apv));  // fallback: old behaviour
}
