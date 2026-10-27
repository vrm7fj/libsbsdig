// fitCM.C
// Crystal Ball fits of the per-APV CM distributions.

#include <TFile.h>
#include <TH1D.h>
#include <TH2.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TTree.h>
#include <TString.h>
#include <TROOT.h>
#include <TMath.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>

static double kLoQuantile  = 0.001;
static double kMinBelow    = 10.;
static double kMaxQuantile = 0.2;
static double kMinTailSig  = 3.0;
static double kNMin        = 1.5;

// Returns fit status (0 = OK), or -1 if not fitted.
static int FitOneAPV(TH1D *h, const char *fname, TF1 *&fcb, double &fitLo, double &fitHi,
                     int &nBelow, int &nAbove, int minEntries)
{
  fcb = nullptr;
  nBelow = 0;
  nAbove = 0;
  const double nent = h->Integral();
  if (nent < minEntries) return -1;

  const int nb = h->GetNbinsX();
  const int bmax = h->GetMaximumBin();
  const double hmax = h->GetBinContent(bmax);
  const double peak = h->GetBinCenter(bmax);
  const double bw = h->GetBinWidth(bmax);

  int bl = bmax, br = bmax;
  while (bl > 1 && h->GetBinContent(bl) > 0.5 * hmax) bl--;
  while (br < nb && h->GetBinContent(br) > 0.5 * hmax) br++;
  const double sig0 = std::max((h->GetBinCenter(br) - h->GetBinCenter(bl)) / 2.355, bw);

  double p = std::max(kLoQuantile, kMinBelow / nent);
  p = std::min(p, kMaxQuantile);
  double q = 0.;
  h->GetQuantiles(1, &q, &p);
  fitLo = std::min(q, peak - kMinTailSig * sig0);
  fitLo = h->GetBinLowEdge(h->FindBin(fitLo));
  fitHi = peak + 4.0 * sig0;

  nBelow = (int)h->Integral(1, h->FindBin(fitLo) - 1);
  nAbove = (int)h->Integral(h->FindBin(fitHi) + 1, nb);

  fcb = new TF1(fname, "crystalball", fitLo, fitHi);
  fcb->SetParNames("Constant", "Mean", "Sigma", "Alpha", "N");
  fcb->SetParameters(hmax, peak, sig0, 1.0, 5.0);
  fcb->SetParLimits(0, 0., 10. * hmax);
  fcb->SetParLimits(1, peak - 3.0 * sig0, peak + 3.0 * sig0);
  fcb->SetParLimits(2, 0.2 * sig0, 5. * sig0);
  fcb->SetParLimits(3, 0.05, 10.);
  fcb->SetParLimits(4, kNMin, 500.);
  fcb->SetLineColor(kRed);
  fcb->SetLineWidth(2);
  fcb->SetNpx(1000);

  TFitResultPtr r = h->Fit(fcb, "LRQS0");
  int status = int(r);

  if (status != 0) {
    fcb->SetParameters(hmax, peak, sig0, 2.0, 2.0);
    r = h->Fit(fcb, "LRQS0");
    status = int(r);
  }
  return status;
}

static double Median(TH1D *h)
{
  double q, p = 0.5;
  h->GetQuantiles(1, &q, &p);
  return q;
}

static void AddShifted(TH1D *dst, const TH1D *src, int shiftBins)
{
  const int nb = src->GetNbinsX();
  for (int b = 1; b <= nb; b++) {
    const double c = src->GetBinContent(b);
    if (c == 0) continue;
    const int tb = b + shiftBins;
    if (tb < 1 || tb > nb) continue;
    dst->AddBinContent(tb, c);
  }
}

struct APVResult {
  TH1D *hp = nullptr;
  TH1D *htpl = nullptr;
  TF1 *fcb = nullptr;
  double fitLo = 0., fitHi = 0.;
  int nBelow = 0;
  int nAbove = 0;
  int status = -1;
  int ownStatus = -1;
  int source = -1;          // 0 own, 1 neighbours, 2 plane-wide, -1 none
  int nNeigh = 0;
  int neighLo = -1, neighHi = -1;
  double mu_t = 0.;
};

void fitCM(const char *rootfile, const char *det = "sbs_gemFT",
           const char *outprefix = "cmfit", bool logy = false,
           int minEntries = 50, int rebin = 1,
           bool useNeighbours = true,
           int maxDist = 2,
           int minEntriesMean = 10,
           bool replaceFailed = false)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile *fin = TFile::Open(rootfile, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "ERROR: cannot open " << rootfile << std::endl;
    return;
  }

  const TString pdfname  = TString::Format("%s.pdf", outprefix);
  const TString txtname  = TString::Format("%s.txt", outprefix);
  const TString rootname = TString::Format("%s.root", outprefix);

  std::ofstream txt(txtname.Data());
  txt << "# Crystal Ball fits of per-APV CM distributions from " << rootfile << "\n";
  txt << "# source: hCommonModeMean_by_APV_{U,V}_" << det << "_m<N>; digitizer plane = 2*module + axis (U=0, V=1)\n";
  txt << "# fit_lo = p-quantile, p = max(" << kLoQuantile << ", " << kMinBelow << "/entries) (<= " << kMaxQuantile
      << "), at most peak - " << kMinTailSig << " sigma0; fit_hi = peak + 4 sigma0; n >= " << kNMin << "\n";
  txt << "# module plane apv fit_lo fit_hi const mean sigma alpha n\n";

  TFile *fout = TFile::Open(rootname, "RECREATE");
  TTree *tree = new TTree("cmfit", "Crystal Ball fits of per-APV CM distributions");
  int t_module, t_axis, t_plane, t_apv, t_entries, t_below, t_above, t_ndf, t_source, t_nneigh;
  double t_par[5], t_err[5], t_chi2, t_fitlo, t_fithi;
  tree->Branch("module", &t_module);
  tree->Branch("axis", &t_axis);
  tree->Branch("plane", &t_plane);
  tree->Branch("apv", &t_apv);
  tree->Branch("entries", &t_entries);
  tree->Branch("below_range", &t_below);
  tree->Branch("above_range", &t_above);
  tree->Branch("fit_lo", &t_fitlo);
  tree->Branch("fit_hi", &t_fithi);
  tree->Branch("par", t_par, "par[5]/D");
  tree->Branch("err", t_err, "err[5]/D");
  tree->Branch("chi2", &t_chi2);
  tree->Branch("ndf", &t_ndf);
  tree->Branch("source", &t_source);
  tree->Branch("n_neighbours", &t_nneigh);

  TCanvas *c = new TCanvas("c", "CM fits", 1600, 1000);
  c->Print(pdfname + "[");

  const char *axisName[2] = {"U", "V"};
  int nPages = 0, nFitted = 0, nSkipped = 0, nFromNeigh = 0, nFromPlane = 0;

  for (int mod = 0; mod < 100; mod++) {
    bool anyAxis = false;

    for (int axis = 0; axis < 2; axis++) {
      TString hname = TString::Format("hCommonModeMean_by_APV_%s_%s_m%d", axisName[axis], det, mod);
      TH2 *h2 = dynamic_cast<TH2*>(fin->Get(hname));
      if (!h2) continue;
      anyAxis = true;

      const int plane = 2 * mod + axis;
      const int nAPV = h2->GetNbinsX();
      const int ncol = (int)std::ceil(std::sqrt((double)nAPV));
      const int nrow = (int)std::ceil((double)nAPV / ncol);

      TDirectory *dir = fout->mkdir(TString::Format("plane%d", plane),
                                    TString::Format("module m%d, %s strips", mod, axisName[axis]));

      // Pass 1: fit every APV on its own data
      std::vector<APVResult> res(nAPV);
      for (int apv = 0; apv < nAPV; apv++) {
        dir->cd();
        APVResult &R = res[apv];
        R.hp = h2->ProjectionY(TString::Format("hCM_apv%d", apv), apv + 1, apv + 1);
        if (rebin > 1) R.hp->Rebin(rebin);
        R.hp->SetTitle(TString::Format("%s APV %d;CM (ADC);entries", hname.Data(), apv));
        R.hp->SetLineColor(kBlue + 1);
        R.hp->SetFillColorAlpha(kAzure - 9, 0.5);

        R.ownStatus = FitOneAPV(R.hp, TString::Format("cb_plane%d_apv%d", plane, apv),
                                R.fcb, R.fitLo, R.fitHi, R.nBelow, R.nAbove, minEntries);
        R.status = R.ownStatus;
        R.source = (R.ownStatus == -1) ? -1 : 0;
      }

      // Pass 2: neighbour templates for low-stat APVs
      if (useNeighbours) {
        for (int apv = 0; apv < nAPV; apv++) {
          APVResult &R = res[apv];
          const bool lowStat = (R.ownStatus == -1);
          const bool failed  = (R.ownStatus > 0);
          if (!lowStat && !(replaceFailed && failed)) continue;

          auto collect = [&](int dist) {
            std::vector<int> nb;
            for (int j = 0; j < nAPV; j++) {
              if (j == apv || res[j].ownStatus != 0) continue;
              if (dist >= 0 && std::abs(j - apv) > dist) continue;
              nb.push_back(j);
            }
            return nb;
          };
          auto sumEntries = [&](const std::vector<int> &nb) {
            double s = 0.;
            for (int j : nb) s += res[j].hp->Integral();
            return s;
          };

          std::vector<int> nbs = collect(maxDist);
          int source = 1;
          if (nbs.empty() || sumEntries(nbs) < minEntries) {
            nbs = collect(-1);
            source = 2;
          }
          if (nbs.empty() || sumEntries(nbs) < minEntries) continue;

          double muAvg = 0., medOff = 0.;
          for (int j : nbs) {
            const double mu_j = res[j].fcb->GetParameter(1);
            muAvg  += mu_j;
            medOff += Median(res[j].hp) - mu_j;
          }
          muAvg  /= nbs.size();
          medOff /= nbs.size();
          const double ownEntries = R.hp->Integral();
          const double mu_t = (ownEntries >= minEntriesMean) ? Median(R.hp) - medOff : muAvg;

          dir->cd();
          TH1D *htpl = (TH1D *)R.hp->Clone(TString::Format("hCM_apv%d_template", apv));
          htpl->Reset();
          htpl->SetTitle(TString::Format("%s APV %d, template from neighbours;CM (ADC);entries",
                                         hname.Data(), apv));
          const double bw = htpl->GetBinWidth(1);
          for (int j : nbs) {
            const int shift = TMath::Nint((mu_t - res[j].fcb->GetParameter(1)) / bw);
            AddShifted(htpl, res[j].hp, shift);
          }
          htpl->Sumw2(kFALSE);
          htpl->Sumw2(kTRUE);
          htpl->SetEntries(htpl->Integral());

          TF1 *ftpl = nullptr;
          double lo = 0., hi = 0.;
          int nBelowTpl = 0, nAboveTpl = 0;
          const int st = FitOneAPV(htpl, TString::Format("cb_plane%d_apv%d_tpl", plane, apv),
                                   ftpl, lo, hi, nBelowTpl, nAboveTpl, minEntries);
          if (!ftpl) { delete htpl; continue; }

          if (ownEntries > 0) {
            const double scale = ownEntries / htpl->Integral();
            ftpl->SetParameter(0, ftpl->GetParameter(0) * scale);
            ftpl->SetParError(0, ftpl->GetParError(0) * scale);
          }
          ftpl->SetLineColor(kGreen + 2);
          ftpl->SetLineStyle(2);

          R.fcb    = ftpl;
          R.htpl   = htpl;
          R.fitLo  = lo;
          R.fitHi  = hi;
          R.status = st;
          R.source = source;
          R.nNeigh = nbs.size();
          R.neighLo = *std::min_element(nbs.begin(), nbs.end());
          R.neighHi = *std::max_element(nbs.begin(), nbs.end());
          R.mu_t = mu_t;
        }
      }

      // Pass 3: write and draw
      c->Clear();
      c->cd();
      TPad *ptitle = new TPad("ptitle", "", 0., 0.95, 1., 1.);
      TPad *pmain  = new TPad("pmain",  "", 0., 0.,   1., 0.95);
      ptitle->Draw();
      pmain->Draw();
      ptitle->cd();
      TLatex title;
      title.SetNDC();
      title.SetTextSize(0.55);
      title.SetTextAlign(22);
      title.DrawLatex(0.5, 0.5, TString::Format("%s  module m%d,  %s strips  (digitizer plane %d)   %s",
                                                det, mod, axisName[axis], plane, rootfile));
      pmain->Divide(ncol, nrow, 0.002, 0.002);

      for (int apv = 0; apv < nAPV; apv++) {
        APVResult &R = res[apv];
        TH1D *hp = R.hp;
        TF1 *fcb = R.fcb;

        dir->cd();
        hp->Write();
        if (R.htpl) R.htpl->Write();
        if (R.source > 0 && fcb) fcb->Write();

        t_module = mod; t_axis = axis; t_plane = plane; t_apv = apv;
        t_entries = (int)hp->Integral();
        t_below = R.nBelow;
        t_above = R.nAbove;
        t_fitlo = R.fitLo;
        t_fithi = R.fitHi;
        t_source = R.source;
        t_nneigh = R.nNeigh;
        for (int k = 0; k < 5; k++) {
          t_par[k] = fcb ? fcb->GetParameter(k) : 0.;
          t_err[k] = fcb ? fcb->GetParError(k) : 0.;
        }
        t_chi2 = fcb ? fcb->GetChisquare() : 0.;
        t_ndf  = fcb ? fcb->GetNDF() : 0;
        tree->Fill();

        if (R.status == 0 && fcb) {
          txt << mod << " " << plane << " " << apv << " "
              << R.fitLo << " " << R.fitHi << " "
              << t_par[0] << " " << t_par[1] << " " << t_par[2] << " "
              << t_par[3] << " " << t_par[4] << "\n";
        }

        if (R.ownStatus == -1) nSkipped++;
        else if (R.ownStatus == 0) nFitted++;
        if (R.source == 1) nFromNeigh++;
        if (R.source == 2) nFromPlane++;

        TVirtualPad *pad = pmain->cd(apv + 1);
        pad->SetLeftMargin(0.14);
        pad->SetRightMargin(0.03);
        pad->SetTopMargin(0.03);
        pad->SetBottomMargin(0.14);

        const bool ownHasData = hp->Integral() > 0;
        TH1D *hdraw = ownHasData ? hp : R.htpl;
        if (R.htpl) {
          R.htpl->SetLineColor(kGreen + 2);
          R.htpl->SetFillColorAlpha(kGreen - 9, 0.4);
        }
        pad->SetLogy(logy && hdraw && hdraw->Integral() > 0);

        if (hdraw && hdraw->Integral() > 0) {
          double xlo, xhi;
          if (fcb) {
            const double sig = std::fabs(fcb->GetParameter(2));
            xlo = R.fitLo - 5. * hdraw->GetBinWidth(1);
            xhi = fcb->GetParameter(1) + 10. * sig;
          } else {
            double q[2], p[2] = {0.001, 0.999};
            hdraw->GetQuantiles(2, q, p);
            xlo = q[0] - 20.;
            xhi = q[1] + 20.;
          }
          hdraw->GetXaxis()->SetRangeUser(xlo - 200, xhi);
          if (logy) hdraw->SetMinimum(0.5);
          hdraw->GetXaxis()->SetLabelSize(0.06);
          hdraw->GetYaxis()->SetLabelSize(0.06);
          hdraw->GetXaxis()->SetNdivisions(505);
          hdraw->Draw("HIST");
          if (fcb) fcb->Draw("SAME");
        } else {
          hp->Draw("AXIS");
        }

        TLatex lat;
        lat.SetNDC();
        lat.SetTextSize(0.075);
        lat.SetTextAlign(13);
        lat.DrawLatex(0.17, 0.95, TString::Format("APV %d  (%d)", apv, t_entries));
        if (R.source == -1) {
          lat.SetTextColor(kGray + 2);
          lat.DrawLatex(0.17, 0.86, "not fitted: too few entries");
          lat.DrawLatex(0.17, 0.78, "(no usable neighbours)");
        } else {
          if (R.source > 0) lat.SetTextColor(kGreen + 3);
          lat.DrawLatex(0.17, 0.86, TString::Format("#mu=%.0f  #sigma=%.1f", t_par[1], t_par[2]));
          lat.DrawLatex(0.17, 0.78, TString::Format("#alpha=%.2f  n=%.1f", t_par[3], t_par[4]));
          if (R.source == 0) {
            lat.DrawLatex(0.17, 0.70, TString::Format("#chi^{2}/ndf=%.0f/%d", t_chi2, t_ndf));
          } else {
            lat.DrawLatex(0.17, 0.70, TString::Format("from %d %s (APV %d-%d)", R.nNeigh,
                                                      R.source == 1 ? "neighbours" : "plane APVs",
                                                      R.neighLo, R.neighHi));
          }
        }
      }

      c->Print(pdfname, TString::Format("Title:m%d %s", mod, axisName[axis]));
      nPages++;
    }

    if (!anyAxis && mod > 0) break;
  }

  c->Print(pdfname + "]");

  fout->cd();
  tree->Write();
  fout->Close();
  txt.close();
  fin->Close();

  if (nPages == 0) {
    std::cerr << "ERROR: no hCommonModeMean_by_APV_{U,V}_" << det << "_m<N> histograms found in " << rootfile << std::endl;
    return;
  }
  std::cout << "Pages: " << nPages << "   APVs fitted OK: " << nFitted
            << "   skipped (< " << minEntries << " entries): " << nSkipped << "\n"
            << "Generated from neighbours (+-" << maxDist << "): " << nFromNeigh
            << "   from whole plane: " << nFromPlane << "\n"
            << "Wrote " << pdfname << ", " << txtname << ", " << rootname << std::endl;
}
