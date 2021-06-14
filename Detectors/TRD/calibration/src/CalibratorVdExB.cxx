
// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CalibratorVdExB.cxx
/// \brief TimeSlot-based calibration of vDrift and ExB
/// \author Ole Schmidt

#include "TRDCalibration/CalibratorVdExB.h"
#include <memory>
//#include "TProfile.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TVirtualFitter.h"
#include "Math/MinimizerOptions.h"
#include "TMinuitMinimizer.h"
#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"
#include "TRDBase/CalVdriftExB.h"

using namespace std;
using namespace o2::ccdb;
using namespace o2::trd;
using namespace o2::trd::constants;

//-------------------------------------------------
static std::array<std::unique_ptr<TProfile>, o2::trd::constants::MAXCHAMBER> arr_tp_angleDiff;
static int iDet;
static double vD_pre_corr_set = 1.546; // 1.546, 1.2
static double LA_pre_corr_set = -0.16133;
static double low_angle_fit = 80.0 * TMath::DegToRad();
static double up_angle_fit = 100.0 * TMath::DegToRad();
//-------------------------------------------------

//----------------------------------------------------------------------------------------
TGraph* calc_Delta_alpha(double vD_fit, double Lorentz_angle_fit, double vD_pre_corr, double Lorentz_angle_pre_corr)
{
  // 06.10.2020 A. Schmah Added comments, replaced drift_vel_ratio by (vD_pre_corr/vD_fit)

  // The first cluster (the one at the bottom, last time time) is drifting to the anode plane. The position there is
  // independent of the true drift velocity and depends only on the Lorentz angle. BUT the time bin depends on the
  // true drift velocity. For large vDtrue the time bin is small, e.g. 18, for small vDtrue the time bin is large, e.g. 24.
  // The cluster is then projected with the fixed drift velocity, e.g. vD_pre_corr = 1.546, along local y. Therefore Delta alpha vs. impact angle
  // depends on vDtrue and vD_pre_corr and the Lorentz angle.

  TGraph* TG_Delta_alpha_vs_impact_angle = new TGraph();

  int i_point = 0;
  for (double impact_angle = 65.0 * TMath::DegToRad(); impact_angle < 115.0 * TMath::DegToRad(); impact_angle += 1.0 * TMath::DegToRad()) {

    double TRD_anode_plane = 0.0335;

    // Direction vector of incoming track
    double x_dir = TMath::Cos(impact_angle);
    double y_dir = TMath::Sin(impact_angle);

    // Slope of incoming track
    double slope = 10000000.0;
    if (x_dir != 0.0)
      slope = y_dir / x_dir;

    // Slope of Lorentz angle
    double Lorentz_tan = TMath::Tan(Lorentz_angle_fit);
    double Lorentz_slope = 10000000.0;
    if (Lorentz_tan != 0.0)
      Lorentz_slope = 1.0 / Lorentz_tan;

    // Hit point of incoming track with anode plane
    double x_anode_hit = TRD_anode_plane / slope;
    double y_anode_hit = TRD_anode_plane;

    // Hit point at anode plane of Lorentz angle shifted cluster from the entrance -> independent of vDtrue
    double x_Lorentz_anode_hit = TRD_anode_plane / Lorentz_slope;
    double y_Lorentz_anode_hit = TRD_anode_plane;

    // Cluster location within drift cell of cluster from entrance after drift velocity ratio is applied
    double x_Lorentz_drift_hit = x_Lorentz_anode_hit;
    double y_Lorentz_drift_hit = TRD_anode_plane - TRD_anode_plane * (vD_pre_corr / vD_fit);

    // Reconstructed hit of first cluster at chamber entrance after pre Lorentz angle correction
    double x_Lorentz_drift_hit_pre_corr = x_Lorentz_anode_hit - (TRD_anode_plane - y_Lorentz_drift_hit) * TMath::Tan(Lorentz_angle_pre_corr);
    double y_Lorentz_drift_hit_pre_corr = TRD_anode_plane - TRD_anode_plane * (vD_pre_corr / vD_fit);

    double impact_angle_track = TMath::ATan2(y_anode_hit, x_anode_hit);

    double Delta_x_Lorentz_drift_hit = x_anode_hit - x_Lorentz_drift_hit_pre_corr;
    double Delta_y_Lorentz_drift_hit = y_anode_hit - y_Lorentz_drift_hit_pre_corr;
    double impact_angle_rec = TMath::ATan2(Delta_y_Lorentz_drift_hit, Delta_x_Lorentz_drift_hit);

    double Delta_angle = -(impact_angle_track - impact_angle_rec);
    TG_Delta_alpha_vs_impact_angle->SetPoint(i_point, impact_angle * TMath::RadToDeg(), Delta_angle * TMath::RadToDeg());
    i_point++;
    //printf("impact_angle: %4.3f, impact_angle_track: %4.3f, impact_angle_rec: %4.3f, Delta: {%4.3f, %4.3f}, x_anode_hit: %4.3f, x_Lorentz_drift_hit: %4.3f \n",impact_angle*TMath::RadToDeg(),impact_angle_track*TMath::RadToDeg(),impact_angle_rec*TMath::RadToDeg(),Delta_x_Lorentz_drift_hit,Delta_y_Lorentz_drift_hit,x_anode_hit,x_Lorentz_drift_hit);
  }

  return TG_Delta_alpha_vs_impact_angle;
}
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
void Chi2_TRD_vDrift(int&, double*, double& sum, double* par, int)
{
  // function to be minimized
  sum = 0;

  //double B_field_use = par[0];
  //double E_field_use = par[1];
  //double v_drift_use = par[2];
  //double vD_use      = par[3];
  //double LA_use      = par[4];

  double LA_fit = par[0];
  double vD_fit = par[1];

  //double vD_fit, double Lorentz_angle_fit, double vD_pre_corr, double Lorentz_angle_pre_corr
  TGraph* tg_Delta_vs_impact_single = calc_Delta_alpha(vD_fit, LA_fit, vD_pre_corr_set, LA_pre_corr_set); // new version
  //TGraph* tg_Delta_vs_impact_single = calc_Delta_alpha(LA_use,vD_ratio_use); // old version
  //TGraph* tg_Delta_vs_impact_single = calc_Delta_alpha(LA_use,vD_ratio_use,Lorentz_angle_pre_corr); // old version

  for (int i_bin = 1; i_bin <= arr_tp_angleDiff[iDet]->GetNbinsX(); i_bin++) {
    //double impact_angle     = vec_tp_Delta_vs_impact_circle[iDet] ->GetBinCenter(i_bin);
    //double Delta_alpha      = vec_tp_Delta_vs_impact_circle[iDet] ->GetBinContent(i_bin);

    double impact_angle = (arr_tp_angleDiff[iDet]->GetBinCenter(i_bin) + 90.0) * TMath::DegToRad();
    double Delta_alpha = (arr_tp_angleDiff[iDet]->GetBinContent(i_bin)) * TMath::DegToRad();

    if (Delta_alpha == 0.0)
      continue;
    if (impact_angle < low_angle_fit)
      continue;
    if (impact_angle > up_angle_fit)
      continue;

    double Delta_alpha_sim = tg_Delta_vs_impact_single->Eval(impact_angle);

    //sum += TMath::Power(Delta_alpha_sim - Delta_alpha,2)/TMath::Power(Delta_alpha_err,2);

    //LOGF(INFO, "iDet: %i, impact_angle: %4.3f, Delta_alpha: %4.3f, Delta_alpha_sim: %4.3f",iDet,impact_angle,Delta_alpha,Delta_alpha_sim);
    sum += TMath::Power(Delta_alpha_sim - Delta_alpha, 2) / 1.0;
  }

  delete tg_Delta_vs_impact_single;
}
//----------------------------------------------------------------------------------------

namespace o2::trd
{

using Slot = o2::calibration::TimeSlot<AngularResidHistos>;

void CalibratorVdExB::initOutput()
{
  // prepare output objects which will go to CCDB

#if 1
  //Connect to CCDB
  //
  printf("Connect to CCDB for writing \n");
  Long64_t Run = 265501;
  o2::ccdb::CcdbApi ccdb;
  map<string, string> metadata;               // do we want to store any meta data?
  ccdb.init("http://ccdb-test.cern.ch:8080"); // or http://localhost:8080 for a local installation

  CalVdriftExB* o2CalVdriftExB = new o2::trd::CalVdriftExB();
  for (Int_t i_det = 0; i_det < 540; i_det++) {
    o2CalVdriftExB->setVDrift(i_det, arr_vD_fit[i_det]);
    o2CalVdriftExB->setExB(i_det, arr_LA_fit[i_det]);
  }

  // Write
  printf("Write to CCDB \n");
  ccdb.storeAsTFileAny(o2CalVdriftExB, "TRD_test/CalVdriftExB", metadata, Run, Run + 1);
  //------------------------------------------
#endif
}

void CalibratorVdExB::finalizeSlot(Slot& slot)
{
  // do actual calibration for the data provided in the given slot
  // TODO!

  LOG(INFO) << "WE ARE IN finalizeSlot";

  o2::trd::AngularResidHistos* residhistos = slot.getContainer();

  //for (int index = 0; index < residhistos->getNEntries(); index++)
  //{
  //    float angleDiff = residhistos->getHistogramEntry(index);
  //    int count = residhistos->getBinCount(index);
  //LOGF(INFO, "Check what we read: angleDiff = %4.3f, count = %d ", angleDiff, count);
  //}

  //----------------------------------------------------------------------
  // Detector loop
  for (iDet = 0; iDet < o2::trd::constants::MAXCHAMBER; ++iDet) {

    //arr_tp_angleDiff[iDet] = std::make_unique<TProfile>(constants::NBINSANGLEDIFF,-25.0,25.0); // xy distribution of TRD space points //moved to init

    //LOG(INFO) << "TEST 1 HERE";
    //if (mIsInitialized == false) arr_tp_angleDiff[iDet] = std::make_unique<TProfile>(Form("arr_tp_angleDiff_%d",iDet),Form("arr_tp_angleDiff_%d",iDet),constants::NBINSANGLEDIFF,-25.0,25.0); // xy distribution of TRD space points

    arr_tp_angleDiff[iDet] = std::make_unique<TProfile>(Form("arr_tp_angleDiff_%d", iDet), Form("arr_tp_angleDiff_%d", iDet), o2::trd::constants::NBINSANGLEDIFF, -25.0, 25.0); // xy distribution of TRD space points

    //LOG(INFO) << "TEST 2";
    //add array with RMSs - on the fly
    //second array with counters?

    //tracklets are calibrated or not?

    //int iPoint = 0;

    //----------------------------------------------------------------------
    // Fill the profiles
    for (int iBin = 0; iBin < o2::trd::constants::NBINSANGLEDIFF; ++iBin) { // note: iBin = constants::NBINSANGLEDIFF - 1 is under-/overflow bin

      //for (int iBin = 5; iBin < 20; ++iBin) {
      //read real values - enable later

      //LOG(INFO) << "TEST 3";
      float angleDiffSum = residhistos->getHistogramEntry(iDet * o2::trd::constants::NBINSANGLEDIFF + iBin);
      short nEntries = residhistos->getBinCount(iDet * o2::trd::constants::NBINSANGLEDIFF + iBin);
      //LOG(INFO) << "TEST 4";
      if (nEntries > 0) {
        //LOGF(INFO, "Found %i entrie(s) in chamber %i, bin %i. Average angular deviation: %f", nEntries, iDet, iBin, angleDiffSum / nEntries);
        //arr_tg_angleDiff[iDet] ->SetPoint(iBin,iBin,angleDiffSum/nEntries);
        //PUT BACK
        //arr_tp_angleDiff[iDet] ->SetPoint(iPoint,(2*iBin-25),angleDiffSum/nEntries); //double check bins

        //LOGF(INFO, "Fill: x = %f, y =  %f, w =  %f",(double)(2*iBin-25),angleDiffSum/nEntries,(double)nEntries);
        arr_tp_angleDiff[iDet]->Fill((double)(2 * iBin - 25), (double)angleDiffSum / nEntries, (double)nEntries); //double check bins
                                                                                                                  //arr_tp_angleDiff[iDet] ->Fill(iBin,2); //test

        //LOGF(INFO, "Fill: x = %f, y =  %f, w =  %f",(2.0*(double)iBin-25.0),angleDiffSum/nEntries,nEntries);
        //arr_tp_angleDiff[iDet] ->Fill((2.0*(double)iBin-25.0),angleDiffSum/nEntries,nEntries); //double check bins

        //iPoint++;
        //LOGF(INFO, "Point written: bin %i, iDet: %i, trkAng: %4.3f, angdiff: %f",iBin,iDet,arr_tp_angleDiff[iDet] ->GetBinCenter(iBin),arr_tp_angleDiff[iDet] ->GetBinContent(iBin));
      }

      //arr_tg_angleDiff[iDet] ->SetPoint(iBin,vec_tp_Delta_vs_impact_circle[iDet]->GetBinCenter(iBin+212),vec_tp_Delta_vs_impact_circle[iDet]->GetBinContent(iBin+212));
      //LOGF(INFO, "Point written: bin %i, x: %f, angdiff: %f",iBin,arr_tg_angleDiff[iDet] ->GetPointX(iBin),arr_tg_angleDiff[iDet] ->GetPointY(iBin));
    }
    //----------------------------------------------------------------------

    //----------------------------------------------------------------------
    //minimization
#if 1

    TVirtualFitter* min = TVirtualFitter::Fitter(0, 2);
    min->SetFCN(Chi2_TRD_vDrift);
    double pStart[2] = {-7.5 * TMath::DegToRad(), 1.0};

    min->SetParameter(0, "LA", pStart[0], 0.01, 0, 0);
    min->SetParameter(1, "vD", pStart[1], 0.01, 0, 0);

    double arglist[2];
    arglist[0] = 1000;  // number of function calls
    arglist[1] = 0.001; // tolerance
    //arglist[0] = -1.0;
    //min->ExecuteCommand("SET PRINT", arglist, 1);

    double arglist_A[1] = {-1};
    double arglist_B[1] = {0};
    double arglist_C[1] = {-1};

    min->ExecuteCommand("SET PRIntout", arglist_A, 1);
    min->ExecuteCommand("SET NOWarnings", arglist_B, 1);

    min->ExecuteCommand("MIGRAD", arglist, 0);

    double parFit[5];

    for (int i = 0; i < 2; ++i) {
      parFit[i] = min->GetParameter(i);
    }

    arr_LA_fit[iDet] = parFit[0];
    arr_vD_fit[iDet] = parFit[1];

    LOGF(INFO, "iDet: %i, arr_vD_fit[%i]: %4.3f, arr_LA_fit[%i]: %4.3f", iDet, iDet, arr_vD_fit[iDet], iDet, arr_LA_fit[iDet] * TMath::RadToDeg());

    delete min;
#endif

    //i_point = i_point+1;
    //----------------------------------------------------------------------

  } // end of detector loop
    //----------------------------------------------------------------------

#if 1

  TH1D* th_vD_out = new TH1D("vD", "vD", 540, 1, 540);
  TH1D* th_LA_out = new TH1D("LA", "LA", 540, 1, 540);

  for (int iDet = 0; iDet < o2::trd::constants::MAXCHAMBER; ++iDet) {
    th_vD_out->Fill(iDet + 1, arr_vD_fit[iDet]);
    th_LA_out->Fill(iDet + 1, arr_LA_fit[iDet]);
  }

#endif

  auto fOut = TFile::Open("trdcalibdummy.root", "recreate");
  fOut->cd();

  for (int iDet = 0; iDet < o2::trd::constants::MAXCHAMBER; ++iDet) {
    arr_tp_angleDiff[iDet]->Write();
    //arr_tp_angleDiff[iDet].reset();
  }
  th_vD_out->Write();
  th_LA_out->Write();

  fOut->Close();
}

Slot& CalibratorVdExB::emplaceNewSlot(bool front, uint64_t tStart, uint64_t tEnd)
{
  auto& container = getSlots();
  auto& slot = front ? container.emplace_front(tStart, tEnd) : container.emplace_back(tStart, tEnd);
  slot.setContainer(std::make_unique<AngularResidHistos>());
  return slot;
}

} // namespace o2::trd
