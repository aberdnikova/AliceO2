// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file   CalibVDrift.cxx
/// \author Ole Schmidt, ole.schmidt@cern.ch

#include "TFile.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TMath.h"
#include "TVirtualFitter.h"
#include "Math/MinimizerOptions.h"
#include "TMinuitMinimizer.h"
#include "CCDB/CcdbApi.h"
//#include "AliTRDCalPad.h"
#include "TRDBase/LocalVDrift.h"
#include "TRDBase/ChamberCalibrations.h"
#include "CCDB/BasicCCDBManager.h"



#include <fairlogger/Logger.h>

#include "TRDCalibration/CalibVDrift.h"

using namespace std;
using namespace o2::ccdb;
using namespace o2::trd;
using namespace o2::trd::constants;

static std::array<TProfile*, constants::MAXCHAMBER> vec_tp_Delta_vs_impact_circle; //for input hists
static int iDet;

//use fixed pre-calib parameters
static Double_t vD_pre_corr_set = 1.546;  // 1.546, 1.2
static Double_t LA_pre_corr_set = -0.16133;

static Double_t low_angle_fit = 80.0;
static Double_t up_angle_fit  = 100.0;
\

#if 1
TGraph* calc_Delta_alpha(Double_t vD_fit, Double_t Lorentz_angle_fit, Double_t vD_pre_corr, Double_t Lorentz_angle_pre_corr)
{
    // 06.10.2020 A. Schmah Added comments, replaced drift_vel_ratio by (vD_pre_corr/vD_fit)

    // The first cluster (the one at the bottom, last time time) is drifting to the anode plane. The position there is
    // independent of the true drift velocity and depends only on the Lorentz angle. BUT the time bin depends on the
    // true drift velocity. For large vDtrue the time bin is small, e.g. 18, for small vDtrue the time bin is large, e.g. 24.
    // The cluster is then projected with the fixed drift velocity, e.g. vD_pre_corr = 1.546, along local y. Therefore Delta alpha vs. impact angle
    // depends on vDtrue and vD_pre_corr and the Lorentz angle.

    TGraph* TG_Delta_alpha_vs_impact_angle = new TGraph();

    Int_t i_point = 0;
    for(Double_t impact_angle = 65.0*TMath::DegToRad(); impact_angle < 115.0*TMath::DegToRad(); impact_angle += 1.0*TMath::DegToRad())
    {

        Double_t TRD_anode_plane = 0.0335;

        // Direction vector of incoming track
        Double_t x_dir = TMath::Cos(impact_angle);
        Double_t y_dir = TMath::Sin(impact_angle);

        // Slope of incoming track
        Double_t slope = 10000000.0;
        if(x_dir != 0.0) slope = y_dir/x_dir;

        // Slope of Lorentz angle
        Double_t Lorentz_tan   = TMath::Tan(Lorentz_angle_fit);
        Double_t Lorentz_slope = 10000000.0;
        if(Lorentz_tan != 0.0) Lorentz_slope = 1.0/Lorentz_tan;

        // Hit point of incoming track with anode plane
        Double_t x_anode_hit = TRD_anode_plane/slope;
        Double_t y_anode_hit = TRD_anode_plane;

        // Hit point at anode plane of Lorentz angle shifted cluster from the entrance -> independent of vDtrue
        Double_t x_Lorentz_anode_hit = TRD_anode_plane/Lorentz_slope;
        Double_t y_Lorentz_anode_hit = TRD_anode_plane;

        // Cluster location within drift cell of cluster from entrance after drift velocity ratio is applied
        Double_t x_Lorentz_drift_hit = x_Lorentz_anode_hit;
        Double_t y_Lorentz_drift_hit = TRD_anode_plane - TRD_anode_plane*(vD_pre_corr/vD_fit);

        // Reconstructed hit of first cluster at chamber entrance after pre Lorentz angle correction
        Double_t x_Lorentz_drift_hit_pre_corr = x_Lorentz_anode_hit - (TRD_anode_plane - y_Lorentz_drift_hit)*TMath::Tan(Lorentz_angle_pre_corr);
        Double_t y_Lorentz_drift_hit_pre_corr = TRD_anode_plane - TRD_anode_plane*(vD_pre_corr/vD_fit);

        Double_t impact_angle_track = TMath::ATan2(y_anode_hit,x_anode_hit);

        Double_t Delta_x_Lorentz_drift_hit = x_anode_hit - x_Lorentz_drift_hit_pre_corr;
        Double_t Delta_y_Lorentz_drift_hit = y_anode_hit - y_Lorentz_drift_hit_pre_corr;
        Double_t impact_angle_rec   = TMath::ATan2(Delta_y_Lorentz_drift_hit,Delta_x_Lorentz_drift_hit);

        Double_t Delta_angle = -(impact_angle_track - impact_angle_rec);
        TG_Delta_alpha_vs_impact_angle ->SetPoint(i_point,impact_angle*TMath::RadToDeg(),Delta_angle*TMath::RadToDeg());
        i_point++;
        //printf("impact_angle: %4.3f, impact_angle_track: %4.3f, impact_angle_rec: %4.3f, Delta: {%4.3f, %4.3f}, x_anode_hit: %4.3f, x_Lorentz_drift_hit: %4.3f \n",impact_angle*TMath::RadToDeg(),impact_angle_track*TMath::RadToDeg(),impact_angle_rec*TMath::RadToDeg(),Delta_x_Lorentz_drift_hit,Delta_y_Lorentz_drift_hit,x_anode_hit,x_Lorentz_drift_hit);
    }

    return TG_Delta_alpha_vs_impact_angle;
}
//----------------------------------------------------------------------------------------
#endif

#if 1
// function to be minimized
void Chi2_TRD_vDrift(Int_t &, Double_t *, Double_t & sum, Double_t * par, Int_t )
{
    sum = 0;

    //Double_t B_field_use = par[0];
    //Double_t E_field_use = par[1];
    //Double_t v_drift_use = par[2];
    //Double_t vD_use      = par[3];
    //Double_t LA_use      = par[4];

    Double_t LA_fit  = par[0];
    Double_t vD_fit  = par[1];
    

    //Double_t vD_fit, Double_t Lorentz_angle_fit, Double_t vD_pre_corr, Double_t Lorentz_angle_pre_corr
    TGraph* tg_Delta_vs_impact_single = calc_Delta_alpha(vD_fit,LA_fit,vD_pre_corr_set,LA_pre_corr_set); // new version
    //TGraph* tg_Delta_vs_impact_single = calc_Delta_alpha(LA_use,vD_ratio_use); // old version
    //TGraph* tg_Delta_vs_impact_single = calc_Delta_alpha(LA_use,vD_ratio_use,Lorentz_angle_pre_corr); // old version

    for(Int_t i_bin = 1; i_bin <= vec_tp_Delta_vs_impact_circle[iDet]->GetNbinsX(); i_bin++)
    {
        Double_t impact_angle     = vec_tp_Delta_vs_impact_circle[iDet] ->GetBinCenter(i_bin);
        Double_t Delta_alpha      = vec_tp_Delta_vs_impact_circle[iDet] ->GetBinContent(i_bin);

        if(Delta_alpha == 0.0) continue;
        if(impact_angle < low_angle_fit) continue;
        if(impact_angle > up_angle_fit)  continue;

        Double_t Delta_alpha_sim = tg_Delta_vs_impact_single ->Eval(impact_angle);

        //sum += TMath::Power(Delta_alpha_sim - Delta_alpha,2)/TMath::Power(Delta_alpha_err,2);
        sum += TMath::Power(Delta_alpha_sim - Delta_alpha,2)/1.0;

    }

    delete tg_Delta_vs_impact_single;
}
//----------------------------------------------------------------------------------------
#endif

void CalibVDrift::process()
{
    LOG(info) << "Started processing for vDrift calibration";


#if 1
    //------------------------------------------
    // Connect to the CCDB
    // Read
    printf("Connect to CCDB for reading \n");
    Long64_t timestamp = 265501;
    LocalVDrift* o2localvdrift_read;
    ChamberCalibrations* o2chambercalibrations_read;
    //auto& ccdbmgr = o2::ccdb::BasicCCDBManager::instance();
    BasicCCDBManager ccdbmgr = o2::ccdb::BasicCCDBManager::instance();
    ccdbmgr.setTimestamp(timestamp);
    //std::string fulltablename = "TRD_test/LocalVDrift/" + tablename;
    std::string fulltablename       = "TRD_test/LocalVDrift";
    std::string fulltablename_calib = "TRD_test/ChamberCalibrations";
    o2localvdrift_read         = ccdbmgr.get<o2::trd::LocalVDrift>(fulltablename);
    o2chambercalibrations_read = ccdbmgr.get<o2::trd::ChamberCalibrations>(fulltablename_calib);
    Double_t value = o2localvdrift_read ->getPadValue(0,0,0); // PadCalibrations.h
    printf("value: %4.3f \n",value);

    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        //Float_t gain = o2chambercalibrations_read->getGainFactor(i_det);
        Float_t vdrift = o2chambercalibrations_read->getVDrift(i_det);
        printf("i_det: %d, vdrift: %4.3f \n",i_det,vdrift);
    }

    //Connect to CCDB
    //

    /*
    printf("Connect to CCDB for writing \n");
    o2::ccdb::CcdbApi ccdb;
    map<string, string> metadata;               // do we want to store any meta data?
    ccdb.init("http://ccdb-test.cern.ch:8080"); // or http://localhost:8080 for a local installation

    //auto o2localvdrift = new o2::trd::LocalVDrift();
    LocalVDrift* o2localvdrift_write = new o2::trd::LocalVDrift();
    ChamberCalibrations* o2chambercalibrations = new o2::trd::ChamberCalibrations();
    Int_t j_ROC = 0;
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        o2localvdrift_write->setPadValue(i_det, j_ROC, 1.56);
        o2chambercalibrations->setVDrift(i_det, 1.82);
    }

    // Write
    printf("Write to CCDB \n");
    Int_t Run = timestamp;
    ccdb.storeAsTFileAny(o2localvdrift_write, "TRD_test/LocalVDrift", metadata, Run, Run + 1);
    ccdb.storeAsTFileAny(o2chambercalibrations, "TRD_test/ChamberCalibrations", metadata, Run, Run + 1);
    */
    //------------------------------------------
#endif


    //std::array<TProfile*, constants::MAXCHAMBER> vec_tp_Delta_vs_impact_circle; //for input hists

    std::array<double, constants::MAXCHAMBER> arr_LA_fit;
    std::array<double, constants::MAXCHAMBER> arr_vD_fit;

    auto input_data = TFile::Open("/home/ceres/berdnikova/TRD-Run3-Calibration/Data/TRD_Calib_on_trkl.root");

    for (int iDet = 0; iDet < constants::MAXCHAMBER; ++iDet)
    {
        vec_tp_Delta_vs_impact_circle[iDet] = (TProfile*)input_data->Get(Form("Delta_impact_circle/vec_th1d_Delta_vs_impact_circle_%d",iDet));
    }

    std::array<std::unique_ptr<TGraphErrors>, constants::MAXCHAMBER> arr_tg_angleDiff;

    for (iDet = 0; iDet < constants::MAXCHAMBER; ++iDet) {

        arr_tg_angleDiff[iDet] = std::make_unique<TGraphErrors>(constants::NBINSANGLEDIFF); // xy distribution of TRD space points

        //add array with RMSs - on the fly
        //second array with counters?

        //tracklets are calibrated or not?

        for (int iBin = 0; iBin < constants::NBINSANGLEDIFF; ++iBin) { // note: iBin = constants::NBINSANGLEDIFF - 1 is under-/overflow bin

        //for (int iBin = 5; iBin < 20; ++iBin) {
        /* //read real values - enable later
            short nEntries = mAngleDiffCounters[iDet * constants::NBINSANGLEDIFF + iBin];
            float angleDiffSum = mAngleDiffSums[iDet * constants::NBINSANGLEDIFF + iBin];
            if (nEntries > 0)
            {
                LOGF(INFO, "Found %i entrie(s) in chamber %i, bin %i. Average angular deviation: %f", nEntries, iDet, iBin, angleDiffSum / nEntries);
                arr_tg_angleDiff[iDet] ->SetPoint(iBin,iBin,angleDiffSum/nEntries);
                LOGF(INFO, "Point written: bin %i, angdiff: %f",iBin, arr_tg_angleDiff[iDet] ->GetPointY(iBin));
            }
            */
            arr_tg_angleDiff[iDet] ->SetPoint(iBin,vec_tp_Delta_vs_impact_circle[iDet]->GetBinCenter(iBin+212),vec_tp_Delta_vs_impact_circle[iDet]->GetBinContent(iBin+212));
            //LOGF(INFO, "Point written: bin %i, x: %f, angdiff: %f",iBin,arr_tg_angleDiff[iDet] ->GetPointX(iBin),arr_tg_angleDiff[iDet] ->GetPointY(iBin));

            //minimization

            TVirtualFitter *min = TVirtualFitter::Fitter(0,2);
            min->SetFCN(Chi2_TRD_vDrift);
            Double_t pStart[2] = {-7.5*TMath::DegToRad(),1.0};

            min->SetParameter(0,"LA",pStart[0],0.01,0,0);
            min->SetParameter(1,"vD",pStart[1],0.01,0,0);

            Double_t arglist[2];
            arglist[0] = 1000; // number of function calls
            arglist[1] = 0.001; // tolerance
            //arglist[0] = -1.0;
            //min->ExecuteCommand("SET PRINT", arglist, 1);

            Double_t arglist_A[1] = {-1};
            Double_t arglist_B[1] = {0};
            Double_t arglist_C[1] = {-1};

            min->ExecuteCommand("SET PRIntout",arglist_A,1);
            min->ExecuteCommand("SET NOWarnings",arglist_B,1);

            min->ExecuteCommand("MIGRAD",arglist,0);


            Double_t parFit[5];

            for(int i = 0; i < 2; ++i)
            {
                parFit[i] = min->GetParameter(i);
            }

            arr_LA_fit[iDet] = parFit[0];
            arr_vD_fit[iDet] = parFit[1];

            LOGF(INFO, "iDet: %i, arr_vD_fit[%i]: %4.3f, arr_LA_fit[%i]: %4.3f",iDet,iDet,arr_vD_fit[iDet],iDet,arr_LA_fit[iDet]);

            delete min;

            //i_point = i_point+1;


        }
    }

    auto fOut = TFile::Open("trdcalibdummy.root", "recreate");
    fOut->cd();

    for (int iDet = 0; iDet < constants::MAXCHAMBER; ++iDet)
    {
        arr_tg_angleDiff[iDet] ->Write();
        arr_tg_angleDiff[iDet].reset();
    }

    fOut->Close();

  /*
  // as an example I loop over the input, create a histogram and write it to a file
  auto fOut = TFile::Open("trdcalibdummy.root", "recreate");
  auto hXY = std::make_unique<TH2F>("histDummy", "foo", 100, -60, 60, 100, 250, 400); // xy distribution of TRD space points
  for (int i = 0; i < mAngulerDeviationProf.size(); ++i) {
    hXY->Fill(mAngulerDeviationProf[i].mX[0], mAngulerDeviationProf[i].mR);
  }
  fOut->cd();
  hXY->Write();
  hXY.reset(); // delete the histogram before closing the output file
  fOut->Close();
  */
}

