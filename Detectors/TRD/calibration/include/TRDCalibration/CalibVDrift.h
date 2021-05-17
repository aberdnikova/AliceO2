// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_TRD_CALIBVDRIFT_H_
#define ALICEO2_TRD_CALIBVDRIFT_H_

/// \file   CalibVDrift.h
/// \author Ole Schmidt, ole.schmidt@cern.ch

#include "DataFormatsTRD/Constants.h"
#include "TProfile.h"
#include "TGraph.h"
#include <array>

namespace o2
{
namespace trd
{

/// \brief VDrift calibration class
///
/// This class is used to determine chamber-wise vDrift values
///
/// origin: TRD
    /// \author Ole Schmidt, ole.schmidt@cern.ch

    static std::array<std::unique_ptr<TProfile>, constants::MAXCHAMBER> arr_tp_angleDiff;
    //static std::array<TProfile*, constants::MAXCHAMBER> arr_tp_angleDiff;

static std::array<double, constants::MAXCHAMBER> arr_LA_fit;
static std::array<double, constants::MAXCHAMBER> arr_vD_fit;

  //bool mIsInitialized{false};

class CalibVDrift
{
 public:
  /// default constructor
  CalibVDrift() = default;

  /// default destructor
  ~CalibVDrift() = default;

  /// set input angular difference sums
  void setAngleDiffSums(float* input)
  {
    for (int i = 0; i < constants::MAXCHAMBER * constants::NBINSANGLEDIFF; ++i) {
      mAngleDiffSums[i] = input[i];
    }
  }

  /// set input angular difference bin counters
  void setAngleDiffCounters(short* input)
  {
    for (int i = 0; i < constants::MAXCHAMBER * constants::NBINSANGLEDIFF; ++i) {
      mAngleDiffCounters[i] = input[i];
    }
  }

  void init();

  /// main processing function
  void process();

  //TGraph* calc_Delta_alpha(Double_t vD_fit, Double_t Lorentz_angle_fit, Double_t vD_pre_corr, Double_t Lorentz_angle_pre_corr);
  //void Chi2_TRD_vDrift(Int_t &, Double_t *, Double_t & sum, Double_t * par, Int_t );

 private:
  std::array<float, constants::MAXCHAMBER * constants::NBINSANGLEDIFF> mAngleDiffSums{};     ///< input TRD track to tracklet angular difference sums per bin
  std::array<short, constants::MAXCHAMBER * constants::NBINSANGLEDIFF> mAngleDiffCounters{}; ///< input bin counters

 //public:

  //static std::array<std::unique_ptr<TProfile>, constants::MAXCHAMBER> arr_tp_angleDiff;

  //static std::array<double, constants::MAXCHAMBER> arr_LA_fit;
  //static std::array<double, constants::MAXCHAMBER> arr_vD_fit;

  bool mIsInitialized{false};

  //short nEntries;
  //float angleDiffSum;

  //define tgraphs
};

} // namespace trd

} // namespace o2
#endif
