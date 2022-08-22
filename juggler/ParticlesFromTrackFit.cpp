// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Sylvester Joosten, Wouter Deconinck

#include <algorithm>

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/RndmGenerators.h"
#include "Gaudi/Property.h"

#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/Surface.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"

// Event Model related classes
#include "eicd/ReconstructedParticleCollection.h"
#include "eicd/TrackerHitCollection.h"
#include "eicd/TrackParametersCollection.h"
#include "JugTrack/IndexSourceLink.hpp"
#include "JugTrack/Track.hpp"
#include "JugTrack/Trajectories.hpp"

#include "Acts/Utilities/Helpers.hpp"

#include "eicd/vector_utils.h"

#include <cmath>

namespace Jug::Reco {

  /** Extract the particles form fit trajectories.
   *
   * \ingroup tracking
   */
   class ParticlesFromTrackFit : public GaudiAlgorithm {
   private:
    DataHandle<TrajectoriesContainer>     m_inputTrajectories{"inputTrajectories", Gaudi::DataHandle::Reader, this};
    DataHandle<eicd::ReconstructedParticleCollection> m_outputParticles{"outputParticles", Gaudi::DataHandle::Writer, this};
    DataHandle<eicd::TrackParametersCollection> m_outputTrackParameters{"outputTrackParameters", Gaudi::DataHandle::Writer, this};

   public:
    //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    ParticlesFromTrackFit(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc) {
          declareProperty("inputTrajectories", m_inputTrajectories,"");
          declareProperty("outputParticles", m_outputParticles, "");
          declareProperty("outputTrackParameters", m_outputTrackParameters, "Acts Track Parameters");
        }

    StatusCode initialize() override {
      if (GaudiAlgorithm::initialize().isFailure()) {
        return StatusCode::FAILURE;
      }
      return StatusCode::SUCCESS;
    }

    StatusCode execute() override {
      // input collection
      const auto* const trajectories = m_inputTrajectories.get();
      // create output collections
      auto* rec_parts = m_outputParticles.createAndPut();
      auto* track_pars = m_outputTrackParameters.createAndPut();

      if (msgLevel(MSG::DEBUG)) {
        debug() << std::size(*trajectories) << " trajectories " << endmsg;
      }

      // Loop over the trajectories
        for (const auto& traj : *trajectories) {

          // Get the entry index for the single trajectory
          // The trajectory entry indices and the multiTrajectory
          const auto& mj        = traj.multiTrajectory();
          const auto& trackTips = traj.tips();
          if (trackTips.empty()) {
            if (msgLevel(MSG::DEBUG)) {
              debug() << "Empty multiTrajectory." << endmsg;
            }
            continue;
          }

          const auto& trackTip = trackTips.front();

          // Collect the trajectory summary info
          auto trajState       = Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);
          int  m_nMeasurements = trajState.nMeasurements;
          int  m_nStates       = trajState.nStates;
          int  m_nCalibrated   = 0;
          int  m_hit_layers[12] = {0}; // reserve 3 vertexing + 7 tracking + 1 pid + 1 ecal surfaces
          float m_chi2_layers[12] = {0};
          int  m_hit_binary = 0;
          if (msgLevel(MSG::DEBUG)) {
            debug() << "n measurement in trajectory " << m_nMeasurements << endmsg;
            debug() << "n state in trajectory " << m_nStates << endmsg;
          }

          mj.visitBackwards(trackTip, [&](auto&& trackstate) {
            // debug() << trackstate.hasPredicted() << endmsg;
            // debug() << trackstate.predicted() << endmsg;
            auto params = trackstate.predicted(); //<< endmsg;
            auto pathlength = trackstate.pathLength();
            auto geoID = trackstate.referenceSurface().geometryId();
            auto volume = geoID.volume();
            auto layer = geoID.layer();
            if (trackstate.hasCalibrated())
            {
              m_nCalibrated++;

              //=========================================================
              //              Vertexing (3 Si layers)
              //=========================================================
              if ((volume==20 || volume==18) && layer==2 && trackstate.pathLength()<130) {m_hit_layers[0] = 1; m_chi2_layers[0] = trackstate.chi2();} // 1st Si vertex (3.3 cm)
              if ((volume==20 || volume==18) && layer==4 && trackstate.pathLength()<130) {m_hit_layers[1] = 1; m_chi2_layers[1] = trackstate.chi2();} // 2nd Si vertex (4.35 cm)
              if ((volume==20 || volume==18) && layer==6 && trackstate.pathLength()<130) {m_hit_layers[2] = 1; m_chi2_layers[2] = trackstate.chi2();} // 3rd Si vertex (5.4 cm)

              //=========================================================
              //       Barrel Tracking (2 Si layers + 4 MM layers)
              //=========================================================
              if ((volume==22 || volume==20) && layer==2 && trackstate.pathLength()>130) {m_hit_layers[3] = 1; m_chi2_layers[3] = trackstate.chi2();} // 1st Si layer (13.34 cm)
              if ((volume==22 || volume==20) && layer==4 && trackstate.pathLength()>130) {m_hit_layers[4] = 1; m_chi2_layers[4] = trackstate.chi2();} // 2nd Si layer (17.96 cm)
              if ((volume==25 || volume==23) && layer==2) {m_hit_layers[5] = 1; m_chi2_layers[5] = trackstate.chi2();} // 1st MPGD layer (47.72 cm)
              if ((volume==25 || volume==23) && layer==4) {m_hit_layers[6] = 1; m_chi2_layers[6] = trackstate.chi2();} // 2nd MPGD layer (49.57 cm)
              if ((volume==28 || volume==26) && layer==2) {m_hit_layers[7] = 1; m_chi2_layers[7] = trackstate.chi2();} // 3rd MPGD layer (75.61 cm)
              if ((volume==28 || volume==26) && layer==4) {m_hit_layers[8] = 1; m_chi2_layers[8] = trackstate.chi2();} // 4th MPGD layer (77.46 cm)

              //=========================================================
              //       Forward Tracking (6 Si disks + 3 GEM disks)
              //=========================================================
              if ((volume==23 || volume==21) && layer==2) {m_hit_layers[3] = 1; m_chi2_layers[3] = trackstate.chi2();} // 1st Si disk (25.0 cm)
              if ((volume==23 || volume==21) && layer==4) {m_hit_layers[4] = 1; m_chi2_layers[4] = trackstate.chi2();} // 2th Si disk (49.0 cm)
              if ((volume==26 || volume==24) && layer==2) {m_hit_layers[5] = 1; m_chi2_layers[5] = trackstate.chi2();} // 3rd Si disk (73.0 cm)
              if ((volume==29 || volume==27) && layer==2) {m_hit_layers[6] = 1; m_chi2_layers[6] = trackstate.chi2();} // 4th Si disk (103.65 cm)
              if (volume==29 && layer==4) {m_hit_layers[6] = 1; m_chi2_layers[6] = trackstate.chi2();} // 1st GEM disk (105.76 cm)
              if ((volume==29 && layer==6) || (volume==27 && layer==4)) {m_hit_layers[7] = 1; m_chi2_layers[7] = trackstate.chi2();} // 5th Si disk (134.22 cm)
              if (volume==29 && layer==8) {m_hit_layers[8] = 1; m_chi2_layers[8] = trackstate.chi2();} // 2nd GEM disk (161.74 cm) ** cannot find it
              if ((volume==29 && layer==10) || (volume==27 && layer==6) ) {m_hit_layers[8] = 1; m_chi2_layers[8] = trackstate.chi2();} // 6th Si disk (165.0 cm)
              if (volume==31 && layer==2) {m_hit_layers[9] = 1; m_chi2_layers[9] = trackstate.chi2();} // 3rd GEM disk (332.0 cm)

              //=========================================================
              //       Backward Tracking (5 Si disks + 3 GEM disks)
              //=========================================================
              if ((volume==16 || volume==14) && layer==4) {m_hit_layers[3] = 1; m_chi2_layers[3] = trackstate.chi2();} // 1st Si disk (-25.0 cm)
              if ((volume==16 || volume==14) && layer==2) {m_hit_layers[4] = 1; m_chi2_layers[4] = trackstate.chi2();} // 2th Si disk (-49.0 cm)
              if ((volume==11 || volume==9) && layer==2) {m_hit_layers[5] = 1; m_chi2_layers[5] = trackstate.chi2();} // 3rd Si disk (-73.0 cm)
              if (volume==6 && layer==8) {m_hit_layers[6] = 1; m_chi2_layers[6] = trackstate.chi2();} // 1st GEM disk (-103.0 cm)
              if ((volume==6 || volume==4) && layer==6) {m_hit_layers[6] = 1; m_chi2_layers[6] = trackstate.chi2();} // 4th Si disk (-109.0 cm)
              if (volume==6 && layer==4) {m_hit_layers[7] = 1; m_chi2_layers[7] = trackstate.chi2();} // 2nd GEM disk (-141.74 cm)
              if ((volume==6 || volume==4) && layer==2) {m_hit_layers[7] = 1; m_chi2_layers[7] = trackstate.chi2();} // 4th Si disk (-145.0 cm)
            } 

            if (msgLevel(MSG::DEBUG)) {
              debug() << "******************************" << endmsg;
              debug() << "predicted variables: \n" << trackstate.predicted() << endmsg;
              debug() << "geoID = " << geoID << endmsg;
              debug() << "volume = " << volume << ", layer = " << layer << endmsg;
              debug() << "pathlength = " << pathlength << endmsg;
              debug() << "hasCalibrated = " << trackstate.hasCalibrated() << endmsg;
              debug() << "track state type = " << trackstate.typeFlags().to_string() << endmsg;

              debug() << "pos 0 = " << params[Acts::eFreePos0] << endmsg;
              debug() << "pos 1 = " << params[Acts::eFreePos1] << endmsg;
              debug() << "pos 2 = " << params[Acts::eFreePos2] << endmsg;
              debug() << "******************************" << endmsg;
            }
          });

          for (int ilayer = 0; ilayer < 10; ++ilayer) {
            if (m_hit_layers[ilayer]) m_hit_binary += pow(2,ilayer);
            if (msgLevel(MSG::DEBUG) && m_hit_layers[ilayer]) debug() << "Hit @ " << ilayer << " with chi2 " << m_chi2_layers[ilayer] << endmsg; 
          }

          if (msgLevel(MSG::DEBUG)) {
            debug() << "n calibrated state in trajectory " << m_nCalibrated << " binary number " << std::bitset<12>(m_hit_binary).to_string() << endmsg;
          }


          // Get the fitted track parameter
          //
          float m_chi2Sum = -9999;
          if (traj.hasTrackParameters(trackTip)) {
            const auto& boundParam = traj.trackParameters(trackTip);
            const auto& parameter  = boundParam.parameters();
            const auto& covariance = *boundParam.covariance();
            m_chi2Sum = trajState.chi2Sum;

            if (msgLevel(MSG::DEBUG)) {
              debug() << "loc 0 = " << parameter[Acts::eBoundLoc0] << endmsg;
              debug() << "loc 1 = " << parameter[Acts::eBoundLoc1] << endmsg;
              debug() << "phi   = " << parameter[Acts::eBoundPhi] << endmsg;
              debug() << "theta = " << parameter[Acts::eBoundTheta] << endmsg;
              debug() << "q/p   = " << parameter[Acts::eBoundQOverP] << endmsg;
              debug() << "p     = " << 1.0 / parameter[Acts::eBoundQOverP] << endmsg;

              debug() << "err phi = " << sqrt(covariance(Acts::eBoundPhi, Acts::eBoundPhi)) << endmsg;
              debug() << "err th  = " << sqrt(covariance(Acts::eBoundTheta, Acts::eBoundTheta)) << endmsg;
              debug() << "err q/p = " << sqrt(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP)) << endmsg;

              debug() << " chi2 = " << trajState.chi2Sum << endmsg;
            }

            const decltype(eicd::TrackParametersData::loc) loc {
              parameter[Acts::eBoundLoc0],
              parameter[Acts::eBoundLoc1]
            };
            const decltype(eicd::TrackParametersData::momentumError) momentumError {
              static_cast<float>(covariance(Acts::eBoundTheta, Acts::eBoundTheta)),
              static_cast<float>(covariance(Acts::eBoundPhi, Acts::eBoundPhi)),
              static_cast<float>(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP)),
              static_cast<float>(covariance(Acts::eBoundTheta, Acts::eBoundPhi)),
              static_cast<float>(covariance(Acts::eBoundTheta, Acts::eBoundQOverP)),
              static_cast<float>(covariance(Acts::eBoundPhi, Acts::eBoundQOverP))};
            const decltype(eicd::TrackParametersData::locError) locError {
              static_cast<float>(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0)),
              static_cast<float>(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1)),
              static_cast<float>(covariance(Acts::eBoundLoc0, Acts::eBoundLoc1))};
            // const float timeError{sqrt(static_cast<float>(covariance(Acts::eBoundTime, Acts::eBoundTime)))};

            eicd::TrackParameters pars{
              0, // type: track head --> 0
              loc,
              locError,
              static_cast<float>(parameter[Acts::eBoundTheta]),
              static_cast<float>(parameter[Acts::eBoundPhi]),
              static_cast<float>(parameter[Acts::eBoundQOverP]),
              momentumError,
              static_cast<float>(parameter[Acts::eBoundTime]),
              m_chi2Sum, // timeError
              static_cast<float>(m_hit_binary)}; // boundParam.charge()
            track_pars->push_back(pars);
          }

          auto tsize = trackTips.size();
          if (msgLevel(MSG::DEBUG)) {
            debug() << "# fitted parameters : " << tsize << endmsg;
          }
          if (tsize == 0) {
            continue;
          }

          mj.visitBackwards(tsize - 1, [&](auto&& trackstate) {
            // debug() << trackstate.hasPredicted() << endmsg;
            // debug() << trackstate.predicted() << endmsg;
            auto params = trackstate.predicted(); //<< endmsg;

            double p0 = (1.0 / params[Acts::eBoundQOverP]) / Acts::UnitConstants::GeV;
            if (msgLevel(MSG::DEBUG)) {
              debug() << "track predicted p = " << p0 << " GeV" << endmsg;
            }
            if (std::abs(p0) > 500) {
              if (msgLevel(MSG::DEBUG)) {
                debug() << "skipping" << endmsg;
              }
              return;
            }

            auto rec_part = rec_parts->create();
            rec_part.setMomentum(
              eicd::sphericalToVector(
                1.0 / std::abs(params[Acts::eBoundQOverP]),
                params[Acts::eBoundTheta],
                params[Acts::eBoundPhi])
            );
            rec_part.setCharge(static_cast<int16_t>(std::copysign(1., params[Acts::eBoundQOverP])));
          });
      }

      return StatusCode::SUCCESS;
    }

  };
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  DECLARE_COMPONENT(ParticlesFromTrackFit)

} // namespace Jug::Reco
