#include <TROOT.h>
#include <TApplication.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <Track.h>
#include <TrackPoint.h>
#include <PlanarMeasurement.h>
#include <DetPlane.h>
#include <RKTrackRep.h>
#include <KalmanFitterRefTrack.h>
#include <FieldManager.h>
#include <ConstField.h>
#include <MaterialEffects.h>
#include <TGeoMaterialInterface.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TMatrixDSym.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <memory>
#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>

// === Global constants ===
constexpr double SMEAR = 0.5;
constexpr int N_POINTS = 15;
constexpr int N_RUNS = 500;

// Function that runs the fit once and appends results to storage containers
void run_one_fit(
    TRandom3& rnd,
    std::vector<double>& residualsU,
    std::vector<double>& residualsV
) {
  // Create random reference direction and points
    TLorentzVector refVec;
    refVec.SetXYZT(
        rnd.Uniform(5, 10),
        rnd.Uniform(5, 10),
        rnd.Uniform(5, 10),
        rnd.Uniform(0.01, 0.2)
    );

    const TVector3 refDir = refVec.Vect().Unit();
    const TVector3 refStart(0, 0, 0);

    // Generate noisy measurements along the track
    std::vector<TVector3> refPoints;
    std::vector<TVector3> measPoints;
    double step = 0.0;
    for (int i = 0; i < N_POINTS; ++i) {
        step += rnd.Uniform(5.0, 15.0);
        TVector3 pos = refStart + step * refDir;
        refPoints.push_back(pos); // reference position for detector plane
        TVector3 meas(rnd.Gaus(pos.X(), SMEAR), rnd.Gaus(pos.Y(), SMEAR), rnd.Gaus(pos.Z(), SMEAR));
        measPoints.push_back(meas); // actual measurement, smeared
    }

    // Create the track rep and the track (use default constructor)
    auto* rep = new genfit::RKTrackRep(13); // muon
    genfit::Track track;
    track.addTrackRep(rep);

    // Fill track with TrackPoints/PlanarMeasurements, each with a valid DetPlane
    for (size_t i = 0; i < measPoints.size(); ++i) {
        const TVector3 planeO = refPoints.at(i);
        const TVector3 planeN = refDir;

        const TVector3 pointOnPlane = measPoints.at(i) - planeO;

        // Local axes for the plane
        TVector3 U = planeN.Orthogonal().Unit();
        TVector3 V = planeN.Cross(U).Unit();

        TVectorD uv(2);
        uv[0] = pointOnPlane.Dot(U);
        uv[1] = pointOnPlane.Dot(V);

        TMatrixDSym cov(2);
        cov.UnitMatrix();
        cov *= (SMEAR*SMEAR);

        // Create PlanarMeasurement
        auto* meas = new genfit::PlanarMeasurement(uv, cov, 0, i, nullptr);
        meas->setPlane(std::make_shared<genfit::DetPlane>(planeO, planeN), 0);

        track.insertPoint(new genfit::TrackPoint(meas, &track));
    }

    // Set initial seed state (position, momentum)
    TVectorD seedState(6);
    seedState[0] = refStart.X();
    seedState[1] = refStart.Y();
    seedState[2] = refStart.Z();
    seedState[3] = refDir.X();
    seedState[4] = refDir.Y();
    seedState[5] = refDir.Z();
    TMatrixDSym seedCov(6);
    seedCov.UnitMatrix();
    seedCov *= (SMEAR*SMEAR);

    track.setStateSeed(seedState);
    track.setCovSeed(seedCov);

    // Fit
    genfit::KalmanFitterRefTrack fitter;
    try {
        fitter.processTrack(&track);
    } catch (const genfit::Exception& e) {
        return;
    }

    // Extrapolate to a plane at the reference point/direction
    auto* fittedRep = track.getCardinalRep();
    auto plane = std::make_shared<genfit::DetPlane>(refStart, refDir); // target plane

    // first point in track; use the 5D state vector obtained from the track and
    // extrapolate the fitted state, not a constructed state from the seed
    auto state = track.getFittedState(0);

    try {
        fittedRep->extrapolateToPlane(state, plane, true, false);
    } catch (const genfit::Exception& e) {
        return;
    }
    TVector3 fitPos = state.getPos();

    // residual: position on readout plane (distance to origin of plane), which
    // will be projected to U/V; histogramming the absolute value pos-O would be
    // a Raleight distribution (no negative values possible)
    const double resU = (fitPos - plane->getO()).Dot(plane->getU());
    const double resV = (fitPos - plane->getO()).Dot(plane->getV());

    residualsU.push_back(resU);
    residualsV.push_back(resV);
}

double calculateMean(const std::vector<double>& values) {
    if (values.empty()) {
        return 0.0;
    }
    double sum = std::accumulate(values.begin(), values.end(), 0.0);
    return sum / values.size();
}

double calculateStandardDeviation(const std::vector<double>& values, double mean) {
    if (values.empty()) {
        return 0.0;
    }
    double squaredSum = 0.0;
    for (double value : values) {
        squaredSum += (value - mean) * (value - mean);
    }
    double variance = squaredSum / values.size();
    return std::sqrt(variance);
}

int main() {
    // Initialize GenFit global services only once
    genfit::FieldManager::getInstance()->init(new genfit::ConstField(0., 0., 0.));
    genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
    genfit::MaterialEffects::getInstance()->setNoEffects();

    TRandom3 rnd(0); // single RNG instance

    std::vector<double> residualsU;
    std::vector<double> residualsV;

    for (int i = 0; i < N_RUNS; ++i) {
        run_one_fit(rnd, residualsU, residualsV);
    }

    // Needed for interactive ROOT
    //TApplication app("app", 0, nullptr);
    gROOT->SetBatch(kTRUE); // no GUI

    // Create histograms
    TH1D* hResU = new TH1D("hResU", "Residuals U", 100, -3*SMEAR, 3*SMEAR);
    TH1D* hResV = new TH1D("hResV", "Residuals V", 100, -3*SMEAR, 3*SMEAR);

    // Fill histograms
    for (double val : residualsU) hResU->Fill(val);
    for (double val : residualsV) hResV->Fill(val);

    // Draw
    TCanvas* c1 = new TCanvas("c1", "Residuals", 1200, 600);
    c1->Divide(2,1);
    c1->cd(1);
    hResU->Draw();
    c1->cd(2);
    hResV->Draw();

    // Update and run
    //c1->Update();
    //app.Run();
    c1->SaveAs("residuals.png");

    // calculate mean and standard deviation
    const double meanU = calculateMean(residualsU);
    const double sdevU = calculateStandardDeviation(residualsU, meanU);
    const double meanV = calculateMean(residualsV);
    const double sdevV = calculateStandardDeviation(residualsV, meanV);

    std::cout << "mean residual(U): " << meanU << " +/- " << sdevU << "\n"
              << "mean residual(V): " << meanV << " +/- " << sdevV << std::endl;


    return 0;
}
