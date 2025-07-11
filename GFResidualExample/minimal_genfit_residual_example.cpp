#include <Track.h>
#include <TrackPoint.h>
#include <SpacepointMeasurement.h>
#include <DetPlane.h>
#include <RKTrackRep.h>
#include <KalmanFitterRefTrack.h>
#include <FieldManager.h>
#include <ConstField.h>
#include <MaterialEffects.h>
#include <TGeoMaterialInterface.h>
#include <Exception.h>

#include <TROOT.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TMatrixDSym.h>
#include <TStyle.h>

#include <iostream>
#include <vector>
#include <memory>

// === Global constants ===
constexpr double SMEAR = .5;
constexpr int N_POINTS = 15;
constexpr int N_RUNS = 5000;


// === Main part ===
// Function that runs the fit once and fills the histograms
void run_one_fit(
    TRandom3& rnd,
    TH1D* hResU,
    TH1D* hResV,
    TH1D* hSigmaU,
    TH1D* hSigmaV,
    TH1D* hPullU,
    TH1D* hPullV,
    TH1D* hPValue,
    TH1D* hChi2,
    TH1D* hNDF,
    TH1D* hChi2NDF,
    int& n_failed
) {
    //--- Generate reference track ---//
    // Create random reference direction and points
    TLorentzVector refVec;
    refVec.SetXYZT(
        rnd.Uniform(5, 10),
        rnd.Uniform(5, 10),
        rnd.Uniform(5, 10),
        rnd.Uniform(0.01, 0.2)
    );

    const TVector3 refDir = refVec.Vect().Unit();
    const TVector3 refStart(0, 0, 5);

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

    // Fill track with TrackPoints/SpacepointMeasurements
    for (size_t i = 0; i < measPoints.size(); ++i) {
        TVectorD xyz(3);
        xyz[0] = measPoints.at(i).X();
        xyz[1] = measPoints.at(i).Y();
        xyz[2] = measPoints.at(i).Z();

        TMatrixDSym cov(3);
        cov.UnitMatrix();
        cov *= (SMEAR*SMEAR);

        // Create SpacepointMeasurement
        auto* meas = new genfit::SpacepointMeasurement(xyz, cov, 0, i, nullptr);

        track.insertPoint(new genfit::TrackPoint(meas, &track));
    }

    //--- Track Fitting ---//
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

    // Run the fit
    genfit::KalmanFitterRefTrack fitter;
    try {
        fitter.processTrack(&track);
        auto* fitStatus = track.getFitStatus(rep);
        if (!fitStatus || !fitStatus->isFitConverged()) {
            n_failed++;
            return;
        }
        // Store GenFit fit quality parameters in histograms
        hChi2->Fill(fitStatus->getChi2());
        hNDF->Fill(fitStatus->getNdf());
        if (fitStatus->getNdf() > 0) {
            hChi2NDF->Fill(fitStatus->getChi2() / fitStatus->getNdf());
        }
    } catch (const genfit::Exception& e) {
        n_failed++;
        return;
    }

    //--- Readout on Reference Plane ---//
    // Get first and last measured points from the track
    auto* firstTrackPoint = track.getPoint(0);
    auto* lastTrackPoint = track.getPoint(track.getNumPoints() - 1);

    const genfit::AbsMeasurement* firstMeas = firstTrackPoint->getRawMeasurement();
    const genfit::AbsMeasurement* lastMeas = lastTrackPoint->getRawMeasurement();
    TVectorD coordsFirst = firstMeas->getRawHitCoords();
    TVectorD coordsLast = lastMeas->getRawHitCoords();
    TVector3 firstPoint(coordsFirst[0], coordsFirst[1], coordsFirst[2]);
    TVector3 lastPoint(coordsLast[0], coordsLast[1], coordsLast[2]);

    // Plane origin: firstPoint (measured), normal: (lastPoint - firstPoint)
    TVector3 planeOrigin = firstPoint;
    TVector3 planeNormal = (lastPoint - firstPoint).Unit();
    auto plane = std::make_shared<genfit::DetPlane>(planeOrigin, planeNormal);

    // Extrapolate fitted state to this plane
    auto* fittedRep = track.getCardinalRep();
    auto state = track.getFittedState(0);
    try {
        fittedRep->extrapolateToPlane(state, plane, true, false);
    } catch (const genfit::Exception& e) {
        return;
    }
    const double pValue = track.getFitStatus(fittedRep)->getPVal();

    const TVector3 fitPos = state.getPos();
    const TMatrixDSym cov6D = fittedRep->get6DCov(state);
    
    //--- Calculation of Residuals and Pulls ---//
    // Calculate the intersection point of the reference (simulated) track with the detector plane.
    // The reference track is defined as: r(t) = refStart + t * refDir
    // The plane is defined by a point (planeOrigin) and a normal vector (planeNormal).
    // The intersection occurs at the value of t where (r(t) - planeOrigin) . planeNormal == 0
    // Solve for t:
    //   (refStart + t * refDir - planeOrigin) . planeNormal = 0
    //   (refStart - planeOrigin) . planeNormal + t * (refDir . planeNormal) = 0
    //   t = (planeOrigin - refStart) . planeNormal / (refDir . planeNormal)
    const double refDir_dot_planeNormal = refDir.Dot(planeNormal);
    if (std::abs(refDir_dot_planeNormal) < 1e-10) {
        // Track is parallel to plane; skip this event
        return;
    }
    const double intersectionParam = (planeOrigin - refStart).Dot(planeNormal) / refDir_dot_planeNormal;
    const TVector3 truePosInPlane = refStart + intersectionParam * refDir; // Intersection point

    // Extract the 3x3 position block from the 6x6 covariance
    TMatrixDSym posCov(3);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            posCov(i, j) = cov6D(i, j);
        }
    }

    // Get the U and V unit vectors from your DetPlane
    const TVector3 planeU = plane->getU().Unit(); // U direction of the plane
    const TVector3 planeV = plane->getV().Unit(); // V direction of the plane

    // Convert to TVectorD for ROOT math
    TVectorD uvec(3); uvec[0] = planeU.X(); uvec[1] = planeU.Y(); uvec[2] = planeU.Z();
    TVectorD vvec(3); vvec[0] = planeV.X(); vvec[1] = planeV.Y(); vvec[2] = planeV.Z();

    // Project residual and uncertainty into plane U/V axes
    // residual: position on readout plane (distance to origin of plane), which
    // will be projected to U/V; histogramming the absolute value pos-O would be
    // a Raleight distribution (no negative values possible)
    TVector3 resVec = fitPos - truePosInPlane;
    const double resU = resVec.Dot(planeU);
    const double resV = resVec.Dot(planeV);

    // Calculate the reconstructed track error of the track
    // Project the 3D position covariance onto the U and V axes of the plane:
    // varU = u^T * posCov * u, where u is the U unit vector of the plane.
    const double varU = uvec * (posCov * uvec); // u^T * Cov * u
    const double varV = vvec * (posCov * vvec); // v^T * Cov * v

    // Standard deviations (error along U/V)
    const double sigmaU = std::sqrt(varU);
    const double sigmaV = std::sqrt(varV);

    // pull: used to assess the consistency of a fitted parameter with its
    // expected value; the difference between the fitted value and the "true"
    // value, divided by the uncertainty of the fitted value
    const double pullU = resU/sigmaU;
    const double pullV = resV/sigmaV;

    hPValue->Fill(pValue);
    hResU->Fill(resU);
    hResV->Fill(resV);
    hSigmaU->Fill(sigmaU);
    hSigmaV->Fill(sigmaV);
    hPullU->Fill(pullU);
    hPullV->Fill(pullV);
}

int main() {
    // Initialize GenFit global services only once
    genfit::FieldManager::getInstance()->init(new genfit::ConstField(0., 0., 0.));
    genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
    genfit::MaterialEffects::getInstance()->setNoEffects();

    TRandom3 rnd(0); // single RNG instance

    // Histograms for GenFit fit quality and output distributions
    TH1D* hResU    = new TH1D("hResU",    "Residuals U", 100, -3*SMEAR, 3*SMEAR);
    TH1D* hResV    = new TH1D("hResV",    "Residuals V", 100, -3*SMEAR, 3*SMEAR);
    TH1D* hSigmaU  = new TH1D("hSigmaU",  "Fit Sigma U", 100, 0, 2*SMEAR);
    TH1D* hSigmaV  = new TH1D("hSigmaV",  "Fit Sigma V", 100, 0, 2*SMEAR);
    TH1D* hPullU   = new TH1D("hPullU",   "Pulls U",     100, -10, 10);
    TH1D* hPullV   = new TH1D("hPullV",   "Pulls V",     100, -10, 10);
    TH1D* hPValue  = new TH1D("hPValue",  "P-Value",     100, 0, 1);
    TH1D* hChi2    = new TH1D("hChi2",    "Chi2",        100, 0, 100);
    TH1D* hNDF     = new TH1D("hNDF",     "NDF",         100, 0, 50);
    TH1D* hChi2NDF = new TH1D("hChi2NDF", "Chi2/NDF",    100, 0, 10);

    int n_failed = 0;

    // --- Main loop ---
    for (int i = 0; i < N_RUNS; ++i) {
        run_one_fit(rnd, hResU, hResV, hSigmaU, hSigmaV, hPullU, hPullV, hPValue, hChi2, hNDF, hChi2NDF, n_failed);
    }

    gROOT->SetBatch(kTRUE); // no GUI
    gStyle->SetOptFit(1); // Show fit parameters in stats box

    // Fit the residual and pull histograms with a Gaussian
    hResU->Fit("gaus", "Q");
    hResV->Fit("gaus", "Q");
    hPullU->Fit("gaus", "Q");
    hPullV->Fit("gaus", "Q");
    TF1* fitResU = hResU->GetFunction("gaus");
    TF1* fitResV = hResV->GetFunction("gaus");
    TF1* fitPullU = hPullU->GetFunction("gaus");
    TF1* fitPullV = hPullV->GetFunction("gaus");
    if (fitResU && fitResV && fitPullU && fitPullV) {
        std::cout << "Residual U: mean = " << fitResU->GetParameter(1)
                  << ", sigma = " << fitResU->GetParameter(2) << std::endl;
        std::cout << "Residual V: mean = " << fitResV->GetParameter(1)
                  << ", sigma = " << fitResV->GetParameter(2) << std::endl;
        std::cout << "Pull U: mean = " << fitPullU->GetParameter(1)
                  << ", sigma = " << fitPullU->GetParameter(2) << std::endl;
        std::cout << "Pull V: mean = " << fitPullV->GetParameter(1)
                  << ", sigma = " << fitPullV->GetParameter(2) << std::endl;
    }
    std::cout << "Number of failed fits: " << n_failed << " out of " << N_RUNS << std::endl;

    // PNG output for quick look
    TCanvas* c1 = new TCanvas("c1", "Residuals", 1200, 600);
    c1->Divide(2,1);
    c1->cd(1);
    hResU->Draw();
    if (fitResU) fitResU->SetLineColor(kRed);
    if (fitResU) fitResU->Draw("same");
    c1->cd(2);
    hResV->Draw();
    if (fitResV) fitResV->SetLineColor(kRed);
    if (fitResV) fitResV->Draw("same");

    c1->SaveAs("residuals.png");

    TCanvas* c2 = new TCanvas("c2", "Pulls", 1200, 600);
    c2->Divide(2,1);
    c2->cd(1);
    hPullU->Draw();
    if (fitPullU) fitPullU->SetLineColor(kRed);
    if (fitPullU) fitPullU->Draw("same");
    c2->cd(2);
    hPullV->Draw();
    if (fitPullV) fitPullV->SetLineColor(kRed);
    if (fitPullV) fitPullV->Draw("same");
    c2->SaveAs("pulls.png");

    // Store all histograms in a ROOT file
    TFile fout("fit_results.root", "RECREATE");
    hResU->Write();
    hResV->Write();
    hSigmaU->Write();
    hSigmaV->Write();
    hPullU->Write();
    hPullV->Write();
    hPValue->Write();
    hChi2->Write();
    hNDF->Write();
    hChi2NDF->Write();
    fout.Close();

    return 0;
}
