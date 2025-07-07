# Documentation: Residuals, Pulls, and Error Analysis in GenFit

The goal of this code is to provide a simulation-based study of track fitting performance using GenFit, focusing on residual and pull distributions.
The implementation mimics realistic detector scenarios in which spacepoint measurements are obtained along a particle trajectory, but deliberately omits specific detector geometry, magnetic field, and material effects for clarity and pedagogical value.
This kind of setup matches the simulation studies presented in the [original GenFit paper][hoeppner:2010], which are designed to validate the mathematical and statistical correctness of Kalman filter-based fitting.

## **Overview**

- **Measurement Generation:**  
  The code generates a set of randomly spaced reference points along a random straight-line trajectory.
  Measurements are then smeared around these points with isotropic Gaussian noise.
  This mimics the process of measurement in a TPC or similar detector, though in a simplified way (for instance, the $x$ and $y$ resolutions are not better than $z$, as would be in a TPC).

- **Track Fitting:**  
  The initial track seed is set to the true starting point and direction.
  The Kalman filter then fits the track to the smeared measurements.

- **Residuals and Pulls:**  
  After fitting, the residuals and pulls are calculated at the virtual detector plane.
  The code iterates this process many times, filling histograms for statistical analysis.

- **Statistical Interpretation:**  
  If the fit is able to correctly determine the track parameters and their errors, the pull distributions will be Gaussians of width $\sigma=1$ and of mean value 0.

---

## **Track Representation and State Vector**

### **Track Representation**

The example uses GenFit’s `RKTrackRep` (Runge-Kutta Track Representation) for the track model, adopted from the COMPASS experiment.

The state vector for a straight-line track (as in this example, since no magnetic field is present) is:
```
x_k = (u, v, du/dw, dv/dw)^T
```
where:
- `u` and `v` span the chosen detector plane,
- `w = u × v` is the plane normal.

Although in this example the propagation is trivial (straight lines), using `RKTrackRep` ensures that the code can be generalized to more realistic scenarios with magnetic fields.
It also ensures that the covariance propagation and extrapolation machinery are fully exercised.


## **Detector Plane**

GenFit’s design supports both physical detector planes (such as silicon trackers) and continuous readout geometries (e.g. for TPCs).
In the latter, measurements do not naturally lie in a fixed plane.
To provide a consistent framework for fitting and residual calculation, **virtual detector planes** are introduced.

In this simulation, each measurement is a spacepoint (with no intrinsic plane), but residuals must be evaluated in a well-defined coordinate system.
Thus, a virtual plane is defined at the first reference point of the track, with its normal along the track direction.
Since there is no defined beamline or target in this example, the first point is a natural choice for the readout plane.
For the direction, the unit vector pointing from the first to the last point (i.e. roughly parallel to the track) is a reasonable approximation.

- **Origin:** The first point along the MC truth track (`track.getPoint(0)`).
- **Normal (`n`):** The unit vector from the first to last point, i.e., the track direction.
- **Axes (`U`, `V`):** Chosen orthogonal vectors spanning the plane, constructed from the normal as in GenFit conventions.

---

## **Residual Calculation**

The **residual** is the difference between the fitted and true positions, projected onto the axes of the virtual detector plane:
```
Residual_U = (P_fit - P_truth) ⋅ U
Residual_V = (P_fit - P_truth) ⋅ V
```
- `P_fit`: The position of the fitted track, extrapolated to the virtual plane.
- `P_truth`: The true (unsmeared) position at the plane.
- `U`, `V`: The orthogonal axes of the plane.

The dot product projects the residual onto the chosen axis, so that it is expressed in the local coordinate system of the plane.
This is necessary because any residual/pull can only be meaningfully calculated on a reference plane ("readout plane").

## **Pull Calculation**

The **pull** is defined as:
```
Pull_{U/V} = Residual_{U/V} / σ_{U/V}
```
where `σ_{U/V}` is the uncertainty of the fitted position projected onto the U or V axis.
The projections of the covariance matrices, and thus, σ, onto the (virtual) detector planes have to be calculated.

**Measurement Error Estimation**  
- The full 6x6 covariance matrix of the fitted state is obtained (`get6DCov`). The first 3x3 block represents position uncertainties.
- The variance along the direction of a unit vector (U or V) is:
  ```
  σ^2 = U^T ⋅ Cov_3x3 ⋅ U
  ```
  This is the standard method for projecting a covariance matrix onto a direction. The square root gives σ.

Projecting the full covariance ensures that all correlations are properly accounted for, which is critical for the statistical interpretation of the pull distribution.

---

## **Parameter Studies**

- **Changing the Number of Points:**  
  More measurement points lead to better constraint of the track and smaller uncertainties; fewer points increase the uncertainty and can degrade the fit.

- **Changing the Smearing Parameter:**  
  Larger smearing increases the measurement noise and widens the residual and pull distributions; smaller smearing tightens them.


## **Technical Caveats and Physical Limitations**

- **No Detector Geometry:**  
  No actual detector geometry or material effects are present, so there is no multiple scattering or energy loss.
- **No Magnetic Field:**  
  The field is set to zero, making the track linear.
- **No Physical Units:**  
  The example is unitless; it serves as a demonstration rather than a model of a real experiment.
- **No Outlier or Failure Analysis:**  
  Fit failures are counted but not analyzed further. In a more realistic study, their causes and impact should be investigated.

---

## **References**

1. [GenFit: A Generic Track-Fitting Toolkit, Höppner et.al., NIM A (620), August 2010, Pages 518-525][hoeppner:2010]
2. [Everything you always wanted to know about pulls: CDF Technical Report 5776, August 2002][demortier:2002]
3. [GenFit Source Code](https://github.com/GenFit/GenFit)
4. [Geant4 Physics Processes](https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/TrackingAndPhysics/physicsProcess.html#g4vprocess)

[hoeppner:2010]: https://doi.org/10.1016/j.nima.2010.03.136 "A novel generic framework for track fitting in complex detector systems"
[demortier:2002]: https://hep-physics.rockefeller.edu/luc/technical_reports/cdf5776_pulls.pdf "Everything you always wanted to know about pulls"
