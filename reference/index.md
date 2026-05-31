# Package index

## Model fitting

The unified formula interface and the experimental-design variant.

- [`scm_fit()`](https://yo5uke.com/coresynth/reference/scm_fit.md) : Fit
  a Synthetic Control Method Model
- [`pred()`](https://yo5uke.com/coresynth/reference/pred.md) : Predictor
  Specification for SCM
- [`scm_design()`](https://yo5uke.com/coresynth/reference/scm_design.md)
  : Experimental Synthetic Control Design

## Inference

Permutation, bootstrap, jackknife, parametric, and conformal inference
for the supported estimators.

- [`conformal_inference()`](https://yo5uke.com/coresynth/reference/conformal_inference.md)
  : Conformal Inference for Synthetic Control Estimators
- [`mspe_ratio_pval()`](https://yo5uke.com/coresynth/reference/mspe_ratio_pval.md)
  : Permutation Inference via MSPE Ratio for SCM
- [`sdid_inference()`](https://yo5uke.com/coresynth/reference/sdid_inference.md)
  : Inference for Synthetic Difference-in-Differences
- [`gsc_boot()`](https://yo5uke.com/coresynth/reference/gsc_boot.md) :
  Parametric Bootstrap Inference for GSC (Xu 2017 §3)
- [`gsc_inference()`](https://yo5uke.com/coresynth/reference/gsc_inference.md)
  : Non-parametric Inference for GSC (Xu 2017)
- [`si_inference()`](https://yo5uke.com/coresynth/reference/si_inference.md)
  : Non-parametric Inference for SI (Agarwal et al. 2025)

## Augmentation & broom tidiers

Augmented SCM and one-row tidiers for inference objects.

- [`augment_scm()`](https://yo5uke.com/coresynth/reference/augment_scm.md)
  : Augmented Synthetic Control Method (Ridge ASCM)
- [`tidy(`*`<coresynth_inference>`*`)`](https://yo5uke.com/coresynth/reference/tidy.coresynth_inference.md)
  : Tidy an inference result
- [`glance(`*`<coresynth_inference>`*`)`](https://yo5uke.com/coresynth/reference/glance.coresynth_inference.md)
  : Glance at an inference result

## Visualization & export

- [`plot(`*`<coresynth>`*`)`](https://yo5uke.com/coresynth/reference/plot.coresynth.md)
  : Plot a coresynth model
- [`plot(`*`<scm_design>`*`)`](https://yo5uke.com/coresynth/reference/plot.scm_design.md)
  : Plot an scm_design object
- [`export_json()`](https://yo5uke.com/coresynth/reference/export_json.md)
  : Export coresynth Results to JSON

## Package overview

- [`coresynth`](https://yo5uke.com/coresynth/reference/coresynth-package.md)
  [`coresynth-package`](https://yo5uke.com/coresynth/reference/coresynth-package.md)
  : coresynth: Fast and Unified Synthetic Control Methods

## Low-level C++ routines

Exported RcppArmadillo workhorses called internally by the wrappers.
Most users do not need these directly.

- [`scm_weights_cpp()`](https://yo5uke.com/coresynth/reference/scm_weights_cpp.md)
  : SCM Outer Weights (Joint Optimization of W and V)
- [`scm_inner_weights_cpp()`](https://yo5uke.com/coresynth/reference/scm_inner_weights_cpp.md)
  : SCM Inner Weights (QP Given V)
- [`scm_placebo_cpp()`](https://yo5uke.com/coresynth/reference/scm_placebo_cpp.md)
  : Fast Leave-One-Out Placebo Test for SCM (Abadie et al. 2010)
- [`sdid_unit_weights_cpp()`](https://yo5uke.com/coresynth/reference/sdid_unit_weights_cpp.md)
  : Calculate SDID Unit Weights (omega)
- [`sdid_time_weights_cpp()`](https://yo5uke.com/coresynth/reference/sdid_time_weights_cpp.md)
  : Calculate SDID Time Weights (lambda)
- [`sdid_estimate_cpp()`](https://yo5uke.com/coresynth/reference/sdid_estimate_cpp.md)
  : Calculate SDID Estimate (tau_sdid)
- [`sdid_placebo_cpp()`](https://yo5uke.com/coresynth/reference/sdid_placebo_cpp.md)
  : Fast Placebo Test for SDID
- [`gsc_ife_cpp()`](https://yo5uke.com/coresynth/reference/gsc_ife_cpp.md)
  : Fast Interactive Fixed Effects (IFE) for Generalized Synthetic
  Control
- [`si_pcr_cpp()`](https://yo5uke.com/coresynth/reference/si_pcr_cpp.md)
  : SI-PCR: Synthetic Interventions via Principal Component Regression
- [`tensor_unfold_cpp()`](https://yo5uke.com/coresynth/reference/tensor_unfold_cpp.md)
  : Tensor Unfolding (Matricization) for Synthetic Interventions
- [`soft_impute_cpp()`](https://yo5uke.com/coresynth/reference/soft_impute_cpp.md)
  : Fast Matrix Completion using Soft-Impute Algorithm
- [`kalman_smoother_cpp()`](https://yo5uke.com/coresynth/reference/kalman_smoother_cpp.md)
  : Kalman Filter and RTS Smoother (TASC)
