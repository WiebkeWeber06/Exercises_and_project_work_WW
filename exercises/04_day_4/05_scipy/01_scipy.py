# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 13:58:45 2026

@author: wiewe372
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar

# ----------------------
# LOAD DATA
# ----------------------
exp = np.load("I_q_IPA_exp.npy")
model = np.load("I_q_IPA_model.npy")

print("exp shape:", exp.shape)
print("model shape:", model.shape)
print("First 5 exp rows:\n", exp[:5])
print("First 5 model rows:\n", model[:5])
print("Arrays exactly equal:", np.array_equal(exp, model))

# ----------------------
# EXTRACT COLUMNS
# ----------------------
q_exp, I_exp = exp[:, 0], exp[:, 1]
q_model, I_model = model[:, 0], model[:, 1]

# ----------------------
# REMOVE NaNs / invalid values
# ----------------------
mask_exp = np.isfinite(q_exp) & np.isfinite(I_exp)
q_exp = q_exp[mask_exp]
I_exp = I_exp[mask_exp]

mask_model = np.isfinite(q_model) & np.isfinite(I_model)
q_model = q_model[mask_model]
I_model = I_model[mask_model]

print("Valid exp points:", len(q_exp))
print("Valid model points:", len(q_model))

# ----------------------
# INTERPOLATE MODEL TO EXPERIMENTAL q
# ----------------------
interp_model = interp1d(
    q_model,
    I_model,
    kind="linear",
    bounds_error=False,
    fill_value=np.nan
)

I_model_interp = interp_model(q_exp)

# Keep only overlapping q-range where interpolation is valid
mask_overlap = np.isfinite(I_model_interp)
q_fit = q_exp[mask_overlap]
I_exp_fit = I_exp[mask_overlap]
I_model_fit = I_model_interp[mask_overlap]

print("Points used for fitting:", len(q_fit))

# ----------------------
# FIT SCALING FACTOR: LINEAR SPACE
# ----------------------
def error_linear(scale):
    return np.sum((I_exp_fit - scale * I_model_fit) ** 2)

result_linear = minimize_scalar(error_linear)
best_scale_linear = result_linear.x
I_model_scaled_linear = best_scale_linear * I_model_fit

print(f"Best linear scaling factor: {best_scale_linear:.4e}")
print("Max abs residual (linear fit):", np.max(np.abs(I_exp_fit - I_model_scaled_linear)))
print("Mean abs residual (linear fit):", np.mean(np.abs(I_exp_fit - I_model_scaled_linear)))

# ----------------------
# FIT SCALING FACTOR: LOG SPACE
# ----------------------
def error_log(scale):
    if scale <= 0:
        return np.inf

    mask = (I_exp_fit > 0) & (I_model_fit > 0)
    return np.sum(
        (np.log(I_exp_fit[mask]) - np.log(scale * I_model_fit[mask])) ** 2
    )

result_log = minimize_scalar(error_log)
best_scale_log = result_log.x
I_model_scaled_log = best_scale_log * I_model_fit

print(f"Best log-space scaling factor: {best_scale_log:.4e}")
print("Max abs residual (log fit):", np.max(np.abs(I_exp_fit - I_model_scaled_log)))
print("Mean abs residual (log fit):", np.mean(np.abs(I_exp_fit - I_model_scaled_log)))

# ----------------------
# DIAGNOSTICS
# ----------------------
residuals_linear = I_exp_fit - I_model_scaled_linear
ratio_linear = I_exp_fit / I_model_scaled_linear

residuals_log = I_exp_fit - I_model_scaled_log
ratio_log = I_exp_fit / I_model_scaled_log

# ----------------------
# PLOTTING (4 PANELS)
# ----------------------
fig, axes = plt.subplots(4, 1, figsize=(8, 13), sharex=False)

# ----------------------
# PANEL 1: DATA + FITS
# ----------------------
ax = axes[0]

ax.scatter(q_exp, I_exp, s=15, label="Experimental")
ax.plot(q_fit, I_model_scaled_linear, linewidth=2, label="Linear fit")
ax.plot(q_fit, I_model_scaled_log, linestyle="--", linewidth=2, label="Log-space fit")

ax.set_ylabel("I(q)")
ax.set_title("Model fit to experimental data")
ax.set_xscale("log")
ax.set_yscale("log")
ax.legend()

# ----------------------
# PANEL 2: RESIDUALS
# ----------------------
ax = axes[1]

ax.plot(q_fit, residuals_linear, marker="o", label="Linear fit")
ax.plot(q_fit, residuals_log, marker="o", linestyle="--", label="Log-space fit")
ax.axhline(0, linestyle="--")

ax.set_ylabel("Residuals\n(exp - model)")
ax.set_xscale("log")
ax.legend()

# ----------------------
# PANEL 3: RATIO
# ----------------------
ax = axes[2]

ax.plot(q_fit, ratio_linear, marker="o", label="Linear fit")
ax.plot(q_fit, ratio_log, marker="o", linestyle="--", label="Log-space fit")
ax.axhline(1, linestyle="--")

ax.set_ylabel("Ratio\n(exp / model)")
ax.set_xlabel("Scattering vector q")
ax.set_xscale("log")
ax.legend()

# ----------------------
# PANEL 4: ORIGINAL MODEL
# ----------------------
ax = axes[3]

ax.plot(q_model, I_model, linestyle=":", linewidth=2, label="Original model")

ax.set_title("Original model (unscaled)")
ax.set_xlabel("Scattering vector q")
ax.set_ylabel("I(q)")
ax.set_xscale("log")
ax.set_yscale("log")
ax.legend()

# ----------------------
# FINAL TOUCHES
# ----------------------
plt.tight_layout()

# SAVE FIGURE
plt.savefig(
    "model_fit_comparison.png",
    dpi=300,           # high resolution
    bbox_inches="tight"
)

plt.show()





















