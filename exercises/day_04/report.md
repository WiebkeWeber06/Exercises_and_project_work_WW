# Model Fitting of Scattering Data

## Overview

In this exercise, an experimental scattering dataset was compared to a theoretical model. The goal was to determine whether the model can describe the experimental data when scaled appropriately.

Two datasets were provided:

* `I_q_IPA_exp.npy`: experimental scattering intensity as a function of scattering vector ( q )
* `I_q_IPA_model.npy`: theoretical model of the same quantity

---

## Methods

### Data Preprocessing

* Both datasets were loaded using `numpy.load()`
* Invalid values (`NaN`) in the experimental dataset were removed
* The model data was interpolated onto the experimental ( q )-grid using `scipy.interpolate.interp1d`

### Model Fitting

A scaling factor ( a ) was introduced:

[
I_{\text{model, scaled}}(q) = a \cdot I_{\text{model}}(q)
]

Two fitting approaches were used:

#### 1. Linear Least Squares Fit

Minimizes:

[
\sum (I_{\text{exp}} - a \cdot I_{\text{model}})^2
]

#### 2. Log-space Fit

Minimizes:

[
\sum (\log I_{\text{exp}} - \log(a \cdot I_{\text{model}}))^2
]

This approach accounts for the large dynamic range of scattering data.

### Optimization

* Both fits were performed using `scipy.optimize.minimize_scalar`

---

## Results

### Scaling Factors

* Linear fit: ( a = 1.155 \times 10^{-4} )
* Log-space fit: ( a = 1.239 \times 10^{-4} )

### Fit Quality

To evaluate the fits, three diagnostic plots were generated:

1. **Data + Fit**
2. **Residuals** (( I_{\text{exp}} - I_{\text{model}} ))
3. **Ratio** (( I_{\text{exp}} / I_{\text{model}} ))

---

## Interpretation

### Agreement Between Model and Data

* The model captures the **overall trend** of the experimental data
* Best agreement is observed in the **mid-q region**

### Systematic Deviations

* **Low-q region**: model strongly underestimates the experimental intensity
* **High-q region**: model fails to reproduce detailed peak structure

These deviations are visible in both residual and ratio plots.

### Linear vs Log-space Fit

* The log-space fit provides a slightly improved agreement across the full ( q )-range
* This is reflected in:

  * lower mean residuals
  * ratio closer to 1 over a broader region

However:

> The improvement is modest, and systematic discrepancies remain.

---

## Conclusion

A simple multiplicative scaling improves agreement between the model and experiment but is insufficient for a complete description of the data.

Key findings:

* The model is **partially valid**, capturing general trends
* It fails to describe:

  * low-q behavior
  * detailed peak structure at higher q
* Log-space fitting provides a more balanced comparison but does not resolve fundamental mismatches

---

## Outlook

Further improvements could include:

* Adding a background offset:
  [
  I(q) = a \cdot I_{\text{model}}(q) + b
  ]
* Refining the theoretical model
* Including additional physical effects in the model

---

## Tools Used

* NumPy (data handling)
* SciPy (interpolation and optimization)
* Matplotlib (visualization)

---

## Author

Wiebke Weber

