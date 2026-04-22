# Droplet Impact on a Superheated Surface

## Overview

This case simulates a liquid droplet impacting a heated wall at a temperature above the Leidenfrost point. The model includes coupling between the liquid and solid regions, accounting for:

- Conjugate Heat Transfer (CHT)
- Interfacial evaporation

To reduce computational cost, symmetry is applied, and only **one quarter of the domain** is simulated.

---

## Computational Cost

With the default setup, the simulation requires approximately:

- **Runtime:** ~11 hours (Dual AMD EPYC 7351)
- **Storage:** ~140 GB

---

## ⚠️ Important Note

At the start of the simulation, a warning appears indicating that the keyword `accelerationModel` is outdated and has been replaced by `accelerationForceModel`.

**Do NOT update this keyword.**
Changing to the newer name causes the droplet to remain static for unknown reasons.

---

## Validation

The impact conditions and thermophysical properties used in this case are consistent with:

> Park, J., & Kim, D. E. (2020).
> *Dynamics of liquid drops levitating on superheated surfaces*.
> International Journal of Thermal Sciences, 152, 106321.

---

## 📄 Related Work

The results of this simulation, along with additional cases, are published in:

> Nabbout, K., & Sommerfeld, M. (2026).
> *Numerical simulations of single droplet impact on a hot substrate with geometric VOF and conjugate heat transfer*.
> International Journal of Heat and Fluid Flow, 120, 110346.
> https://doi.org/10.1016/j.ijheatfluidflow.2026.110346

---

## 📄 OpenFOAM version

The original case was provided by Kaissar Nabbout and ran with OpenFOAM-v2112.

The case provided here and the automatic mesh refinement, multiDimAMR, has been
updated to work with OpenFOAM-v2512.
