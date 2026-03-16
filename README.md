# Estimating Future State Water Project Deliveries to Southern California

> A simplified reservoir-network framework for projecting SWP water deliveries
> under climate-driven inflow reduction scenarios (2026–2100).

**Author:** Bharadwaj Vangipuram  
**Course:** Water Resource Systems Engineering  
**Date:** March 2026  

---

## Overview

Southern California receives roughly one-quarter of its total water supply
from the State Water Project (SWP), which draws on eight major reservoirs in
the Sacramento–San Joaquin River Basin (SSJRB). Climate change is projected
to reduce Sierra Nevada snowpack by 33–70%, directly threatening the
reliability of these upstream inflows.

This repository implements a **two-step analytical framework** to estimate how
potential future inflow reductions might affect SWP deliveries to Southern
California:

- **Step 1** establishes the historical statistical relationship between
  reservoir supply conditions (median inflow + September carryover storage)
  and observed reservoir outflows using 46 water years of CDEC data
  (1980–2026). The system-level Pearson correlation is **r = 0.961**
  (R² = 0.924), with the fitted regression equation **Y = 33.13 + 0.6148 X**.

- **Step 2** generates synthetic future inflow sequences using the
  **Thomas–Fiering stochastic streamflow model**, applies inflow-reduction
  scenarios ranging from −10% to −50%, and feeds the resulting supply values
  through the Step-1 regression to predict system outflows for **2026–2100**.

The entire analysis — from raw CDEC download to publication-ready figures —
runs from a single Python script.

---

## Methodology

### Step 1 — Establishing Correlation

| Symbol | Description | Units |
|--------|-------------|-------|
| I_t | Water-year median daily inflow | AF/day |
| C_t | September carryover storage / 365.25 | AF/day |
| Q_out,t | Water-year median daily outflow | AF/day |
| X | Supply variable = I_t + C_t | AF/day |
| Y | Demand proxy = Q_out,t | AF/day |

**Eq. 1 — Pearson Correlation:**

```
r = Cov(I_t + C_t, Q_out,t) / (σ_{I+C} · σ_{Q_out})
```

**Eq. 2 — System-Level Linear Fit:**

```
Y = c + a·X      where c = 33.13,  a = 0.6148
```

Reservoir operations are approximated by a linear decision rule
(Revelle, Joeres & Kirby, 1969), which justifies the use of Pearson
correlation and ordinary least-squares regression.

---

### Step 2 — Scenario-Based Future Deliveries

**Eq. 3 — Thomas–Fiering Lognormal Recursion:**

```
log(I_{t+1}) = μ_{t+1}
               + ρ · (σ_{t+1}/σ_t) · (log(I_t) − μ_t)
               + σ_{t+1} · √(1 − ρ²) · Z_t
```

The lognormal formulation guarantees non-negative synthetic flows.
Parameters (μ, σ, ρ) are estimated from all available monthly observations
per reservoir, handling incomplete records gracefully.

**Eq. 4 — Inflow Reduction Scenarios:**

```
I_scenario = (1 − δ) · I_t
```

Scenarios: Baseline (δ = 0%), −10%, −20%, −40%, −50%.

**Eq. 5 — Aggregated Future System Supply:**

```
X_t = Σ_r ( I_{r,t}^scenario + C_{r,t} )
```

Carryover storage C is held at its historical median (treated as a
management-controlled parameter independent of climate).

---

## Study Reservoirs

| Code | Reservoir | System Role |
|------|-----------|-------------|
| SHA | Shasta | Largest CVP reservoir; primary Sacramento River storage |
| CLE | Trinity | Trinity River diversion; feeds Sacramento River |
| ORO | Oroville | Largest SWP reservoir on the Feather River |
| FOL | Folsom | American River CVP storage |
| BER | Berryessa | Putah Creek; short inflow record ⚠️ |
| NML | New Melones | Stanislaus River; joint SWP/CVP use |
| DNP | Don Pedro | Tuolumne River; CVP ties |
| SNL | San Luis | Off-stream reservoir; pumped from the Delta ⚠️ |

> ⚠️ BER and SNL have incomplete historical inflow records. Thomas–Fiering
> parameters for these two reservoirs carry higher uncertainty; results should
> be interpreted with caution.

**CDEC Sensor IDs:** Inflow = 76 (CFS, daily) · Outflow = 23 (CFS, daily) · Storage = 15 (AF, monthly)

---

## Repository Structure

```
swp_southern_california_analysis/
│
├── swp_southern_california_analysis.py   # Main analysis script (all steps)
├── README.md                             # This file
├── requirements.txt                      # Python dependencies
│
└── outputs/                              # Created automatically on first run
    ├── figures/
    │   ├── figure01_SWP_allocation.png / .pdf
    │   ├── figure02_daily_inflows_8panel.png
    │   ├── figure03_daily_outflows_8panel.png
    │   ├── figure04_storage_8panel.png
    │   ├── figure05_pearson_correlation_matrices.png
    │   ├── figure06_individual_reservoir_fits.png
    │   ├── figure07_system_level_fit.png
    │   ├── figure08_TF_historical_vs_synthetic_WY_medians.png
    │   ├── figure09_TF_monthly_means_validation.png
    │   ├── figure10_future_outflow_timeseries.png
    │   └── figure11_future_outflow_statistical_summary.png
    │
    ├── data/
    │   ├── WY_reservoir_summary_all.csv
    │   ├── system_supply_outflow_WY_summary.csv
    │   ├── aggregated_system_fit_summary.csv
    │   ├── individual_reservoir_supply_outflow_fits.csv
    │   ├── pearson_inflow_correlation_matrix.csv
    │   ├── pearson_outflow_correlation_matrix.csv
    │   ├── pearson_carryover_correlation_matrix.csv
    │   ├── thomas_fiering_synthetic_summary_all_reservoirs.csv
    │   ├── thomas_fiering_WY_daily_medians_historical_and_synthetic.csv
    │   ├── historical_median_carryover_storage.csv
    │   ├── future_reservoir_supply_outflow_scenarios.csv
    │   ├── future_system_outflows_scenarios_2026_2100.csv
    │   ├── future_outflow_scenario_summary.csv
    │   └── {STA}_thomas_fiering_synthetic_monthly_inflow.csv  (×8)
    │
    └── workbooks/
        ├── INFLOWS_DAILY_1980-01-01_2026-01-31.xlsx
        ├── OUTFLOWS_DAILY_1980-01-01_2026-01-31.xlsx
        └── STORAGE_MONTHLY_1980-01-01_2026-01-31.xlsx
```

---

## Figures Produced

| Figure | Description | Paper Reference |
|--------|-------------|-----------------|
| 01 | SWP Table A, requested & approved allocations 1996–2025 | Fig. 1 |
| 02 | Daily inflows for all 8 reservoirs with median | Fig. 4 |
| 03 | Daily outflows for all 8 reservoirs with median | Fig. 3 |
| 04 | Monthly storage with September carryover median | Fig. 5 |
| 05 | Pearson correlation matrices (inflow / outflow / storage) | Fig. 6 |
| 06 | Individual reservoir linear fits (4×2 panel) | Fig. 7 |
| 07 | System-level linear fit (r = 0.961, R² = 0.924) | Fig. 8 |
| 08 | TF historical vs synthetic WY-median daily inflows | Fig. 9 |
| 09 | TF monthly mean seasonal validation | — |
| 10 | Future outflow time series 2026–2100 by scenario | Fig. 10 |
| 11 | Statistical summary: fit, mean bar chart, boxplot | Fig. 11 |

---

## Installation

**Python 3.10 or higher is required.**

```bash
# Clone the repository
git clone https://github.com/your-username/swp_southern_california_analysis.git
cd swp_southern_california_analysis

# Install dependencies
pip install -r requirements.txt
```

---

## Usage

### Full run (downloads CDEC and produces all outputs)

```bash
python swp_southern_california_analysis.py
```

Outputs are written to `~/Downloads/SWP_Analysis` by default.

### Custom output folder

```bash
python swp_southern_california_analysis.py --output /path/to/results
```

### Regenerate figures without re-downloading from CDEC

Useful after the first full run when you want to adjust a figure style
without waiting for CDEC again.

```bash
python swp_southern_california_analysis.py --skip-download
```

> **Note:** `--skip-download` requires that `WY_reservoir_summary_all.csv`
> and `system_supply_outflow_WY_summary.csv` already exist in the output folder
> from a previous full run.

### Custom date window

```bash
python swp_southern_california_analysis.py --start 1985-01-01 --end 2025-12-31
```

### All options

```
usage: swp_southern_california_analysis.py [-h] [--output OUTPUT]
                                           [--start START] [--end END]
                                           [--skip-download]

options:
  -h, --help        Show this help message and exit
  --output OUTPUT   Output folder (default: ~/Downloads/SWP_Analysis)
  --start START     Historical start date YYYY-MM-DD (default: 1980-01-01)
  --end END         Historical end date YYYY-MM-DD   (default: 2026-01-31)
  --skip-download   Reuse cached WY summary CSVs instead of re-downloading
```

---

## Data Source

All hydrological data are downloaded from the **California Data Exchange
Center (CDEC)** maintained by the California Department of Water Resources:

```
https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet
```

CDEC is a public portal and no API key is required. The script uses
automatic retry with exponential back-off (up to 8 attempts per annual
chunk) to handle transient server errors.

The SWP Table A allocation data (Figure 1) is hardcoded in the script from
public DWR records and does not require a network request.

---

## Known Limitations

1. **BER and SNL data quality.** Berryessa and San Luis have incomplete
   historical inflow and outflow records. Thomas–Fiering parameters for
   these two reservoirs carry higher uncertainty.

2. **Simplified operational representation.** The Step-1 linear fit does
   not capture the full operational complexity of the SWP–CVP system.
   Environmental flow requirements, Delta export restrictions, flood-control
   rules, and conjunctive-use agreements are not modelled.

3. **Uniform inflow reduction assumption.** Inflow reductions are applied as
   constant fractions across all future years, months, and reservoirs. Real
   climate-driven changes will be spatially and temporally heterogeneous.

4. **Static carryover storage.** Carryover storage is held at historical
   medians. Temperature-driven evaporation increases under warming are not
   accounted for.

---

## Potential Future Work

- Replace constant percentage reductions with regionally downscaled GCM
  projections to represent non-uniform climate responses across watersheds.
- Reconstruct BER and SNL inflows via mass-balance using hourly reservoir
  stage data to improve TF parameter reliability for those reservoirs.
- Incorporate temperature-driven evaporation losses and surface-area dynamics
  to make carryover storage climate-responsive rather than fixed.
- Extend the framework to include Delta export constraints and regulatory
  flow requirements for a more operationally realistic simulation.

---

## Requirements

```
pandas>=1.5
numpy>=1.23
matplotlib>=3.6
requests>=2.28
urllib3>=1.26
scipy>=1.9
openpyxl>=3.0
```

Install all at once:

```bash
pip install -r requirements.txt
```

---

## References

Arnold, W., 2021. *The economic value of carryover storage in California's
water supply system.* University of California, Davis.

Berman, J.J., 2016. *Data simplification: taming information with open source
tools.* Morgan Kaufmann. https://doi.org/10.1016/B978-0-12-803781-2.00004-7

CDEC, 2026. California Data Exchange Center.
https://cdec.water.ca.gov

Draper, A.J. et al., 2004. CalSim: Generalized model for reservoir system
analysis. *J. Water Resour. Plan. Manage.*, 130(6), 480–489.

Li, D. et al., 2024. Uncovering historical reservoir operation rules and
patterns. *Water Resources Research*, 60, e2023WR036686.

Pagán, B.R. et al., 2016. Extreme hydrological changes in the southwestern
US drive reductions in water supply to Southern California by mid-century.
*Environ. Res. Lett.*, 11(9), 094026.

Ray, P. et al., 2020. Vulnerability and risk: climate change and water supply
from California's Central Valley water system. *Clim. Change*, 161, 177–199.

Revelle, C., Joeres, E. & Kirby, W., 1969. The linear decision rule in
reservoir management and design. *Water Resour. Res.*, 5(4), 767–777.

SCAG, 2025. *Southern California economic update.* Southern California
Association of Governments.

Steinschneider, S. et al., 2023. Uncertainty decomposition to understand
the influence of water systems model error. *Water Resour. Res.*, 59,
e2022WR032349.

Sunding, D., Browne, D. & Zhu, A., 2023. *The economy of the State Water
Project.* California DWR.

---

## License

This project is released for academic and research purposes.  
If you use this code or methodology in your own work, please cite the paper:

> Vangipuram, B., 2026. Estimating Future State Water Project Deliveries to
> Southern California Using Synthetic Inflow Simulation and System-Level
> Supply–Outflow Relationships. *Water Resource Systems Engineering.*
