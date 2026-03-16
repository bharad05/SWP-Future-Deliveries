#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
swp_southern_california_analysis.py
=====================================
Estimating Future State Water Project Deliveries to Southern California
Using Synthetic Inflow Simulation and System-Level Supply-Outflow Relationships

Author  : Bharadwaj Vangipuram
Course  : Water Resource Systems Engineering
Date    : March 2026

-------------------------------------------------------------------------------
PAPER CONTEXT
-------------------------------------------------------------------------------
Southern California relies on the State Water Project (SWP) for roughly a
quarter of its total water supply.  Climate change is projected to reduce
Sierra Nevada snowpack by 33-70%, which will directly reduce inflows to the
eight major reservoirs in the SWP-CVP system.  This script implements a
two-step simplified framework to estimate how those inflow reductions might
affect future water deliveries.

Step 1 establishes a statistical relationship between upstream supply conditions
(inflow + carryover storage) and observed reservoir outflows using 46 historical
water years (1980-2026).  The system-level Pearson correlation is r = 0.961
(R² = 0.924), with the fitted linear equation Y = 33.13 + 0.6148 X.

Step 2 uses the Thomas-Fiering stochastic streamflow model to generate synthetic
inflow sequences that preserve the statistical character of the historical record,
then applies inflow-reduction scenarios (-10% to -50%) and feeds the resulting
supply values through the Step-1 regression to predict future system outflows
for the period 2026-2100.

-------------------------------------------------------------------------------
METHODOLOGY EQUATIONS 
-------------------------------------------------------------------------------

Eq. 1 -- Pearson Correlation Coefficient
    r = Cov(I_t + C_t, Q_out,t) / (sigma_{I+C} * sigma_{Q_out})

    where I_t is median reservoir inflow, C_t is September carryover storage,
    Q_out,t is median reservoir outflow, and sigma denotes standard deviation.
    Used to quantify linear relationships between supply and outflow at both
    the individual reservoir and system levels.

Eq. 2 -- System-Level Linear Fit  (Step-1 result, used in Step 2)
    Y = c + a*X

    Y = aggregated system demand (sum of median outflows, all reservoirs) [AF/day]
    X = aggregated system supply (sum of median inflows + carryover)      [AF/day]
    c = 33.13   (regression intercept)
    a = 0.6148  (regression slope)

    Assumption: Reservoir operations can be approximated by a linear decision
    rule (Revelle, Joeres & Kirby, 1969), implying that supply and outflow
    exhibit an approximately linear relationship.

Eq. 3 -- Thomas-Fiering Lognormal Recursion
    log(I_{t+1}) = mu_{t+1}
                   + rho * (sigma_{t+1}/sigma_t) * (log(I_t) - mu_t)
                   + sigma_{t+1} * sqrt(1 - rho^2) * Z_t

    I_{t+1}     = synthetic inflow at time t+1         [AF/month]
    mu_{t+1}    = historical mean of log-inflow, month t+1
    mu_t        = historical mean of log-inflow, month t
    rho         = lag-1 autocorrelation between months t and t+1
    sigma_{t+1} = historical std dev of log-inflow, month t+1
    sigma_t     = historical std dev of log-inflow, month t
    Z_t         = standard normal random variate

    The lognormal formulation guarantees non-negative synthetic flows, which is
    a physical requirement for streamflow.  Parameters are fitted from all
    available monthly observations (incomplete water years are handled gracefully
    by estimating each month's statistics independently).

    Validation: A two-sample Kolmogorov-Smirnov test checks that the synthetic
    and historical daily-equivalent inflow distributions are statistically
    indistinguishable (p > 0.05 indicates no significant difference).

Eq. 4 -- Inflow Reduction Scenarios
    I_scenario = (1 - delta) * I_t

    I_scenario  = reduced synthetic inflow under a given scenario   [AF/day]
    I_t         = baseline synthetic WY-median daily inflow         [AF/day]
    delta       = fractional reduction (0.00, 0.10, 0.20, 0.40, 0.50)

    Scenarios evaluated: Baseline (0%), -10%, -20%, -40%, -50%.
    Assumption: Reductions are uniform across all future years, months, and
    reservoirs.  This is a simplification; real climate-driven changes will be
    spatially and temporally heterogeneous.

Eq. 5 -- Aggregated Future System Supply Variable
    X_t = SUM_r ( I_{r,t}^scenario + C_{r,t} )

    X_t              = system supply for water year t              [AF/day]
    I_{r,t}^scenario = scenario-reduced synthetic inflow, reservoir r [AF/day]
    C_{r,t}          = historical median Sept carryover / 365.25,  reservoir r [AF/day]

    Assumption: Carryover storage (C) is treated as a human-controlled,
    management-driven parameter that remains at its historical median.  Although
    carryover can respond to climate variability (Arnold, 2021; Appendix H1a1,
    2025), this study holds it constant to isolate the effect of inflow changes.

-------------------------------------------------------------------------------
STUDY RESERVOIRS
-------------------------------------------------------------------------------
All eight reservoirs are in the Sacramento-San Joaquin River Basin (SSJRB) and
together represent the major storage and conveyance nodes of the SWP-CVP system:

    Code  Reservoir        Role in system
    ----  ---------------  -------------------------------------------------------
    SHA   Shasta           Largest CVP reservoir; primary Sacramento River storage
    CLE   Trinity          Trinity River diversion; feeds Sacramento River system
    ORO   Oroville         Largest SWP reservoir on the Feather River
    FOL   Folsom           American River CVP storage
    BER   Berryessa        Putah Creek; Lake Berryessa (short inflow record)
    NML   New Melones      Stanislaus River; joint SWP/CVP use
    DNP   Don Pedro        Tuolumne River; primarily local use with CVP ties
    SNL   San Luis         Off-stream reservoir; pumped from the Delta
                            (operates differently from upstream reservoirs --
                             see paper Section 3.1 for discussion of negative
                             inflow/outflow correlations with other reservoirs)

CDEC Sensor IDs (California Data Exchange Center):
    Inflow  = 76  (daily, CFS)
    Outflow = 23  (daily, CFS)
    Storage = 15  (monthly, AF)

-------------------------------------------------------------------------------
KNOWN LIMITATIONS
-------------------------------------------------------------------------------
1.  BER and SNL have incomplete historical inflow and outflow records.
    TF parameters for these two reservoirs carry higher uncertainty, which
    is reflected in more erratic synthetic behavior (see Figure 9).

2.  The Step-1 linear fit does not represent the full operational complexity
    of the SWP-CVP system.  Environmental flow requirements, Delta export
    restrictions, flood-control rules, and conjunctive-use agreements are not
    modelled.

3.  Inflow reductions are applied as constant fractions; real climate responses
    will vary across seasons and decades.

4.  Carryover storage is held at historical medians; temperature-driven
    evaporation increases under warming are not accounted for.

-------------------------------------------------------------------------------
OUTPUTS
-------------------------------------------------------------------------------
All files are written to SAVE_FOLDER (configurable at top of script).

Figures (PNG, 300 dpi; Figure 1 also saved as PDF at 600 dpi):
    figure01_SWP_allocation.png / .pdf
    figure02_daily_inflows_8panel.png
    figure03_daily_outflows_8panel.png
    figure04_storage_8panel.png
    figure05_pearson_correlation_matrices.png
    figure06_individual_reservoir_fits.png
    figure07_system_level_fit.png
    figure08_TF_historical_vs_synthetic_WY_medians.png
    figure09_TF_monthly_means_validation.png
    figure10_future_outflow_timeseries.png
    figure11_future_outflow_statistical_summary.png

CSVs:
    WY_reservoir_summary_all.csv
    system_supply_outflow_WY_summary.csv
    aggregated_system_fit_summary.csv
    individual_reservoir_supply_outflow_fits.csv
    pearson_inflow_correlation_matrix.csv
    pearson_outflow_correlation_matrix.csv
    pearson_carryover_correlation_matrix.csv
    thomas_fiering_synthetic_summary_all_reservoirs.csv
    thomas_fiering_WY_daily_medians_historical_and_synthetic.csv
    historical_median_carryover_storage.csv
    future_reservoir_supply_outflow_scenarios.csv
    future_system_outflows_scenarios_2026_2100.csv
    future_outflow_scenario_summary.csv
    {STA}_thomas_fiering_synthetic_monthly_inflow.csv  (one per reservoir)

Excel workbooks (raw CDEC downloads, one sheet per reservoir):
    INFLOWS_DAILY_{start}_{end}.xlsx
    OUTFLOWS_DAILY_{start}_{end}.xlsx
    STORAGE_MONTHLY_{start}_{end}.xlsx

-------------------------------------------------------------------------------
USAGE
-------------------------------------------------------------------------------
    # Full run (downloads CDEC data and produces all figures)
    python swp_southern_california_analysis.py

    # Write outputs to a custom folder
    python swp_southern_california_analysis.py --output /path/to/results

    # Re-run figures without re-downloading CDEC (uses cached CSVs)
    python swp_southern_california_analysis.py --skip-download

    # Custom date window
    python swp_southern_california_analysis.py --start 1985-01-01 --end 2025-12-31

-------------------------------------------------------------------------------
DEPENDENCIES
-------------------------------------------------------------------------------
    pip install pandas numpy matplotlib requests urllib3 scipy openpyxl

    Tested with Python 3.10+.  No non-standard system libraries required.

-------------------------------------------------------------------------------
REFERENCES
-------------------------------------------------------------------------------
Arnold, W., 2021.  The economic value of carryover storage in California's water
    supply system.  University of California, Davis.

Bardeen, S., 2023.  Can we capture more water in the Delta?
    Public Policy Institute of California.

Berman, J.J., 2016.  Data simplification: taming information with open source
    tools.  Morgan Kaufmann.  doi:10.1016/B978-0-12-803781-2.00004-7

CDEC, 2026.  California Data Exchange Center.  https://cdec.water.ca.gov

Draper, A.J. et al., 2004.  CalSim: Generalized model for reservoir system
    analysis.  J. Water Resour. Plan. Manage., 130(6), 480-489.

Li, D. et al., 2024.  Uncovering historical reservoir operation rules and
    patterns.  Water Resources Research, 60, e2023WR036686.

Pagan, B.R. et al., 2016.  Extreme hydrological changes in the southwestern US
    drive reductions in water supply to Southern California by mid-century.
    Environ. Res. Lett., 11(9), 094026.

Ray, P. et al., 2020.  Vulnerability and risk: climate change and water supply
    from California's Central Valley water system.  Clim. Change, 161, 177-199.

Revelle, C., Joeres, E. & Kirby, W., 1969.  The linear decision rule in
    reservoir management and design.  Water Resour. Res., 5(4), 767-777.

SCAG, 2025.  Southern California economic update.  Southern California
    Association of Governments.

Steinschneider, S. et al., 2023.  Uncertainty decomposition to understand the
    influence of water systems model error.  Water Resour. Res., 59,
    e2022WR032349.

Sunding, D., Browne, D. & Zhu, A., 2023.  The economy of the State Water
    Project.  California DWR.
"""

# ===========================================================================
# Standard library imports
# ===========================================================================
import os
import argparse
from datetime import date
from io import StringIO

# ===========================================================================
# Third-party imports
# ===========================================================================
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from scipy.stats import ks_2samp


# ===========================================================================
# ==========================================================================
#
#   SECTION 1 -- CONFIGURATION
#   All user-editable parameters live here.  Nothing below this section
#   should need to be changed for a standard run.
#
# ===========================================================================
# ===========================================================================

# Output directory.  Created automatically if it does not exist.
SAVE_FOLDER = os.path.join(os.path.expanduser("~"), "Downloads", "SWP_Analysis")

# Historical data download window.
# CDEC records generally begin around 1980 for most study reservoirs.
HIST_START = "1980-01-01"
HIST_END   = "2026-01-31"

# Future projection period (water years).
# The synthetic series is generated to cover this window.
FUTURE_START_WY = 2026
FUTURE_END_WY   = 2100

# Step-1 regression coefficients from the paper (Y = c + a*X, Eq. 2).
# These are hard-coded as fallback values in case the live CDEC download
# does not return enough data to re-fit the model.  When a full download
# is available, these constants are re-estimated from the data and the
# computed values replace these defaults.
STEP1_INTERCEPT_C = 33.13
STEP1_SLOPE_A     = 0.6148
STEP1_R           = 0.961
STEP1_R2          = 0.924

# Thomas-Fiering settings.
# N_SYN_YEARS is the number of years in the synthetic series; the series is
# then trimmed to FUTURE_START_WY through FUTURE_END_WY.
# MIN_POSITIVE_FLOW is the floor applied before the log transform in Eq. 3
# (lognormal TF requires strictly positive values).
N_SYN_YEARS       = 120     # generates enough years to cover 2026-2100
RANDOM_SEED       = 42      # fixed for reproducibility
MIN_POSITIVE_FLOW = 10.0    # AF/month

# Inflow-reduction scenarios (delta values as fractions, Eq. 4).
SCENARIOS = {
    "Baseline":    0.00,
    "Inflow -10%": 0.10,
    "Inflow -20%": 0.20,
    "Inflow -40%": 0.40,
    "Inflow -50%": 0.50,
}

# SWP Table A allocation data (1996-2025).
# Hardcoded so Figure 1 can be produced without any network request.
SWP_DATA = {
    "Year": list(range(1996, 2026)),
    "Table_A": [
        2510200, 2480200, 2505200, 2508200, 2553200, 2554200, 2557200, 2558200,
        2559200, 2569600, 2582800, 2593100, 2593100, 2593100, 2623100, 2623100,
        2623100, 2623100, 2626544, 2629544, 2629544, 2629544, 2629544, 2629544,
        2633544, 2633544, 2633544, 2633544, 2633544, 2633544,
    ],
    "Initial_Request": [
        1087750, 1385150, 1767400, 1807900, 2048836, 2554200, 2345867, 2558200,
        2559200, 2571600, 2582800, 2593100, 2593100, 2593100, 2623100, 2623100,
        2623100, 2623100, 2626544, 2629544, 2629544, 2629544, 2629544, 2629544,
        2633544, 2633544, 2633544, 2633544, 2633544, 2633544,
    ],
    "Approved_Allocation": [
        1062950, 1385150, 1469110, 1475640, 1781261,  510840, 1150740, 1279100,
        1279600, 1541760, 1807960,  648275,  907585,  518620,  393464, 1573860,
        1311550,  918086,       0,  394433,  394433, 1577726,  525909,  394433,
         395033,  131678,  395033,  790064,  395033,  921742,
    ],
}


# ===========================================================================
# ==========================================================================
#
#   SECTION 2 -- PHYSICAL CONSTANTS AND RESERVOIR DEFINITIONS
#
# ===========================================================================
# ===========================================================================

# Unit conversion: 1 CFS of flow sustained for one day equals 1.983471074 AF.
CFS_TO_AF_PER_DAY = 1.983471074

# Conversion used to express carryover storage (AF) on the same AF/day basis
# as median daily inflow, so both can be added in Eq. 5.
DAYS_PER_YEAR = 365.25

# California water year runs October through September.
# The list below maps calendar months to their sequential slot in the
# water year (October = slot 0, September = slot 11).
WY_MONTHS = [10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9]
WY_LABELS = ["Oct", "Nov", "Dec", "Jan", "Feb",
             "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep"]
M_TO_SLOT = {m: i for i, m in enumerate(WY_MONTHS)}

# CDEC data endpoint
CDEC_BASE = "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet"

# Study reservoirs with CDEC station codes and sensor IDs.
#   Sensor 76 = inflow (CFS, daily)
#   Sensor 23 = outflow (CFS, daily)
#   Sensor 15 = storage (AF, monthly)
RESERVOIRS = {
    "SHA": {"name": "Shasta",      "inflow": 76, "outflow": 23, "storage": 15},
    "CLE": {"name": "Trinity",     "inflow": 76, "outflow": 23, "storage": 15},
    "ORO": {"name": "Oroville",    "inflow": 76, "outflow": 23, "storage": 15},
    "FOL": {"name": "Folsom",      "inflow": 76, "outflow": 23, "storage": 15},
    "BER": {"name": "Berryessa",   "inflow": 76, "outflow": 23, "storage": 15},
    "NML": {"name": "New Melones", "inflow": 76, "outflow": 23, "storage": 15},
    "DNP": {"name": "Don Pedro",   "inflow": 76, "outflow": 23, "storage": 15},
    "SNL": {"name": "San Luis",    "inflow": 76, "outflow": 23, "storage": 15},
}
# Ordered list used wherever we iterate and need consistent left-to-right ordering.
ORDER = ["SHA", "CLE", "ORO", "FOL", "BER", "NML", "DNP", "SNL"]


# ===========================================================================
# ==========================================================================
#
#   SECTION 3 -- MATPLOTLIB STYLE
#   A single rcParams block keeps figure styling consistent across all outputs.
#
# ===========================================================================
# ===========================================================================

mpl.rcParams.update({
    "font.family":       "serif",
    "font.size":         11,
    "axes.titlesize":    12,
    "axes.labelsize":    11,
    "xtick.labelsize":   10,
    "ytick.labelsize":   10,
    "legend.fontsize":   9,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.linewidth":    0.8,
})

# Colour palette used consistently across all figures.
C_BLUE   = "#4C72B0"
C_RED    = "#C44E52"
C_GREEN  = "#55A868"
C_PURPLE = "#8172B2"
C_GREY   = "#4D4D4D"
C_SHADE  = "0.85"    # drought-period shading in Figure 1


# ===========================================================================
# ==========================================================================
#
#   SECTION 4 -- NETWORK UTILITIES
#   Functions for downloading data from the CDEC API with retry logic.
#
# ===========================================================================
# ===========================================================================

def make_session():
    """
    Create a requests.Session with automatic retry and exponential back-off.

    CDEC returns transient HTTP errors (429, 5xx) under load.  The session
    retries up to 8 times per chunk, waiting progressively longer between
    attempts (0.8, 1.6, 3.2 ... seconds).  This avoids hammering the server
    while still recovering from short outages.
    """
    session = requests.Session()
    retry = Retry(
        total=8, connect=8, read=8,
        backoff_factor=0.8,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET",),
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=10, pool_maxsize=10)
    session.mount("https://", adapter)
    session.mount("http://",  adapter)
    return session


# A single session is shared across all downloads so connections are reused.
SESSION = make_session()


def _cdec_url(station, sensor, start, end, dur):
    """Build a CDEC CSVDataServlet query URL."""
    return (
        f"{CDEC_BASE}?Stations={station}"
        f"&SensorNums={int(sensor)}"
        f"&dur_code={dur}"
        f"&Start={start}&End={end}"
    )


def _year_chunks(start, end):
    """
    Split a date range into one-calendar-year chunks.

    CDEC becomes unreliable for multi-year daily requests.  Splitting into
    annual chunks keeps each HTTP request small and easy to retry if one
    year fails without losing all the data.
    """
    s = pd.Timestamp(start).date()
    e = pd.Timestamp(end).date()
    for y in range(s.year, e.year + 1):
        cs = max(date(y, 1, 1), s)
        ce = min(date(y, 12, 31), e)
        yield pd.Timestamp(cs).strftime("%Y-%m-%d"), pd.Timestamp(ce).strftime("%Y-%m-%d")


def fetch_cdec_series(station, sensor, start, end, dur, chunk_yearly=True):
    """
    Download a CDEC time series for one station and sensor.

    Parameters
    ----------
    station      : CDEC station code (e.g. 'SHA').
    sensor       : CDEC sensor number (e.g. 76 for inflow in CFS).
    start / end  : ISO date strings bounding the download window.
    dur          : 'D' for daily data, 'M' for monthly data.
    chunk_yearly : Split the request into per-year chunks (recommended for
                   daily data to avoid CDEC timeouts).

    Returns
    -------
    pd.DataFrame indexed by datetime with a single float column 'value'.
    Returns an empty DataFrame if CDEC returns no usable data.
    """
    parts  = []
    chunks = list(_year_chunks(start, end)) if chunk_yearly else [(start, end)]

    for cs, ce in chunks:
        try:
            r = SESSION.get(_cdec_url(station, sensor, cs, ce, dur), timeout=60)
        except requests.RequestException:
            continue
        if r.status_code != 200 or not r.text.strip():
            continue
        try:
            tmp = pd.read_csv(
                StringIO(r.text),
                usecols=["DATE TIME", "VALUE"],
                parse_dates=["DATE TIME"],
                na_values=["---", "M", ""],
            )
        except Exception:
            continue
        if tmp.empty:
            continue
        tmp = tmp.rename(columns={"DATE TIME": "datetime", "VALUE": "value"})
        tmp["value"] = pd.to_numeric(tmp["value"], errors="coerce")
        tmp = tmp.dropna(subset=["datetime"])
        parts.append(tmp)

    if not parts:
        return pd.DataFrame()

    out = (pd.concat(parts, ignore_index=True)
           .drop_duplicates(subset=["datetime"], keep="last")
           .sort_values("datetime")
           .set_index("datetime"))
    out.index.name = "datetime"
    return out


# ===========================================================================
# ==========================================================================
#
#   SECTION 5 -- UNIT CONVERSION AND WATER-YEAR HELPERS
#
# ===========================================================================
# ===========================================================================

def cfs_to_af_day(df, src_col, tgt_col):
    """
    Add a column converting CFS to AF/day in-place on a copy of df.
    The conversion factor 1.983471074 follows from: 1 CFS * 86400 s/day
    / 43560 ft^3/AF = 1.9835 AF/day.
    """
    out = df.copy()
    out[tgt_col] = pd.to_numeric(out[src_col], errors="coerce") * CFS_TO_AF_PER_DAY
    return out


def water_year_scalar(dt):
    """
    Return the California water year label for a single datetime object.
    The water year runs Oct 1 through Sep 30 and carries the label of the
    calendar year in which it ends (e.g., Oct 2024 is WY 2025).
    """
    return dt.year + 1 if dt.month >= 10 else dt.year


def water_year_array(idx):
    """Vectorised water year for a DatetimeIndex."""
    return np.where(idx.month >= 10, idx.year + 1, idx.year)


def monthly_af_to_daily_equiv(series):
    """
    Convert monthly volumes [AF/month] to a daily-equivalent rate [AF/day]
    by dividing each month's total by the number of days in that month.
    This puts monthly inflow and monthly storage on the same time basis as
    median daily outflow, enabling the supply variable in Eq. 5.
    """
    s = series.copy().dropna()
    return s / s.index.days_in_month


def wy_median_daily_inflow(monthly_af_series):
    """
    Compute water-year median daily inflow [AF/day] from a monthly [AF/month]
    series.

    Steps:
        1. Convert each monthly total to a daily equivalent (divide by days
           in month).
        2. Group by water year.
        3. Take the median across the 12 monthly daily-equivalent values.

    Using the median rather than the mean reduces the influence of extreme
    flood months (the same reasoning used throughout the paper for outflow).
    """
    daily = monthly_af_to_daily_equiv(monthly_af_series)
    df    = pd.DataFrame({"q": daily}).dropna()
    df["wy"] = df.index.map(water_year_scalar)
    return df.groupby("wy")["q"].median()


# ===========================================================================
# ==========================================================================
#
#   SECTION 6 -- CDEC DATA DOWNLOAD PIPELINE
#
# ===========================================================================
# ===========================================================================

def download_all_cdec_data(hist_start, hist_end):
    """
    Download daily inflow, daily outflow, and monthly storage from CDEC for
    all eight study reservoirs.

    Inflow and outflow are converted from CFS to AF/day.  All three data
    types are returned as dictionaries keyed by station code.  If CDEC
    returns no data for a particular reservoir/variable combination, that
    station is omitted from the dictionary (with a warning message printed).

    Returns
    -------
    inflows_daily   : {station: DataFrame(inflow_af_day)}
    outflows_daily  : {station: DataFrame(outflow_af_day)}
    storage_monthly : {station: DataFrame(storage_af)}
    """
    inflows_daily   = {}
    outflows_daily  = {}
    storage_monthly = {}

    for sta in ORDER:
        s    = RESERVOIRS[sta]
        name = s["name"]
        print(f"  [{sta}] {name} ...")

        # Daily inflow (sensor 76, CFS)
        dfi = fetch_cdec_series(sta, s["inflow"], hist_start, hist_end,
                                "D", chunk_yearly=True)
        if not dfi.empty:
            dfi = dfi.rename(columns={"value": "inflow_cfs"})
            dfi = cfs_to_af_day(dfi, "inflow_cfs", "inflow_af_day")
            inflows_daily[sta] = dfi
            print(f"    inflow  : {len(dfi):>8,} daily rows")
        else:
            print(f"    inflow  : EMPTY -- check CDEC sensor {s['inflow']}")

        # Daily outflow (sensor 23, CFS)
        dfo = fetch_cdec_series(sta, s["outflow"], hist_start, hist_end,
                                "D", chunk_yearly=True)
        if not dfo.empty:
            dfo = dfo.rename(columns={"value": "outflow_cfs"})
            dfo = cfs_to_af_day(dfo, "outflow_cfs", "outflow_af_day")
            outflows_daily[sta] = dfo
            print(f"    outflow : {len(dfo):>8,} daily rows")
        else:
            print(f"    outflow : EMPTY -- check CDEC sensor {s['outflow']}")

        # Monthly storage (sensor 15, AF)
        # Monthly data is reliable as a single request without year chunking.
        dfs = fetch_cdec_series(sta, s["storage"], hist_start, hist_end,
                                "M", chunk_yearly=False)
        if not dfs.empty:
            dfs = dfs.rename(columns={"value": "storage_af"})
            dfs["storage_af"] = pd.to_numeric(dfs["storage_af"], errors="coerce")
            storage_monthly[sta] = dfs
            print(f"    storage : {len(dfs):>8,} monthly rows")
        else:
            print(f"    storage : EMPTY -- check CDEC sensor {s['storage']}")

    return inflows_daily, outflows_daily, storage_monthly


def save_raw_workbooks(inflows_daily, outflows_daily, storage_monthly,
                       hist_start, hist_end, save_folder):
    """
    Save raw downloaded time series to Excel workbooks, one sheet per reservoir.

    These workbooks serve as the audit trail between the CDEC raw data and the
    computed statistics.  The README sheet in each workbook describes the sensor
    and units so the files are self-contained.
    """
    specs = [
        (inflows_daily,   f"INFLOWS_DAILY_{hist_start}_{hist_end}.xlsx",
         "Daily inflows (AF/day converted from CFS, sensor 76). One sheet per reservoir."),
        (outflows_daily,  f"OUTFLOWS_DAILY_{hist_start}_{hist_end}.xlsx",
         "Daily outflows (AF/day converted from CFS, sensor 23). One sheet per reservoir."),
        (storage_monthly, f"STORAGE_MONTHLY_{hist_start}_{hist_end}.xlsx",
         "Monthly storage (AF, sensor 15). September values used as carryover storage."),
    ]
    for data, filename, note in specs:
        path  = os.path.join(save_folder, filename)
        label = filename.split("_")[0].lower()
        with pd.ExcelWriter(path, engine="openpyxl") as w:
            pd.DataFrame({"note": [note]}).to_excel(w, "README", index=False)
            for sta, df in data.items():
                df.to_excel(w, sheet_name=f"{sta}_{label}"[:31])
        print(f"    Saved workbook: {filename}")


# ===========================================================================
# ==========================================================================
#
#   SECTION 7 -- STEP 1: WATER-YEAR SUMMARIES AND LINEAR FITS
#
# ===========================================================================
# ===========================================================================

def build_wy_summaries(inflows_daily, outflows_daily, storage_monthly, save_folder):
    """
    Build per-reservoir and system-level water-year summary tables.

    For each reservoir and each water year where all three variables are
    available, compute:
        median_inflow_af_day  -- WY median of daily inflow (AF/day)
        median_outflow_af_day -- WY median of daily outflow (AF/day)
        sept_carryover_af     -- September end-of-month storage (AF),
                                  shifted forward one WY to represent the
                                  carryover *into* the following water year.
        carryover_equiv_af_day -- sept_carryover_af / 365.25  (AF/day)
        X_supply_af_day       -- median_inflow + carryover_equiv  (Eq. 5 inputs)
        Y_outflow_af_day      -- median_outflow  (demand proxy in Eq. 1 & 2)

    The shift of September storage by one WY is intentional: September storage
    at end of WY t becomes the carryover available at the start of WY t+1.

    The system-level table sums X and Y across all reservoirs for each water
    year, giving the inputs to the system-level fit (Eq. 2).

    Returns
    -------
    wy_all    : per-reservoir WY summary (long format, station column)
    system_wy : system-level WY summary (one row per water year)
    """
    rows = []

    for sta in ORDER:
        if sta not in inflows_daily or sta not in outflows_daily \
                or sta not in storage_monthly:
            print(f"  [WARN] {sta}: missing data -- excluded from WY summary")
            continue

        dfi = inflows_daily[sta].copy()
        dfi["WY"] = water_year_array(dfi.index)

        dfo = outflows_daily[sta].copy()
        dfo["WY"] = water_year_array(dfo.index)

        dfs = storage_monthly[sta].copy()
        dfs["WY"] = water_year_array(dfs.index)

        wy_in  = dfi.groupby("WY")["inflow_af_day"].median().rename("median_inflow_af_day")
        wy_out = dfo.groupby("WY")["outflow_af_day"].median().rename("median_outflow_af_day")

        # Carryover: September storage of WY t -> shifted to WY t+1.
        sept   = dfs[dfs.index.month == 9]
        wy_st  = (sept.groupby("WY")["storage_af"]
                  .last()          # end-of-month September value
                  .shift(1)        # shift: Sept WY t becomes carryover into WY t+1
                  .rename("sept_carryover_storage_af"))

        wy = pd.concat([wy_in, wy_out, wy_st], axis=1).dropna().reset_index()
        if wy.empty:
            print(f"  [WARN] {sta}: no overlapping WY rows after alignment")
            continue

        wy["carryover_equiv_af_day"] = wy["sept_carryover_storage_af"] / DAYS_PER_YEAR
        wy["X_supply_af_day"]        = wy["median_inflow_af_day"] + wy["carryover_equiv_af_day"]
        wy["Y_outflow_af_day"]       = wy["median_outflow_af_day"]
        wy["station"]                = sta
        wy["reservoir"]              = RESERVOIRS[sta]["name"]
        rows.append(wy)

    if not rows:
        raise RuntimeError(
            "No water-year summaries could be built.  "
            "Check CDEC sensor IDs and network connectivity."
        )

    wy_all = pd.concat(rows, ignore_index=True)
    wy_all.to_csv(os.path.join(save_folder, "WY_reservoir_summary_all.csv"), index=False)

    # System-level: sum X and Y across all reservoirs for each WY.
    system_wy = (
        wy_all.groupby("WY")[
            ["median_inflow_af_day", "carryover_equiv_af_day", "median_outflow_af_day"]
        ]
        .sum()
        .reset_index()
        .rename(columns={
            "median_inflow_af_day":   "system_inflow_af_day",
            "carryover_equiv_af_day": "system_carryover_af_day",
            "median_outflow_af_day":  "system_outflow_af_day",
        })
    )
    system_wy["X_supply_af_day"]  = (system_wy["system_inflow_af_day"]
                                     + system_wy["system_carryover_af_day"])
    system_wy["Y_outflow_af_day"] = system_wy["system_outflow_af_day"]
    system_wy.to_csv(
        os.path.join(save_folder, "system_supply_outflow_WY_summary.csv"), index=False)

    print(f"  WY_reservoir_summary_all.csv  ({len(wy_all)} rows, "
          f"{wy_all['station'].nunique()} reservoirs)")
    print(f"  system_supply_outflow_WY_summary.csv  ({len(system_wy)} water years)")
    return wy_all, system_wy


def fit_individual_reservoirs(wy_all, save_folder):
    """
    Fit the linear model Y = c + a*X at the individual reservoir level (Eq. 2)
    and save the coefficient table.

    Returns a DataFrame with one row per reservoir containing:
        pearson_r, r2, slope_a, intercept_c, n
    """
    fit_rows = []
    for sta in ORDER:
        sub = wy_all[wy_all["station"] == sta].dropna(
            subset=["X_supply_af_day", "Y_outflow_af_day"])
        name = RESERVOIRS[sta]["name"]

        if len(sub) < 3:
            fit_rows.append({
                "station": sta, "reservoir": name, "n": len(sub),
                "pearson_r": np.nan, "r2": np.nan,
                "slope_a": np.nan, "intercept_c": np.nan,
            })
            continue

        x = sub["X_supply_af_day"].values
        y = sub["Y_outflow_af_day"].values
        a, c = np.polyfit(x, y, 1)
        r    = pd.Series(x).corr(pd.Series(y))
        fit_rows.append({
            "station": sta, "reservoir": name, "n": len(sub),
            "pearson_r": round(r, 4), "r2": round(r**2, 4),
            "slope_a": round(a, 4),   "intercept_c": round(c, 2),
        })

    fit_df = pd.DataFrame(fit_rows)
    fit_df.to_csv(
        os.path.join(save_folder, "individual_reservoir_supply_outflow_fits.csv"),
        index=False)
    return fit_df


def fit_system_level(system_wy, save_folder):
    """
    Fit the system-level linear model Y = c + a*X (Eq. 2).

    The fitted coefficients are used in Step 2 to predict future system
    outflows from synthetic supply values.

    Returns
    -------
    slope_a, intercept_c, pearson_r, r2  (all floats)
    Falls back to the paper constants if fewer than 3 water years are available.
    """
    fd = system_wy.dropna(subset=["X_supply_af_day", "Y_outflow_af_day"])

    if len(fd) < 3:
        print("  [WARN] Insufficient system data -- using paper fallback constants.")
        return STEP1_SLOPE_A, STEP1_INTERCEPT_C, STEP1_R, STEP1_R2

    x = fd["X_supply_af_day"].values
    y = fd["Y_outflow_af_day"].values
    a, c = np.polyfit(x, y, 1)
    r    = pd.Series(x).corr(pd.Series(y))
    r2   = r ** 2

    pd.DataFrame([{
        "n": len(fd), "pearson_r": r, "r2": r2,
        "slope_a": a, "intercept_c": c,
    }]).to_csv(os.path.join(save_folder, "aggregated_system_fit_summary.csv"), index=False)

    print(f"  System fit  n={len(fd)}  r={r:.3f}  R2={r2:.3f}")
    print(f"  Y = {c:.2f} + {a:.4f} * X")
    return a, c, r, r2


def compute_pearson_matrices(inflows_daily, outflows_daily, storage_monthly, save_folder):
    """
    Compute pairwise Pearson correlation matrices (Eq. 1) for daily inflow,
    daily outflow, and September carryover storage across all eight reservoirs.

    Each matrix is 8 x 8.  Correlations are computed using only the overlapping
    dates between each pair of stations (pairwise complete cases), which handles
    the different record lengths gracefully.

    Saves three CSV files and returns the three correlation DataFrames.
    """
    # Build wide DataFrames where each column is one reservoir's time series.
    in_wide  = pd.DataFrame({
        s: inflows_daily[s]["inflow_af_day"]
        for s in ORDER if s in inflows_daily
    })
    out_wide = pd.DataFrame({
        s: outflows_daily[s]["outflow_af_day"]
        for s in ORDER if s in outflows_daily
    })

    # For carryover, index by calendar year of the September observation
    # (one value per year per reservoir).
    car_wide = pd.DataFrame()
    for s in ORDER:
        if s not in storage_monthly:
            continue
        sept = storage_monthly[s][storage_monthly[s].index.month == 9]
        if sept.empty:
            continue
        car_wide[s] = pd.Series(sept["storage_af"].values,
                                index=sept.index.year, name=s)

    def pairwise_r(wide, cols):
        cols = [c for c in cols if c in wide.columns]
        M = pd.DataFrame(index=cols, columns=cols, dtype=float)
        for c1 in cols:
            for c2 in cols:
                pair = pd.concat([wide[c1], wide[c2]], axis=1,
                                 keys=["x", "y"]).dropna()
                M.loc[c1, c2] = pair["x"].corr(pair["y"]) if len(pair) >= 3 else np.nan
        return M

    in_corr  = pairwise_r(in_wide,  ORDER)
    out_corr = pairwise_r(out_wide, ORDER)
    car_corr = pairwise_r(car_wide, ORDER)

    in_corr.to_csv(os.path.join(save_folder, "pearson_inflow_correlation_matrix.csv"))
    out_corr.to_csv(os.path.join(save_folder, "pearson_outflow_correlation_matrix.csv"))
    car_corr.to_csv(os.path.join(save_folder, "pearson_carryover_correlation_matrix.csv"))

    return in_corr, out_corr, car_corr


# ===========================================================================
# ==========================================================================
#
#   SECTION 8 -- STEP 2: THOMAS-FIERING MODEL
#
# ===========================================================================
# ===========================================================================

def _daily_to_monthly_af(df_daily):
    """
    Convert a daily inflow DataFrame (CFS already in af_day column) to
    monthly totals in AF/month by summing within each calendar month.
    min_count=1 ensures months with at least one valid observation produce
    a result rather than NaN.
    """
    monthly = df_daily["inflow_af_day"].resample("MS").sum(min_count=1)
    monthly.name = "monthly_af"
    return monthly


def _prepare_tf_frame(monthly_series):
    """
    Prepare a monthly inflow pd.Series for TF parameter estimation.

    Adds water-year (wy) and WY-month-slot (mi, 0=Oct...11=Sep) columns.
    Values below MIN_POSITIVE_FLOW are clipped to that floor before the log
    transform so the lognormal recursion (Eq. 3) remains numerically stable.
    The original values are retained in 'q_raw' for KS-test validation.
    """
    df = pd.DataFrame({"q": monthly_series.copy()}).dropna()
    df["wy"]    = df.index.map(water_year_scalar)
    df["mi"]    = df.index.month.map(M_TO_SLOT)
    df["q_raw"] = df["q"]
    df["q"]     = df["q"].clip(lower=MIN_POSITIVE_FLOW)
    return df


def _fit_tf_parameters(df_m):
    """
    Estimate lognormal Thomas-Fiering parameters from incomplete monthly data.

    For each of the 12 WY month slots the function computes:
        mu[i]    -- mean of log(q) across all available observations
        sigma[i] -- std dev of log(q) across all available observations
        rho[i]   -- lag-1 autocorrelation between month slot i and slot i-1

    Handles missing months gracefully: if fewer than three paired observations
    are available for computing rho, that month's rho is set to zero (the model
    then draws from the marginal monthly distribution for that transition).

    Any parameter that cannot be estimated (e.g., only one observation in a
    month) is filled with the cross-month mean to allow the simulation to proceed.

    Returns
    -------
    mu, sigma, rho : each a length-12 numpy array
    counts         : length-12 int array of observation counts per slot
    """
    df = df_m.copy()
    df["logq"] = np.log(df["q"])

    mu     = np.full(12, np.nan)
    sigma  = np.full(12, np.nan)
    rho    = np.zeros(12)
    counts = np.zeros(12, dtype=int)

    for mi in range(12):
        vals = df.loc[df["mi"] == mi, "logq"].dropna().values
        counts[mi] = len(vals)
        if len(vals) >= 1:
            mu[mi] = np.mean(vals)
        sigma[mi] = np.std(vals, ddof=1) if len(vals) >= 2 else 0.0

    # Cross-month mean fill for any month with no observations.
    mu    = np.where(np.isfinite(mu),    mu,    np.nanmean(mu))
    sigma = np.where(np.isfinite(sigma), sigma, np.nanmean(sigma))
    sigma = np.clip(sigma, 0.0, None)

    # Lag-1 autocorrelation.
    # For the October->September transition (mi=0 follows mi=11 from the
    # prior WY), we adjust the water-year label of September records.
    for mi in range(12):
        prev = (mi - 1) % 12
        cur  = df.loc[df["mi"] == mi,   ["wy", "logq"]].rename(columns={"logq": "c"})
        prv  = df.loc[df["mi"] == prev,  ["wy", "logq"]].rename(columns={"logq": "p"})
        if mi == 0:
            prv = prv.copy()
            prv["wy"] = prv["wy"] + 1
        pair = pd.merge(cur, prv, on="wy", how="inner").dropna()
        if len(pair) >= 3:
            rho[mi] = np.corrcoef(pair["c"], pair["p"])[0, 1]

    rho = np.clip(rho, -0.99, 0.99)
    return mu, sigma, rho, counts


def _generate_tf_series(mu, sigma, rho, month_index, seed=42):
    """
    Generate a synthetic monthly inflow series using the lognormal TF
    recursion (Eq. 3).

    The recursion propagates in log-space so that back-transforming via
    exp() always yields positive inflows.  The series is generated to cover
    exactly the months in `month_index`, preserving the calendar structure.

    Parameters
    ----------
    mu, sigma, rho : length-12 parameter arrays from _fit_tf_parameters.
    month_index    : pd.DatetimeIndex of monthly start-of-month timestamps.
    seed           : integer random seed for reproducibility.

    Returns
    -------
    pd.Series of synthetic AF/month values indexed by month_index.
    """
    rng = np.random.default_rng(seed)
    N   = len(month_index)
    X   = np.zeros(N)

    first_mi = M_TO_SLOT[month_index[0].month]
    X[0] = rng.normal(mu[first_mi],
                      sigma[first_mi] if sigma[first_mi] > 0 else 0.0)

    for t in range(1, N):
        mi   = M_TO_SLOT[month_index[t].month]
        prev = (mi - 1) % 12
        r    = rho[mi]
        # Innovation term: the conditional standard deviation scaled by a
        # standard normal draw (Eq. 3 last term).
        eps  = (sigma[mi] * np.sqrt(max(0.0, 1.0 - r**2)) * rng.standard_normal()
                if sigma[mi] > 0 else 0.0)
        X[t] = mu[mi] + r * (X[t - 1] - mu[prev]) + eps

    return pd.Series(np.exp(X), index=month_index, name="synthetic_af_month")


def run_thomas_fiering(inflows_daily, save_folder):
    """
    Fit TF parameters from historical daily inflows (aggregated to monthly)
    and generate synthetic inflow series for each reservoir.

    The synthetic series is generated over N_SYN_YEARS starting from the
    beginning of the historical record, long enough to extend through 2100.
    The water-year median daily inflow [AF/day] is then computed from both
    the historical and synthetic series for comparison (Figure 8) and for
    use in the scenario calculations (Step 2).

    Validation: A two-sample KS test compares the historical and synthetic
    distributions of monthly daily-equivalent inflows.  Reservoirs with
    p < 0.05 fail the test; results should be interpreted with caution.
    Reservoirs with fewer than 60 observed months are flagged 'low_confidence'.

    Saves per-reservoir synthetic CSVs and a summary CSV.

    Returns
    -------
    synthetic_dict : dict with per-reservoir results (see inline comments)
    summary_df     : DataFrame with TF diagnostics, one row per reservoir
    """
    synthetic_dict = {}
    summary_rows   = []

    for i, sta in enumerate(ORDER):
        if sta not in inflows_daily:
            print(f"  [{sta}] skipped (no inflow data)")
            continue

        # Aggregate daily AF/day -> monthly AF/month.
        monthly = _daily_to_monthly_af(inflows_daily[sta]).dropna()

        if len(monthly) < 12:
            print(f"  [{sta}] skipped (< 12 months of data)")
            continue

        # Reindex over the full historical span so we have a continuous
        # month_index even if some months have missing values.
        month_index  = pd.date_range(monthly.index.min(), monthly.index.max(), freq="MS")
        monthly_full = monthly.reindex(month_index)

        df_tf = _prepare_tf_frame(monthly_full)
        if len(df_tf) < 12:
            print(f"  [{sta}] skipped (too few usable months)")
            continue

        mu, sigma, rho, counts = _fit_tf_parameters(df_tf)

        # Generate N_SYN_YEARS starting from the historical start date so
        # the water-year labels are real calendar years up to 2100+.
        syn_index  = pd.date_range(
            month_index[0],
            periods=12 * N_SYN_YEARS,
            freq="MS",
        )
        syn_monthly = _generate_tf_series(mu, sigma, rho, syn_index,
                                          seed=RANDOM_SEED + i)

        # KS test (Kolmogorov-Smirnov): compare daily-equivalent AF/day
        # distributions for the overlapping historical period only.
        hist_daily = monthly_af_to_daily_equiv(df_tf.set_index(df_tf.index)["q_raw"])
        syn_daily  = monthly_af_to_daily_equiv(
            syn_monthly[syn_monthly.index.isin(month_index)]
        )
        _, p_ks = ks_2samp(hist_daily.values.astype(float),
                           syn_daily.values.astype(float))

        synthetic_dict[sta] = {
            "monthly_full":        monthly_full,       # historical monthly AF/month
            "synthetic_monthly":   syn_monthly,        # synthetic monthly AF/month (full span)
            "hist_wy_median":      wy_median_daily_inflow(monthly_full),
            "syn_wy_median":       wy_median_daily_inflow(syn_monthly),
            "mu": mu, "sigma": sigma, "rho": rho,
        }

        # Save per-reservoir synthetic series for audit and reproducibility.
        syn_out = pd.DataFrame({
            "datetime":                  syn_monthly.index,
            "synthetic_inflow_af_month": syn_monthly.values,
            "station":                   sta,
            "reservoir":                 RESERVOIRS[sta]["name"],
        })
        syn_out.to_csv(
            os.path.join(save_folder,
                         f"{sta}_thomas_fiering_synthetic_monthly_inflow.csv"),
            index=False,
        )

        low_conf = len(df_tf) < 60
        summary_rows.append({
            "station":          sta,
            "reservoir":        RESERVOIRS[sta]["name"],
            "hist_start":       str(monthly.index.min().date()),
            "hist_end":         str(monthly.index.max().date()),
            "observed_months":  len(df_tf),
            "ks_p_value":       round(p_ks, 4),
            "low_confidence":   "YES" if low_conf else "NO",
            "min_slot_count":   int(counts.min()),
            "max_slot_count":   int(counts.max()),
        })
        flag = " *** LOW CONFIDENCE (short record)" if low_conf else ""
        print(f"  [{sta}] done  obs_months={len(df_tf)}  KS_p={p_ks:.4f}{flag}")

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(
        os.path.join(save_folder,
                     "thomas_fiering_synthetic_summary_all_reservoirs.csv"),
        index=False,
    )

    # Save combined historical + synthetic WY-median table.
    wy_rows = []
    for sta, d in synthetic_dict.items():
        for stype, series in [("historical", d["hist_wy_median"]),
                               ("synthetic",  d["syn_wy_median"])]:
            tmp = series.reset_index()
            tmp.columns = ["WY", "median_daily_inflow_af_day"]
            tmp["series_type"] = stype
            tmp["station"]     = sta
            tmp["reservoir"]   = RESERVOIRS[sta]["name"]
            wy_rows.append(tmp)
    if wy_rows:
        (pd.concat(wy_rows, ignore_index=True)
         .to_csv(
             os.path.join(save_folder,
                          "thomas_fiering_WY_daily_medians_historical_and_synthetic.csv"),
             index=False,
         ))

    return synthetic_dict, summary_df


# ===========================================================================
# ==========================================================================
#
#   SECTION 9 -- STEP 2: SCENARIO CALCULATIONS
#
# ===========================================================================
# ===========================================================================

def download_carryover_storage(hist_start, hist_end):
    """
    Download monthly storage from CDEC for all reservoirs and compute the
    historical median September carryover storage.

    September storage is used as the carryover variable (C_r in Eq. 5) because
    September marks the end of the dry season and represents the water a reservoir
    carries into the next water year.  The median is taken across all available
    September observations to produce a single representative carryover value
    per reservoir that is treated as fixed under all future scenarios.

    Returns a DataFrame with:
        station, reservoir,
        median_sept_storage_af   [AF]
        carryover_equiv_af_day   [AF/day] = median_sept_storage_af / 365.25
    """
    rows = []
    for sta in ORDER:
        sensor = RESERVOIRS[sta]["storage"]
        df = fetch_cdec_series(sta, sensor, hist_start, hist_end,
                               "M", chunk_yearly=False)
        if df.empty:
            print(f"  [{sta}] no storage data for carryover")
            continue
        df = df.rename(columns={"value": "storage_af"})
        df["storage_af"] = pd.to_numeric(df["storage_af"], errors="coerce")
        sept = df[df.index.month == 9].dropna(subset=["storage_af"])
        if sept.empty:
            print(f"  [{sta}] no September storage records")
            continue
        med = sept["storage_af"].median()
        rows.append({
            "station":                  sta,
            "reservoir":                RESERVOIRS[sta]["name"],
            "median_sept_storage_af":   med,
            "carryover_equiv_af_day":   med / DAYS_PER_YEAR,
        })
        print(f"  [{sta}] Sept records={len(sept)}  "
              f"median={med:,.0f} AF  "
              f"equiv={med/DAYS_PER_YEAR:,.1f} AF/day")

    if not rows:
        raise RuntimeError("No carryover storage computed for any reservoir.")
    return pd.DataFrame(rows)


def build_future_scenarios(synthetic_dict, carryover_df, slope_a, intercept_c,
                           save_folder):
    """
    Apply inflow-reduction scenarios (Eq. 4), assemble the aggregated system
    supply variable (Eq. 5), and predict future system outflows (Eq. 2).

    The future period is FUTURE_START_WY through FUTURE_END_WY.  Only water
    years where all eight reservoirs have valid synthetic inflow values are
    included in the system aggregation (incomplete years are dropped).

    Carryover storage is held at the historical median for all future years
    (see Eq. 5 assumption in the paper header docstring above).

    Returns
    -------
    system_future : DataFrame with predicted outflow per WY and scenario
    summary       : descriptive statistics per scenario
    """
    # Collect future WY-median synthetic inflows for all reservoirs.
    wy_rows = []
    for sta, d in synthetic_dict.items():
        wy_med = d["syn_wy_median"]
        future_slice = wy_med[
            (wy_med.index >= FUTURE_START_WY) &
            (wy_med.index <= FUTURE_END_WY)
        ].reset_index()
        future_slice.columns = ["WY", "median_inflow_af_day"]
        future_slice["station"]   = sta
        future_slice["reservoir"] = RESERVOIRS[sta]["name"]
        wy_rows.append(future_slice)

    if not wy_rows:
        raise RuntimeError("No future synthetic inflow data found.  "
                           "Increase N_SYN_YEARS or check date ranges.")

    syn_wy = pd.concat(wy_rows, ignore_index=True)

    # Drop water years missing any reservoir (to keep the system sum consistent).
    wy_counts = syn_wy.groupby("WY")["station"].nunique()
    valid_wy  = wy_counts[wy_counts == len(ORDER)].index
    syn_wy    = syn_wy[syn_wy["WY"].isin(valid_wy)].copy()

    syn_wy.to_csv(
        os.path.join(save_folder,
                     "future_synthetic_WY_daily_median_inflows_all_reservoirs.csv"),
        index=False,
    )

    # Merge historical carryover onto the inflow table.
    merged = syn_wy.merge(
        carryover_df[["station", "carryover_equiv_af_day"]],
        on="station", how="left",
    ).dropna(subset=["carryover_equiv_af_day"])

    # Apply Eq. 4 and Eq. 5 for each scenario.
    res_rows = []
    for scen, delta in SCENARIOS.items():
        df_s = merged.copy()
        df_s["scenario"]           = scen
        df_s["reduction_fraction"] = delta
        # Eq. 4: I_scenario = (1 - delta) * I_t
        df_s["scenario_inflow_af_day"]  = df_s["median_inflow_af_day"] * (1.0 - delta)
        # Eq. 5 (per reservoir): supply = I_scenario + C
        df_s["scenario_supply_af_day"]  = (df_s["scenario_inflow_af_day"]
                                           + df_s["carryover_equiv_af_day"])
        res_rows.append(df_s)

    res_all = pd.concat(res_rows, ignore_index=True)
    res_all.to_csv(
        os.path.join(save_folder, "future_reservoir_supply_outflow_scenarios.csv"),
        index=False,
    )

    # Aggregate to system level (Eq. 5 system sum).
    system_future = (
        res_all
        .groupby(["WY", "scenario", "reduction_fraction"], as_index=False)[
            ["scenario_inflow_af_day", "carryover_equiv_af_day",
             "scenario_supply_af_day"]
        ]
        .sum()
    )

    # Eq. 2: Y = c + a * X_t
    system_future["predicted_outflow_af_day"] = (
        intercept_c + slope_a * system_future["scenario_supply_af_day"]
    ).clip(lower=0)   # physical lower bound: outflow cannot be negative

    system_future.to_csv(
        os.path.join(save_folder, "future_system_outflows_scenarios_2026_2100.csv"),
        index=False,
    )

    summary = (
        system_future.groupby("scenario")["predicted_outflow_af_day"]
        .agg(["mean", "std", "min", "max", "median"])
        .reindex(list(SCENARIOS.keys()))
        .reset_index()
    )
    summary.to_csv(
        os.path.join(save_folder, "future_outflow_scenario_summary.csv"),
        index=False,
    )

    return system_future, summary


# ===========================================================================
# ==========================================================================
#
#   SECTION 10 -- FIGURES
#   One function per figure.  All figures follow the same style block
#   defined in Section 3 and are saved at 300 dpi PNG.
#   Figure 1 is additionally saved as 600 dpi PDF.
#
# ===========================================================================
# ===========================================================================

def _save(fig, name, save_folder, dpi=300, also_pdf=False):
    """Save a figure to PNG (and optionally PDF) and close it."""
    fig.savefig(os.path.join(save_folder, name), dpi=dpi, bbox_inches="tight")
    if also_pdf:
        fig.savefig(os.path.join(save_folder, name.replace(".png", ".pdf")),
                    dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {name}")


def _heatmap(ax, matrix, title):
    """Draw an annotated Pearson correlation heatmap on a given Axes."""
    arr = matrix.astype(float).values
    im  = ax.imshow(arr, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    n   = len(matrix.columns)
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right")
    ax.set_yticklabels(matrix.index)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xticks(np.arange(-0.5, n, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1)
    ax.tick_params(which="minor", bottom=False, left=False)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            v = arr[i, j]
            c = "white" if pd.notna(v) and abs(v) > 0.55 else "black"
            ax.text(j, i, f"{v:.2f}" if pd.notna(v) else "NA",
                    ha="center", va="center", color=c, fontsize=8)
    return im


def plot_figure01_swp_allocation(save_folder):
    """
    Figure 1 -- SWP Table A, requested, and approved allocations (1996-2025).

    Panel (a): Three time series in MAF -- contractual Table A entitlements,
    initial contractor requests, and actual approved allocations -- together
    with the 3-year centred moving average of approved allocations.  Long-run
    averages for each series are shown as dotted reference lines.  Shaded
    bands highlight the four major drought periods (2001, 2009, 2014, 2021)
    that drove the most severe allocation shortfalls.

    Panel (b): Approved allocation expressed as a percentage of the initial
    request.  The 3-year moving average (dashed) and long-run mean of 40.9%
    (dotted) are overlaid.

    Since 1996 approved aggregate allocations to Southern California
    contractors have declined by approximately 70%.
    """
    df = pd.DataFrame(SWP_DATA)
    df["Pct"]          = df["Approved_Allocation"] / df["Initial_Request"] * 100
    df["TableA_MAF"]   = df["Table_A"]             / 1e6
    df["Request_MAF"]  = df["Initial_Request"]     / 1e6
    df["Approved_MAF"] = df["Approved_Allocation"] / 1e6
    df["App_MAF_r3"]   = df["Approved_MAF"].rolling(3, center=True).mean()
    df["Pct_r3"]       = df["Pct"].rolling(3, center=True).mean()

    fig, axes = plt.subplots(2, 1, figsize=(11, 7.5), sharex=True,
                              gridspec_kw={"height_ratios": [2.2, 1.2]})

    ax = axes[0]
    ax.plot(df["Year"], df["TableA_MAF"],   color=C_BLUE,   lw=2.2, label="Table A entitlement")
    ax.plot(df["Year"], df["Request_MAF"],  color=C_GREEN,  lw=2.2, label="Initial request")
    ax.plot(df["Year"], df["Approved_MAF"], color=C_RED,    lw=2.5, label="Approved allocation")
    ax.plot(df["Year"], df["App_MAF_r3"],   color=C_PURPLE, lw=2.0, ls="--",
            label="Approved allocation (3-yr mean)")
    for val, col in zip([df["TableA_MAF"].mean(), df["Request_MAF"].mean(),
                         df["Approved_MAF"].mean()], [C_BLUE, C_GREEN, C_RED]):
        ax.axhline(val, color=col, ls=":", lw=1.2, alpha=0.75)
    for yr in [2001, 2009, 2014, 2021]:
        ax.axvspan(yr - 0.5, yr + 0.5, color=C_SHADE, alpha=0.35, zorder=0)
    ax.text(2014, 0.15, "Zero allocation", fontsize=9, color=C_GREY, ha="center")
    ax.set_ylabel("Water volume (MAF)")
    ax.set_title(
        "SWP Table A, Requested, and Approved Allocations to Southern California",
        pad=16)
    ax.grid(True, axis="y", ls="--", lw=0.5, alpha=0.4)
    ax.legend(loc="upper right", frameon=False, ncol=2)
    ax.text(0.01, 0.96, "(a)", transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="top")

    avg_pct = df["Pct"].mean()
    ax2 = axes[1]
    ax2.plot(df["Year"], df["Pct"],    color=C_RED,    lw=2.3,
             label="Approved allocation (%)")
    ax2.plot(df["Year"], df["Pct_r3"], color=C_PURPLE, lw=2.0, ls="--",
             label="Approved allocation (%) (3-yr mean)")
    ax2.axhline(avg_pct, color=C_GREY, ls=":", lw=1.3,
                label=f"Long-run average = {avg_pct:.1f}%")
    for yr in [2001, 2009, 2014, 2021]:
        ax2.axvspan(yr - 0.5, yr + 0.5, color=C_SHADE, alpha=0.35, zorder=0)
    ax2.set_ylabel("Approved allocation (%)"); ax2.set_xlabel("Year")
    ax2.set_ylim(0, 110)
    ax2.grid(True, axis="y", ls="--", lw=0.5, alpha=0.4)
    ax2.legend(loc="upper right", frameon=False)
    ax2.text(0.01, 0.96, "(b)", transform=ax2.transAxes,
             fontsize=12, fontweight="bold", va="top")
    ax2.set_xticks(df["Year"][::2])
    ax2.set_xticklabels(df["Year"][::2], rotation=45)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    _save(fig, "figure01_SWP_allocation.png", save_folder, also_pdf=True)


def _eight_panel(data_dict, value_col, median_label, ylabel, suptitle, filename,
                 save_folder, sharex=True):
    """
    Generic 4x2 panel figure helper used for Figures 2, 3, and 4.

    For each reservoir panel the raw time series is plotted in blue and the
    overall median is shown as a dashed red horizontal line.
    """
    fig, axes = plt.subplots(4, 2, figsize=(14, 16), sharex=sharex)
    axes = axes.flatten()
    for i, sta in enumerate(ORDER):
        ax = axes[i]
        if sta not in data_dict:
            ax.text(0.5, 0.5, f"{sta}\nNo data",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off(); continue
        df  = data_dict[sta]
        med = df[value_col].median()
        ax.plot(df.index, df[value_col], lw=0.5, color=C_BLUE)
        ax.axhline(med, ls="--", lw=1.8, color=C_RED,
                   label=f"{median_label} = {med:,.0f}")
        ax.set_title(f"{sta} ({RESERVOIRS[sta]['name']})",
                     fontsize=11, fontweight="bold")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.2)
        ax.legend(frameon=False, fontsize=8, loc="upper right")
    fig.suptitle(suptitle, fontsize=14, fontweight="bold")
    fig.supxlabel("Date")
    plt.tight_layout(rect=[0, 0.02, 1, 0.97])
    _save(fig, filename, save_folder)


def plot_figure02_daily_inflows(inflows_daily, save_folder):
    """
    Figure 2 -- Daily inflows for all 8 reservoirs (4x2 panel).
    Shows the complete daily record with the long-run median inflow [AF/day]
    marked as a dashed line for reference.
    """
    _eight_panel(
        inflows_daily, "inflow_af_day",
        "Median",
        "Inflow (AF/day)",
        "Daily Inflows for Major Study Reservoirs with Median Inflow Highlighted",
        "figure02_daily_inflows_8panel.png",
        save_folder,
    )


def plot_figure03_daily_outflows(outflows_daily, save_folder):
    """
    Figure 3 -- Daily outflows for all 8 reservoirs (4x2 panel).
    Median outflow is used as the demand proxy Y in Eq. 1 and Eq. 2.
    Note the abbreviated and/or noisy records for SNL and BER, which are
    discussed in the paper as a limitation.
    """
    _eight_panel(
        outflows_daily, "outflow_af_day",
        "Median",
        "Outflow (AF/day)",
        "Daily Outflows for Major Study Reservoirs with Median Outflow Highlighted",
        "figure03_daily_outflows_8panel.png",
        save_folder,
    )


def plot_figure04_storage(storage_monthly, save_folder):
    """
    Figure 4 -- Monthly storage for all 8 reservoirs (4x2 panel).
    The dashed line marks the median September storage, which serves as the
    carryover variable C_t in Eq. 5.
    """
    # Build a version of storage_monthly that only holds September values
    # for the median label, but plots the full monthly series.
    sept_medians = {}
    for sta in ORDER:
        if sta not in storage_monthly:
            continue
        df   = storage_monthly[sta].copy()
        sept = df[df.index.month == 9]
        sept_medians[sta] = {
            "storage_af":      df["storage_af"],
            "_sept_med":       sept["storage_af"].median() if not sept.empty else np.nan,
        }

    fig, axes = plt.subplots(4, 2, figsize=(14, 16), sharex=True)
    axes = axes.flatten()
    for i, sta in enumerate(ORDER):
        ax = axes[i]
        if sta not in storage_monthly:
            ax.text(0.5, 0.5, f"{sta}\nNo data",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off(); continue
        df  = storage_monthly[sta]
        med = sept_medians[sta]["_sept_med"]
        ax.plot(df.index, df["storage_af"], lw=0.9, color=C_BLUE,
                label="Monthly storage")
        if pd.notna(med):
            ax.axhline(med, ls="--", lw=1.8, color=C_RED,
                       label=f"Median Sept = {med:,.0f} AF")
        ax.set_title(f"{sta} ({RESERVOIRS[sta]['name']})",
                     fontsize=11, fontweight="bold")
        ax.set_ylabel("Storage (AF)")
        ax.grid(True, alpha=0.2)
        ax.legend(frameon=False, fontsize=8, loc="upper right")
    fig.suptitle(
        "Reservoir Storage with September Carryover Storage Median Highlighted",
        fontsize=14, fontweight="bold")
    fig.supxlabel("Date")
    plt.tight_layout(rect=[0, 0.02, 1, 0.97])
    _save(fig, "figure04_storage_8panel.png", save_folder)


def plot_figure05_correlation_matrices(in_corr, out_corr, car_corr, save_folder):
    """
    Figure 5 -- Pearson correlation matrices (1x3 panel, Figure 6 in paper).

    Panel (a): Inflow correlations -- SHA, ORO, FOL, NML, DNP cluster together
    (high positive r), while SNL shows negative correlations because it is an
    off-stream reservoir fed by Delta pumping rather than direct watershed runoff.

    Panel (b): Outflow correlations -- similar pattern but weaker, reflecting
    operational decisions layered on top of hydrologic variability.

    Panel (c): September carryover storage correlations -- generally high
    positive correlations because storage tracks multi-year wet/dry cycles
    that affect all northern Sierra reservoirs simultaneously.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    for ax, mat, title in zip(
        axes,
        [in_corr,  out_corr,  car_corr],
        ["(a) Inflow correlation",
         "(b) Outflow correlation",
         "(c) Carryover storage correlation"],
    ):
        _heatmap(ax, mat, title)
    im3 = ax.images[0]  # re-use last image for shared colorbar

    fig.suptitle(
        "Figure 5 (Paper Fig. 6): Pearson Correlation Matrices for Inflow, Outflow, "
        "and September Carryover Storage\n"
        "SHA=Shasta  CLE=Trinity  ORO=Oroville  FOL=Folsom  "
        "BER=Berryessa  NML=New Melones  DNP=Don Pedro  SNL=San Luis",
        fontsize=11, fontweight="bold")
    plt.tight_layout(rect=[0, 0.08, 1, 0.92])
    cax = fig.add_axes([0.25, 0.03, 0.50, 0.025])
    fig.colorbar(im3, cax=cax, orientation="horizontal").set_label("Pearson r")
    _save(fig, "figure05_pearson_correlation_matrices.png", save_folder)


def plot_figure06_individual_fits(wy_all, save_folder):
    """
    Figure 6 -- Individual reservoir linear fits (4x2 panel, Figure 7 in paper).

    Each sub-panel shows the scatter of water-year median outflow vs supply
    (median inflow + carryover equivalent) with the ordinary least-squares fit.
    Annotated with Pearson r, R², sample size n, and the regression equation
    y = intercept + slope * x.  Reservoirs with fewer than 3 observations are
    shown as 'Insufficient data'.
    """
    fig, axes = plt.subplots(4, 2, figsize=(14, 16))
    axes = axes.flatten()
    for i, sta in enumerate(ORDER):
        ax   = axes[i]
        name = RESERVOIRS[sta]["name"]
        sub  = wy_all[wy_all["station"] == sta].dropna(
            subset=["X_supply_af_day", "Y_outflow_af_day"])
        if len(sub) < 3:
            ax.text(0.5, 0.5, f"{sta}\nInsufficient data",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off(); continue
        x = sub["X_supply_af_day"].values
        y = sub["Y_outflow_af_day"].values
        a, c = np.polyfit(x, y, 1)
        r    = pd.Series(x).corr(pd.Series(y))
        ax.scatter(x, y, s=32, color=C_BLUE, zorder=3)
        xl = np.linspace(x.min(), x.max(), 100)
        ax.plot(xl, a * xl + c, lw=1.8, color=C_RED)
        ax.set_title(f"{sta} ({name})", fontsize=11, fontweight="bold")
        ax.set_xlabel("Median inflow + carryover (AF/day)")
        ax.set_ylabel("Median outflow (AF/day)")
        ax.grid(True, alpha=0.25)
        ax.text(0.03, 0.97,
                f"r = {r:.2f}   R\u00b2 = {r**2:.2f}   n = {len(sub)}\n"
                f"y = {c:.1f} + {a:.3f}x",
                transform=ax.transAxes, va="top", ha="left", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.7", alpha=0.9))
    fig.suptitle(
        "Figure 6 (Paper Fig. 7): Individual Reservoir Linear Relationships "
        "Between Median Supply and Median Outflow",
        fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    _save(fig, "figure06_individual_reservoir_fits.png", save_folder)


def plot_figure07_system_fit(system_wy, slope_a, intercept_c, r, r2, save_folder):
    """
    Figure 7 -- System-level linear fit (Figure 8 in paper).

    Each point is one historical water year.  X is the sum of all reservoir
    supply indices (inflow + carryover equivalent, AF/day); Y is the sum of
    all reservoir median outflows (AF/day).  The strong correlation (r = 0.961)
    supports using this relationship to predict future system outflows from
    synthetic supply values in Step 2.
    """
    fd = system_wy.dropna(subset=["X_supply_af_day", "Y_outflow_af_day"])
    x  = fd["X_supply_af_day"].values
    y  = fd["Y_outflow_af_day"].values

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(x, y, s=45, color=C_BLUE, zorder=3, label="Water years")
    xl = np.linspace(0, x.max() * 1.05, 200)
    ax.plot(xl, intercept_c + slope_a * xl, lw=2.2, color=C_RED, label="Linear fit")
    ax.set_xlabel("Aggregated median inflow + carryover storage (AF/day)")
    ax.set_ylabel("Aggregated median outflow (AF/day)")
    ax.set_title(
        "Figure 7 (Paper Fig. 8): System-Level Linear Relationship\n"
        "Aggregated Reservoir Supply vs Aggregated Outflows",
        fontweight="bold")
    ax.set_xlim(left=0)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    ax.text(0.04, 0.96,
            f"r = {r:.3f}   R\u00b2 = {r2:.3f}   n = {len(fd)}\n"
            f"Y = {intercept_c:.2f} + {slope_a:.4f} X",
            transform=ax.transAxes, va="top", ha="left", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", alpha=0.9))
    plt.tight_layout()
    _save(fig, "figure07_system_level_fit.png", save_folder)


def plot_figure08_tf_wy_medians(synthetic_dict, save_folder):
    """
    Figure 8 -- Comparison of historical and TF synthetic water-year median
    daily inflows (Figure 9 in paper).

    For each reservoir the historical WY-median daily inflow [AF/day] (blue,
    filled circles) is plotted against the synthetic WY-median daily inflow
    (orange, filled circles) over the same calendar years.  Good visual
    agreement indicates the TF model reproduces the magnitude and interannual
    variability of the historical record.

    Note the larger discrepancies for BER and SNL (short, incomplete records).
    """
    fig, axes = plt.subplots(4, 2, figsize=(16, 14))
    axes = axes.flatten()
    for ax, sta in zip(axes, ORDER):
        if sta not in synthetic_dict:
            ax.text(0.5, 0.5, f"{sta}\nNo data",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off(); continue
        d     = synthetic_dict[sta]
        h_wy  = d["hist_wy_median"]
        s_wy  = d["syn_wy_median"]
        # Only show the synthetic series over the historical period for comparison.
        s_wy_hist = s_wy[s_wy.index.isin(h_wy.index)]
        ax.plot(h_wy.index, h_wy.values, marker="o", ms=4, lw=1.4,
                color=C_BLUE, label="Historical WY median")
        ax.plot(s_wy_hist.index, s_wy_hist.values, marker="o", ms=4, lw=1.4,
                color=C_RED, label="Synthetic WY median")
        ax.set_title(f"{sta} ({RESERVOIRS[sta]['name']})",
                     fontsize=10, fontweight="bold")
        ax.set_ylabel("WY median daily inflow (AF/day)")
        ax.grid(True, alpha=0.25)
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle(
        "Figure 8 (Paper Fig. 9): Historical vs Thomas-Fiering Synthetic "
        "Water-Year Daily Median Inflows",
        fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    _save(fig, "figure08_TF_historical_vs_synthetic_WY_medians.png", save_folder)


def plot_figure09_tf_monthly_means(synthetic_dict, save_folder):
    """
    Figure 9 -- Thomas-Fiering monthly mean validation (4x2 panel).

    For each reservoir the long-run monthly mean inflow [AF/month] is computed
    from the historical monthly series (blue) and from the synthetic monthly
    series (orange dashed) and plotted over the Oct-Sep water-year cycle.

    This figure validates that the TF model preserves the seasonal pattern
    (the shape of the hydrograph) in addition to the interannual statistics
    checked by the KS test.
    """
    fig, axes = plt.subplots(4, 2, figsize=(14, 14))
    axes = axes.flatten()
    x = np.arange(12)

    for ax, sta in zip(axes, ORDER):
        if sta not in synthetic_dict:
            ax.text(0.5, 0.5, f"{sta}\nNo data",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off(); continue

        d     = synthetic_dict[sta]
        hist_m = d["monthly_full"].dropna()
        syn_m  = d["synthetic_monthly"]

        # Compute monthly means from both series.
        def slot_means(series):
            df = pd.DataFrame({"q": series.dropna()})
            df["slot"] = df.index.month.map(M_TO_SLOT)
            return df.groupby("slot")["q"].mean().reindex(range(12)).values

        hist_means = slot_means(hist_m)
        syn_means  = slot_means(syn_m)

        ax.plot(x, hist_means, marker="o", ms=4, lw=1.6,
                color=C_BLUE, label="Historical monthly mean")
        ax.plot(x, syn_means,  marker="s", ms=4, lw=1.6, ls="--",
                color=C_RED, label="Synthetic monthly mean")
        ax.set_xticks(x)
        ax.set_xticklabels(WY_LABELS, rotation=45, fontsize=8)
        ax.set_title(f"{sta} ({RESERVOIRS[sta]['name']})",
                     fontsize=10, fontweight="bold")
        ax.set_ylabel("Mean monthly inflow (AF/month)")
        ax.grid(True, alpha=0.25)
        ax.legend(frameon=False, fontsize=8)

    fig.suptitle(
        "Figure 9: Thomas-Fiering Monthly Mean Validation\n"
        "Historical vs Synthetic Mean Monthly Inflows (Oct through Sep)",
        fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    _save(fig, "figure09_TF_monthly_means_validation.png", save_folder)


def plot_figure10_future_timeseries(system_future, save_folder):
    """
    Figure 10 -- Predicted system outflow time series 2026-2100 (Figure 10 in paper).

    Each line represents one inflow-reduction scenario.  The vertical spread
    between lines reflects sensitivity of predicted outflows to inflow reductions.
    Year-to-year variability within each line reflects the interannual variability
    embedded in the TF synthetic inflow sequences -- the synthetic series preserves
    historical variance, so wet and dry years still appear even under reduced-inflow
    scenarios.  The pattern shows that while inflow reductions shift the mean level
    downward, the nature of interannual variability is broadly maintained.
    """
    fig, ax = plt.subplots(figsize=(13, 6))
    for scen in SCENARIOS:
        sub = system_future[system_future["scenario"] == scen].sort_values("WY")
        ax.plot(sub["WY"], sub["predicted_outflow_af_day"], lw=1.8, label=scen)
    ax.set_title(
        "Figure 10: Projected Future SWP-CVP System Outflows Under Synthetic "
        "Inflow Reduction Scenarios\nPeriod: 2026-2100",
        fontweight="bold")
    ax.set_xlabel("Future water year")
    ax.set_ylabel("Predicted aggregated outflow (AF/day)")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    plt.tight_layout()
    _save(fig, "figure10_future_outflow_timeseries.png", save_folder)


def plot_figure11_statistical_summary(system_future, summary, system_wy,
                                      slope_a, intercept_c, r, r2, save_folder):
    """
    Figure 11 -- Three-panel statistical summary of scenario outflows
    (Figure 11 in paper).

    Panel (a): The Step-1 regression line overlaid with future scenario
    supply-outflow points.  Because all scenario supply values are calculated
    from synthetic inflows derived from the same historical statistics, the
    scenario points cluster tightly along the regression line.  This confirms
    that the methodology consistently applies the empirical system relationship.

    Panel (b): Mean predicted outflow per scenario (bar chart with +/- 1 SD
    error bars).  The progressive decline from Baseline to -50% illustrates
    the structural dependence of Southern California deliveries on upstream
    inflow availability even when carryover storage is included in the supply.

    Panel (c): Boxplot of the predicted outflow distribution per scenario.
    The box spans the interquartile range; the whiskers extend to 1.5*IQR;
    circles beyond the whiskers are outliers driven by extreme synthetic years.
    As inflow reductions increase, the distributions shift toward lower values
    and become more concentrated, reflecting the narrowing of the range of
    possible outcomes when the supply is systematically constrained.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # -- Panel (a): regression fit + scenario scatter --------------------------
    ax = axes[0]
    fd = system_wy.dropna(subset=["X_supply_af_day", "Y_outflow_af_day"])
    ax.scatter(fd["X_supply_af_day"], fd["Y_outflow_af_day"],
               s=35, color=C_BLUE, zorder=4, label="Historical water years")
    xmin = min(system_future["scenario_supply_af_day"].min(),
               fd["X_supply_af_day"].min())
    xmax = max(system_future["scenario_supply_af_day"].max(),
               fd["X_supply_af_day"].max())
    xl   = np.linspace(xmin * 0.98, xmax * 1.02, 300)
    ax.plot(xl, intercept_c + slope_a * xl, lw=2.2, color=C_GREY,
            label="Step 1 linear fit")
    for scen in SCENARIOS:
        sub = system_future[system_future["scenario"] == scen]
        ax.scatter(sub["scenario_supply_af_day"], sub["predicted_outflow_af_day"],
                   s=14, alpha=0.55, label=scen)
    ax.set_title("(a) Step 1 fit with future scenario estimates",
                 fontweight="bold")
    ax.set_xlabel("Aggregated system supply (AF/day)")
    ax.set_ylabel("Predicted / observed system outflow (AF/day)")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=7, ncol=2)
    ax.text(0.03, 0.97,
            f"Y = {intercept_c:.2f} + {slope_a:.4f} X\n"
            f"r = {r:.3f},  R\u00b2 = {r2:.3f}",
            transform=ax.transAxes, va="top", ha="left", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", alpha=0.9))

    # -- Panel (b): mean +/- SD bar chart -------------------------------------
    ax = axes[1]
    xi  = np.arange(len(summary))
    ax.bar(xi, summary["mean"], yerr=summary["std"],
           capsize=5, color=C_BLUE, edgecolor="white", error_kw={"lw": 1.5})
    ax.set_xticks(xi)
    ax.set_xticklabels(summary["scenario"], rotation=22, ha="right")
    ax.set_ylabel("Predicted aggregated outflow (AF/day)")
    ax.set_title("(b) Mean future outflow by scenario", fontweight="bold")
    ax.grid(True, axis="y", alpha=0.25)

    # -- Panel (c): boxplot distribution -------------------------------------
    ax = axes[2]
    box_data = [
        system_future.loc[system_future["scenario"] == s,
                          "predicted_outflow_af_day"].values
        for s in SCENARIOS
    ]
    bp = ax.boxplot(box_data, tick_labels=list(SCENARIOS.keys()),
                    patch_artist=True, medianprops={"color": "white", "lw": 2})
    for patch in bp["boxes"]:
        patch.set_facecolor(C_BLUE)
        patch.set_alpha(0.65)
    ax.set_ylabel("Predicted aggregated outflow (AF/day)")
    ax.set_title("(c) Distribution of future outflow scenarios",
                 fontweight="bold")
    ax.grid(True, axis="y", alpha=0.25)
    for tick in ax.get_xticklabels():
        tick.set_rotation(22)
        tick.set_ha("right")

    fig.suptitle(
        "Figure 11: Statistical Evaluation of Scenario-Based Future SWP-CVP "
        "System Outflows\n"
        "Derived from Thomas-Fiering Synthetic Inflows and Step-1 "
        "Supply-Outflow Relationship",
        fontsize=12, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    _save(fig, "figure11_future_outflow_statistical_summary.png", save_folder)


# ===========================================================================
# ==========================================================================
#
#   SECTION 11 -- MAIN PIPELINE
#   Orchestrates all steps in order.  Each major block prints progress so
#   a user watching the terminal can follow what is happening.
#
# ===========================================================================
# ===========================================================================

def main(save_folder=SAVE_FOLDER,
         hist_start=HIST_START,
         hist_end=HIST_END,
         skip_download=False):
    """
    Execute the complete analysis pipeline from raw CDEC data to all figures.

    Parameters
    ----------
    save_folder    : Path where all outputs are written.  Created if needed.
    hist_start     : Start date of the historical CDEC download window.
    hist_end       : End date of the historical CDEC download window.
    skip_download  : If True and cached WY summary CSVs exist, skip the CDEC
                     download and reuse them.  Useful for regenerating figures
                     after a full run without waiting for CDEC again.
    """
    os.makedirs(save_folder, exist_ok=True)

    print("=" * 70)
    print("Estimating Future SWP Deliveries to Southern California")
    print(f"Output : {save_folder}")
    print(f"Period : {hist_start}  to  {hist_end}")
    print(f"Future : WY {FUTURE_START_WY} -- {FUTURE_END_WY}")
    print("=" * 70)

    # -----------------------------------------------------------------------
    # FIGURE 1 (no CDEC download needed -- data is hardcoded above)
    # -----------------------------------------------------------------------
    print("\n[Fig 1] SWP allocation trends ...")
    plot_figure01_swp_allocation(save_folder)

    # -----------------------------------------------------------------------
    # STEP 1 -- DOWNLOAD OR LOAD CACHED DATA
    # -----------------------------------------------------------------------
    wy_csv  = os.path.join(save_folder, "WY_reservoir_summary_all.csv")
    sys_csv = os.path.join(save_folder, "system_supply_outflow_WY_summary.csv")

    if skip_download and os.path.exists(wy_csv) and os.path.exists(sys_csv):
        print("\n[Step 1] Loading cached WY summaries (--skip-download) ...")
        wy_all    = pd.read_csv(wy_csv)
        system_wy = pd.read_csv(sys_csv)
        inflows_daily   = {}
        outflows_daily  = {}
        storage_monthly = {}
        in_corr = out_corr = car_corr = None
    else:
        print("\n[Step 1] Downloading CDEC data ...")
        inflows_daily, outflows_daily, storage_monthly = \
            download_all_cdec_data(hist_start, hist_end)

        print("\n[Step 1] Saving raw Excel workbooks ...")
        save_raw_workbooks(inflows_daily, outflows_daily, storage_monthly,
                           hist_start, hist_end, save_folder)

        print("\n[Step 1] Building water-year summaries ...")
        wy_all, system_wy = build_wy_summaries(
            inflows_daily, outflows_daily, storage_monthly, save_folder)

        print("\n[Fig 2-4] Daily inflow / outflow / storage panels ...")
        plot_figure02_daily_inflows(inflows_daily, save_folder)
        plot_figure03_daily_outflows(outflows_daily, save_folder)
        plot_figure04_storage(storage_monthly, save_folder)

        print("\n[Step 1] Computing Pearson correlation matrices ...")
        in_corr, out_corr, car_corr = compute_pearson_matrices(
            inflows_daily, outflows_daily, storage_monthly, save_folder)

        print("\n[Fig 5] Pearson correlation matrices ...")
        plot_figure05_correlation_matrices(in_corr, out_corr, car_corr, save_folder)

    print("\n[Step 1] Fitting individual reservoir relationships ...")
    fit_df = fit_individual_reservoirs(wy_all, save_folder)
    print(fit_df[["station", "pearson_r", "r2", "n"]].to_string(index=False))

    print("\n[Fig 6] Individual reservoir linear fits ...")
    plot_figure06_individual_fits(wy_all, save_folder)

    print("\n[Step 1] Fitting system-level relationship ...")
    slope_a, intercept_c, r, r2 = fit_system_level(system_wy, save_folder)

    print("\n[Fig 7] System-level linear fit ...")
    plot_figure07_system_fit(system_wy, slope_a, intercept_c, r, r2, save_folder)

    # -----------------------------------------------------------------------
    # STEP 2 -- THOMAS-FIERING SYNTHETIC INFLOW GENERATION
    # -----------------------------------------------------------------------
    # If download was skipped, try to reload the raw inflow workbook or the
    # previously generated per-reservoir synthetic CSVs so TF can still run.
    if not inflows_daily:
        wb = os.path.join(save_folder, f"INFLOWS_DAILY_{hist_start}_{hist_end}.xlsx")
        if os.path.exists(wb):
            print("\n[Step 2] Reloading inflows from saved workbook ...")
            xl = pd.ExcelFile(wb)
            for sta in ORDER:
                sheet = f"{sta}_inflow"
                if sheet in xl.sheet_names:
                    df = xl.parse(sheet, index_col=0, parse_dates=True)
                    if "inflow_af_day" not in df.columns and "inflow_cfs" in df.columns:
                        df = cfs_to_af_day(df, "inflow_cfs", "inflow_af_day")
                    inflows_daily[sta] = df
        # Final fallback: rebuild from per-reservoir synthetic CSVs if present.
        if not inflows_daily:
            print("  [INFO] Building proxy inflows from existing synthetic CSVs ...")
            for sta in ORDER:
                p = os.path.join(
                    save_folder,
                    f"{sta}_thomas_fiering_synthetic_monthly_inflow.csv")
                if os.path.exists(p):
                    df = pd.read_csv(p, parse_dates=["datetime"]).set_index("datetime")
                    df["inflow_af_day"] = (df["synthetic_inflow_af_month"]
                                           / df.index.days_in_month)
                    inflows_daily[sta] = df

    print("\n[Step 2] Thomas-Fiering: fitting parameters and generating synthetic inflows ...")
    synthetic_dict, tf_summary = run_thomas_fiering(inflows_daily, save_folder)
    print("\n  TF summary:")
    print(tf_summary[["station", "observed_months", "ks_p_value",
                       "low_confidence"]].to_string(index=False))

    print("\n[Fig 8] TF historical vs synthetic WY medians ...")
    plot_figure08_tf_wy_medians(synthetic_dict, save_folder)

    print("\n[Fig 9] TF monthly mean validation ...")
    plot_figure09_tf_monthly_means(synthetic_dict, save_folder)

    # -----------------------------------------------------------------------
    # STEP 2 -- CARRYOVER STORAGE AND SCENARIO CALCULATIONS
    # -----------------------------------------------------------------------
    # If storage_monthly is empty (skip-download path), re-download it now.
    if not storage_monthly:
        print("\n[Step 2] Downloading monthly storage for carryover calculation ...")
        for sta in ORDER:
            s   = RESERVOIRS[sta]["storage"]
            dfs = fetch_cdec_series(sta, s, hist_start, hist_end, "M",
                                    chunk_yearly=False)
            if not dfs.empty:
                dfs = dfs.rename(columns={"value": "storage_af"})
                dfs["storage_af"] = pd.to_numeric(dfs["storage_af"], errors="coerce")
                storage_monthly[sta] = dfs

    print("\n[Step 2] Computing historical median September carryover storage ...")
    carryover_df = download_carryover_storage(hist_start, hist_end, storage_monthly)
    carryover_df.to_csv(
        os.path.join(save_folder, "historical_median_carryover_storage.csv"),
        index=False)
    print(carryover_df[["station", "median_sept_storage_af",
                         "carryover_equiv_af_day"]].to_string(index=False))

    print("\n[Step 2] Building future scenarios and predicting outflows ...")
    system_future, scenario_summary = build_future_scenarios(
        synthetic_dict, carryover_df, slope_a, intercept_c, save_folder)

    print("\n  Scenario summary (AF/day):")
    print(scenario_summary.to_string(index=False))

    print("\n[Fig 10] Future outflow time series ...")
    plot_figure10_future_timeseries(system_future, save_folder)

    print("\n[Fig 11] Statistical summary of future scenarios ...")
    plot_figure11_statistical_summary(
        system_future, scenario_summary, system_wy,
        slope_a, intercept_c, r, r2, save_folder)

    # -----------------------------------------------------------------------
    # FINAL SUMMARY
    # -----------------------------------------------------------------------
    print()
    print("=" * 70)
    print("Complete.  All outputs written to:")
    print(f"  {save_folder}")
    print()
    print("FIGURES")
    figs = [
        "figure01_SWP_allocation.png / .pdf",
        "figure02_daily_inflows_8panel.png",
        "figure03_daily_outflows_8panel.png",
        "figure04_storage_8panel.png",
        "figure05_pearson_correlation_matrices.png",
        "figure06_individual_reservoir_fits.png",
        "figure07_system_level_fit.png",
        "figure08_TF_historical_vs_synthetic_WY_medians.png",
        "figure09_TF_monthly_means_validation.png",
        "figure10_future_outflow_timeseries.png",
        "figure11_future_outflow_statistical_summary.png",
    ]
    for f in figs:
        print(f"  {f}")
    print()
    print("KEY CSVs")
    csvs = [
        "WY_reservoir_summary_all.csv",
        "system_supply_outflow_WY_summary.csv",
        "aggregated_system_fit_summary.csv",
        "individual_reservoir_supply_outflow_fits.csv",
        "thomas_fiering_synthetic_summary_all_reservoirs.csv",
        "historical_median_carryover_storage.csv",
        "future_system_outflows_scenarios_2026_2100.csv",
        "future_outflow_scenario_summary.csv",
    ]
    for c in csvs:
        print(f"  {c}")
    print("=" * 70)


def download_carryover_storage(hist_start, hist_end, storage_monthly):
    """
    Compute historical median September carryover storage from an already-
    downloaded storage_monthly dict (avoids a second CDEC request).

    Returns a DataFrame with:
        station, reservoir,
        median_sept_storage_af,
        carryover_equiv_af_day
    """
    rows = []
    for sta in ORDER:
        if sta not in storage_monthly:
            continue
        df   = storage_monthly[sta]
        df["storage_af"] = pd.to_numeric(df["storage_af"], errors="coerce")
        sept = df[df.index.month == 9].dropna(subset=["storage_af"])
        if sept.empty:
            continue
        med = sept["storage_af"].median()
        rows.append({
            "station":                 sta,
            "reservoir":               RESERVOIRS[sta]["name"],
            "median_sept_storage_af":  med,
            "carryover_equiv_af_day":  med / DAYS_PER_YEAR,
        })
    if not rows:
        raise RuntimeError("No carryover storage computed.  "
                           "Check CDEC sensor IDs and storage data.")
    return pd.DataFrame(rows)


# ===========================================================================
# ==========================================================================
#
#   SECTION 12 -- COMMAND-LINE INTERFACE
#
# ===========================================================================
# ===========================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="swp_southern_california_analysis.py",
        description=(
            "Estimate future SWP deliveries to Southern California using "
            "Thomas-Fiering synthetic inflows and a system-level supply-outflow "
            "regression.  Produces 11 figures and all supporting CSVs."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples\n"
            "--------\n"
            "  # Full run (downloads CDEC and produces everything)\n"
            "  python swp_southern_california_analysis.py\n\n"
            "  # Write to a custom folder\n"
            "  python swp_southern_california_analysis.py --output ~/results\n\n"
            "  # Regenerate figures without re-downloading CDEC\n"
            "  python swp_southern_california_analysis.py --skip-download\n\n"
            "  # Custom date window\n"
            "  python swp_southern_california_analysis.py "
            "--start 1985-01-01 --end 2025-12-31\n"
        ),
    )
    parser.add_argument(
        "--output", default=SAVE_FOLDER,
        help=f"Output folder (default: {SAVE_FOLDER})",
    )
    parser.add_argument(
        "--start", default=HIST_START,
        help=f"Historical start date YYYY-MM-DD (default: {HIST_START})",
    )
    parser.add_argument(
        "--end", default=HIST_END,
        help=f"Historical end date YYYY-MM-DD (default: {HIST_END})",
    )
    parser.add_argument(
        "--skip-download", action="store_true",
        help="Reuse cached WY summary CSVs instead of re-downloading from CDEC",
    )
    args = parser.parse_args()
    main(
        save_folder=args.output,
        hist_start=args.start,
        hist_end=args.end,
        skip_download=args.skip_download,
    )
