Faradaic Efficiency Analysis Tool for Polarization Decoupling Experiments
JVD_PEC_compiler_v1.1, Updated March 02, 2023
Used in the production of: https://pubs.acs.org/doi/10.1021/acs.jpclett.3c03051

A Python toolkit for analyzing electrochemical CO2 reduction data, combining chronoamperometry (CA) measurements with gas chromatography (GC) product analysis to calculate Faradaic efficiencies.

Overview

This tool processes electrochemical data from Gamry potentiostats and SRI gas chromatographs to:

    Extract and synchronize CA and GC data
    Calculate Faradaic efficiencies for H2, CO, CH4, C2H4, and C2H6
    Generate publication-quality plots and data compilations
    Analyze voltage-dependent product selectivity

Features
Data Processing

    Automatic extraction of CA data from Gamry .DTA files
    GC peak integration from SRI .log files (TCD for H2, FID for carbon products)
    Time-synchronized injection point matching between CA and GC data
    Automatic data cleaning and outlier removal

Analysis Capabilities

    Faradaic efficiency calculations for all major CO2 reduction products
    Current density normalization
    Power density calculations
    CO2RR/H2 selectivity ratio analysis
    Statistical averaging with standard deviation calculations

Visualization

    Stacked bar charts of product distributions
    Voltage-dependent Faradaic efficiency plots
    Power density vs. selectivity relationships
    Time-resolved current density traces

Requirements

    Python 2.7
    NumPy
    matplotlib
    pandas
    SciPy
    extractor_v2 module (custom data extraction utilities)

Installation

bash

pip install numpy matplotlib pandas scipy

Ensure the extractor_v2.py module is in your Python path.
Usage
Basic Analysis

python

from faradaic_efficiency_analyzer import *

# Analyze a single dataset
Analyze(directory0, 
        headerCA=2,           # CA file header lines
        headerGC=0,           # GC file header lines
        flowrate=4.85,        # Gas flow rate (sccm)
        t_inject_start=345,   # First injection time (s)
        GC_runtime=11.5,      # Time between injections (min)
        num_injections=3,     # Injections per condition
        datatype='separate',  # Data format
        saveas='experiment_1')

Batch Processing

python

# Process multiple folders with voltage series
compile_arrays(headerCA=58, 
               headerGC=None,
               flowrate=4.8,
               t_inject_start=390,
               GC_runtime=13,
               num_injections=2,
               datatype='separate')

Input Data Structure

bash

input_data/
├── experiment_1/
│   ├── *.DTA          # Gamry CA data files
│   ├── *TCD*.log      # H2 GC data
│   └── *FID*.log      # CO, CH4, C2H4, C2H6 GC data
├── experiment_2/
│   └── ...

Output Files

The tool generates several output files in the output/ directory:

    *_Compiled_injection_timepoints+FEs.csv: Complete dataset with all calculations
    *_Compiled_raw_CA_data.csv: Raw chronoamperometry data
    *_A1_avg.csv: Averaged Faradaic efficiencies
    *_A1_std.csv: Standard deviations
    *_volt_dependence_FE.png: Stacked bar chart visualization
    product_ratio.png: Selectivity analysis plots

Key Parameters
Electrochemical Setup

    electrode_area: 2.5 cm²
    Zreal: 57 Ω (for iR compensation)

GC Calibration

The tool includes built-in calibration curves for:

    H2: Linear calibration via TCD
    CO, CH4, C2H4, C2H6: FID calibrations

Data Quality Control

    Automatic removal of data points with:
        H2 FE > 100% or < 0%
        Total FE > 114% or < 70%
    Outlier detection and removal

Customization
Modifying Calibrations

Update the calibration() function with your GC-specific calibration curves:

python

def calibration(H2_integral, CO_integral, CH4_integral, C2H4_integral, C2H6_integral):
    H2_ppm = slope_H2 * H2_integral + intercept_H2
    # ... etc for other products

Adjusting Voltage Ranges

Modify the Vall array in compile_arrays():

python

Vall = np.arange(-3, -5, -0.5)  # Forward scan
# or
Vall = np.arange(-4.0, -1.75, 0.25)  # Reverse scan

Author

Created by pcagbo (February 2020)
Notes

    Ensure proper file naming conventions for automatic TCD/FID detection
    The tool assumes specific column positions in GC log files - verify these match your instrument output
    For best results, maintain consistent injection timing across experiments

