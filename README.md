# pysims

Magnetic-sector secondary ion mass spectroscopy (MS-SIMS) data analysis package.
It includes:
- processing routines for mass spectra, depth profile and energy spectra
- efficient loading procedures for CAMECA data files (mass spectra, depth profile and energy spectra) exported in the .txt format
- Mass_Spectrum_Monitor: a user interface for displaying mass spectra together with natural isotopic abundance of single ions and molecular ions combining up to 10 atoms (asking for more is possible but it comes at the expense of long calculation times).

## Installation:

Requires python >= 3.10

(in the `pysims` directory)
```
pip install -r requirements.txt
pip install .
```

## To create the corresponding jupyter ipykernel :
 (in the `pysims` directory)
```
venv/bin/python3 -m ipykernel install --name=pysims-venv
```
if permission errors occur, 
```
sudo venv/bin/python3 -m ipykernel install --name=pysims-venv
```

## Using Mass_Spectrum_Monitor

Run the interface from the MassSpectrumMonitor directory:
```
python3 mass_spectrum_monitor.py
```