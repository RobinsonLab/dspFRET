dspFRET
=======
Matlab code for analyzing TTTR data from PicoQuant. Typical application is analysis of single pair FRET measurements with pulsed interleaved excitation (PIE-FRET).

File Structure
==============
- analyzePIE. Main analysis engine
- readHeader. Function that reads the header of a Symphotime TTTR binary file
- readCounts2. Speedup of over 100x from PicoQuant parser.
- calcMCS. calculates the MCS binning. Returns MCS data and the maps of photons into bins.

Output
======
TCSPC histogram
Photon counting histogram (PCH)
MCS traces (also called burst traces)
1D histograms of TE, S, and tau_DD
2D histograms of TE-S and TE and tau_DD