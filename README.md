dspFRET
=======
Matlab code for analyzing TTTR data from PicoQuant. Typical application is analysis of single pair FRET measurements with pulsed interleaved excitation (PIE-FRET).

File Structure
==============
analyzePIE. Main analysis engine
readHeader. Function that reads the header of a Symphotime TTTR binary file
readCounts2. Speedup of over 100x from PicoQuant parser.
calcMCS. calculates the MCS binning. Returns MCS data and the maps of photons into bins.