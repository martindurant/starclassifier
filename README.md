starclassifier
==============

Find best-match to UVBRIJHK photometry from library specta

This is a simple program which compares input optical/NIR broad-band photometry to a set of model spectra covering the whol M-K system.
The spectral model with closest correspondance is returned (together with the residuals).

General inputs: magnitudes in UVBRIJHK and their uncertainties.
Any band not provided is ignored, and equal uncertainties are assumed for all bands if not given.

Command line version
--------------------

Requirements:
- python 2.x
- numpy

Installation:
put star_phot.py and grid.npz in the same place, and execute. You may put it on your PATH if you wish.
If your python lives in a non-standard location, you can always invoke as

> python star_phot.py

GUI version
-----------

Includes plot of best-fit spectrum

Requires Qt and pyQt (or could be converted to pyside). 
For linux this is simple, for mac/windows, please google.

The main.py file can be made executable, or run with
> python main.py
