A Demosaicking Algorithm with Adaptive Inter-channel Correlation (version 1)

Joan Duran, joan.duran@uib.es, Universitat de les Illes Balears (Spain)
Antoni Buades, toni.buades@uib.es, Universitat de les Illes Balears (Spain)


# OVERVIEW

This C source code is derived from the demo code from the Image Processing On
Line (IPOL) article "A Demosaicking Algorithm with Adaptive Inter-channel
Correlation" at

    http://www.ipol.im/pub/art/

The demo code is located at

    http://demo.ipol.im/demo/

This program reads a gray-scale TIFF image of the Bayer array and writes out
16-bit RGB TIFF.

# USAGE

Usage: duran-buades bayer.tiff output.tiff beta

bayer.tiff    :: gray-scale CFA image.
output.tiff   :: full color demosaicked image.
beta          :: fixed channel-correlation parameter.

The following parameters are fixed in the main demo function:
epsilon   : thresholding parameter avoiding numerical intrincacies when
            computing local variation of chromatic components.
M         : bounding parameter above which a discontinuity of the luminance
            gradient is considered.\n");
halfL     : half-size of the support zone where the variance of the chromatic
            components is computed.
reswind   : half-size of research window.
compwind  : half-size of comparison window.
N         : number of most similar pixels for filtering.
redx redy : coordinates of the first red value in CFA.


#LICENSE

All files are distributed under the terms of the LGPLv3 license.


# REQUIREMENTS

The code is written in ANSI C and C++, and should compile on any system with an
ANSI C/C++ compiler.

The libtiff header and libraries are required on the system for compilation and
execution.


# COMPILATION

Simply use the provided makefile, with the command 'make'.

