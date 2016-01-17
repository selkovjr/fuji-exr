fuji-exr-ssd
============

Self-similarity Driven Demosaicking of Fuji EXR images

after [www.ipol.im/pub/art/2011/bcms-ssdd/](http://www.ipol.im/pub/art/2011/bcms-ssdd/)

## ABOUT

This software was originally written by A. Buades <toni.buades@uib.es>
with contributions from Nicolas Limare.

Adapted to Fuji EXR by Gene Selkov <selkovjr@gmail.com>.

All files are distributed under the terms of the LGPLv3 license.


## OVERVIEW

This source code provides an implementation of the Self Similar
demosaicking algorithm, as described in IPOL
  http://www.ipol.im/pub/algo/bcms_self_similarity_driven_demosaicking/


## REQUIREMENTS

The code is written in ANSI C, and should compile on any system with
an ANSI C compiler.

The libtiff header and libraries are required on the system for
compilation and execution.



## COMPILATION

Simply use the provided makefile, with the command `make`.



## USAGE

```
dcraw -v -w -d -s all -6 -T -b 0.7 raw.RAF
./fuji-exr-ssd raw_[01].tiff out.tiff
```

* `raw.RAF`: Fuji raw image, shot in EXR high-res mode, or in P-mode
* `raw_[01].tiff`:  camera sensor data in 16-bit grayscale (Bayer), two frames
* `out.tiff`:  demosaicked RGB output, rotated 45 degrees

Presently supported camera orientations: landscape (horizontal), portrait (270CCW). Other orientations need more work (interleaving rules are different for each).
