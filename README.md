fuji-exr-ssd
============

Self-similarity Driven Demosaicking of Fuji EXR images

after [www.ipol.im/pub/art/2011/bcms-ssdd/](http://www.ipol.im/pub/art/2011/bcms-ssdd/)

Fuji EXR filter layout:

<img src="doc/image/fuji-cfm.png" width="405" height="357" />


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
./fuji-exr-ssd raw.tiff out.tiff
```

* `raw.tiff`  :  camera sensor data in 16-bit grayscale
* `output.tiff` :  demosaicked RGB output
