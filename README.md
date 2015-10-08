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
dcraw -v -w -d -s -all -6 -T -b 0.7 raw.RAF
./fuji-exr-ssd raw_[01].tiff out.tiff
```

* `raw.RAF`: Fuji raw image, shot in EXR high-res mode, or in P-mode
* `raw_[01].tiff`:  camera sensor data in 16-bit grayscale (Bayer), two frames
* `out.tiff`:  demosaicked RGB output, rotated 45 degrees

Presently supported camera orientations: landscape (horizontal), portrait (270CCW). Other orientations need more work (interleaving rules are different for each).


## REFERENCES / FOOD FOR THOUGHT

date | authors | title
---|---|---
2015-07-07 | Joan Duran, Antoni Buades | [A Demosaicking Algorithm with Adaptive Inter-channel Correlation](http://www.ipol.im/pub/pre/145/preprint.pdf)
2014 | Alan Gibson | [Demosaicing with ImageMagic](http://im.snibgo.com/demosaic.htm)
2012-08-01 | Eric Dubois and Gwanggil Jeon | [Demosaicking of Noisy Bayer-Sampled Color Images with Least-Squares Luma-Chroma Demultiplexing and Noise Level Estimation: Additional Results](http://www.site.uottawa.ca/~edubois/lslcd_ne/)
2011-06-01 | Antoni Buades, Bartomeu Coll, Jean-Michel Morel, and Catalina Sbert | [Self-similarity Driven Demosaicking, Image Processing On Line, 1 (2011)](http://dx.doi.org/10.5201/ipol.2011.bcms-ssdd)
2010 | Brian Leung, Gwanggil Jeon, and Eric Dubois | [Least-Squares Luma-Chroma Demultiplexing Algorithm for Bayer Demosaicking](http://www.site.uottawa.ca/~edubois/lslcd/article/TIP-06195-2010.R1_2col.pdf)
 | Craig Stark | [Debayering Demystified](http://www.stark-labs.com/craig/resources/Articles-&-Reviews/Debayering_API.pdf)
 | Henrique S. Malvar, Li-wei He, and Ross Cutler | [High-Quality Linear Interpolation for Demosaicing of Bayer-Patterned Color Images](http://research.microsoft.com/pubs/102068/Demosaicing_ICASSP04.pdf)
 | Sander Pool, Zbynek Vrastil | [PixInsight Reference Documentation: Debayer](http://pixinsight.com/doc/tools/Debayer/Debayer.html)
2010-01 | Robert A. Maschal Jr., S. Susan Young, Joe Reynolds, Keith Krapels, Jonathan Fanning, and Ted Corbin | [Review of Bayer Pattern Color Filter Array (CFA) Demosaicing with New Quality Assessment Algorithms](http://www.arl.army.mil/arlreports/2010/ARL-TR-5061.pdf)
2007 | Sira Ferradans, Marcelo Bertalm ́ıo and Vicent Caselles | [Geometry based Demosaicking](http://www.gpi.upf.edu/static/sira/Sira_Ferradans/Demosaicking_files/GeometrybasedDemosaicking.pdf)
2007-09-05 | Antoni Buades, Bartomeu Coll, Jean-Michel Morel, Catalina Sbert | [Non local demosaicing](http://dmi.uib.es/~tami/publicacions/CMLA2007-15.pdf)
2004 | Eric P. Bennett, Matthew Uyttendaele, C. Lawrence Zitnick, Richard Szeliski, and Sing Bing Kang | [Video and Image Bayesian Demosaicing with a Two Color Image Prior](http://research.microsoft.com/en-us/um/people/larryz/bennett-eccv06.pdf)
