//This decodes DICOM images with Transfer Syntax 1.2.840.10008.1.2.4.70
// "JPEG Lossless, Nonhierarchical, First- Order Prediction"
// This format is described in http://www.w3.org/Graphics/JPEG/itu-t81.pdf
// It is identified with the 'Start of Frame' (SOF) code 0xC3
// It appears unique to medical imaging, and is not supported by most JPEG libraries
// http://www.dicomlibrary.com/dicom/transfer-syntax/
#ifndef _WRITE_TIFF_H_
#define _WRITE_TIFF_H_

#ifdef  __cplusplus
extern "C" {
#endif

int write_tiff_img (const char *fn, unsigned char *img, int nx, int ny, int bits, int frames, int isPlanar);

#ifdef  __cplusplus
}
#endif

#endif