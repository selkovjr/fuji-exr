from ij import IJ, ImagePlus
from ij.process import ShortProcessor  
from ij.io import FileSaver
from array import zeros  
import os
import subprocess

#basename = '141006_163724'
basename = '141016_155913'
frame = []

for fn in ['0', '1']:
    filename = os.environ['HOME'] + '/' + basename + '_' + fn + '.pgm'
    imp = IJ.openImage(filename)
    if imp is None:
        print "Could not open image from file:", filename
    frame.append(imp)
print 'loaded frames'

width = frame[0].width
height = frame[0].height
bpix = zeros('h', (width + height)**2)

for fn in [0, 1]:
    pix = frame[fn].getProcessor().getPixels()
    for i in range(len(pix)):
        x = i % width + i / width + fn # the second frame (fn == 1) is shifted 1px to the right
        y = (width - i % width - 1) + (i / width)
        bpix[y * (width + height) + x] = pix[i]
    print 'wrote frame %s' % fn

ip = ShortProcessor(width + height, width + height, bpix, None)  
bayer = ImagePlus("Sensor", ip)

fs = FileSaver(bayer)  
fs.saveAsTiff(os.path.dirname(filename) + '/' + 'raw.tiff')

subprocess.call(['/opt/local/bin/tiffset', '-s', '270', 'width = ' + str(width) + ', height = ' + str(height), os.environ['HOME'] + '/' + 'raw.tiff'])

bayer.show()
 