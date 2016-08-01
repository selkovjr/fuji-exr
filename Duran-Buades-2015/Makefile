#COBJ = io_tiff.o
CXXOBJ = io_tiff.o libAuxiliary.o libdemosaic.o duran-buades.o

BIN = duran-buades

hdrdir= -I/usr/local/include/ -I/usr/include/
libdir= -L/usr/local/lib/ -L/usr/lib/

COPT = -g -O3 -fopenmp -funroll-loops -fomit-frame-pointer  -Wall -Wextra
CFLAGS  += $(COPT) $(hdrdir)

LDFLAGS +=  $(libdir) -lgomp -ltiff

# use openMP with `make OMP=1`
ifdef OMP
CFLAGS += -fopenmp
LDFLAGS += -lgomp
endif

default: $(COBJ) $(CXXOBJ)  $(BIN)

$(COBJ) : %.o : %.c
	$(CC) -c $(CFLAGS)   $< -o $@

$(CXXOBJ) : %.o : %.cpp
	$(CXX) -std=c++98 -c $(CFLAGS)   $< -o $@

$(BIN) : % : %.o  io_tiff.o libAuxiliary.o libdemosaic.o
	$(CXX)  -std=c++98  -o $@  $^ $(LDFLAGS)

.PHONY : clean
clean:
	$(RM) $(COBJ) $(CXXOBJ) ; rm -f $(BIN)

sample-normal: 151018_180303.RAF
	dcraw -v -w -W -d -s all -4 -T 151018_180303.RAF
	convert 151018_180303_0.tiff -crop 276x252+1756+1286 +repage fragment-0.tiff
	convert 151018_180303_1.tiff -crop 276x252+1756+1286 +repage fragment-1.tiff
	convert 151018_180303.JPG -crop 276x252+1756+1278 fragment.jpeg

sample-270: 151018_180328.RAF
	dcraw -v -w -W -d -s all -4 -T 151018_180328.RAF
	convert 151018_180328_0.tiff -crop 252x276+1116+1416 +repage fragment-0.tiff
	convert 151018_180328_1.tiff -crop 252x276+1116+1416 +repage fragment-1.tiff
	convert 151018_180328.JPG -crop 276x252+1576+1106 fragment.jpeg

test: duran-buades
	./duran-buades fragment-0.tiff out-0.tiff `exiftool -Orientation -n 151018_180328.RAF | cut -f2 -d:` 0
	./duran-buades fragment-1.tiff out-1.tiff `exiftool -Orientation -n 151018_180328.RAF | cut -f2 -d:` 0
	hdr_create -4 -o test.exr out-*tiff
	#hdr_denoise -r 1 -d 2 -i test.exr -o denoised.exr
	hdr_denoise -r 2 -d 4 -i test.exr -o denoised.exr
	hdr_convert -s 0.2 -i test.exr -o test.tiff
	hdr_convert -s 0.2 -i denoised.exr -o denoised.tiff
	convert -average test.tiff denoised.tiff denoised.tiff average.tiff
	convert average.tiff -level 3% -modulate 100,140,106 -gamma 1.14 -sharpen 0x1 color.tiff

fulltest: duran-buades
	time ./duran-buades 151018_180328_0.tiff decoded-0.tiff `exiftool -Orientation -n 151018_180328.RAF | cut -f2 -d:` 0
	time ./duran-buades 151018_180328_1.tiff decoded-1.tiff `exiftool -Orientation -n 151018_180328.RAF` 0
	time hdr_create -4 -o fulltest.exr decoded-*tiff
	time hdr_denoise -r 2 -d 4 -i fulltest.exr -o fulltest-denoised.exr
	time hdr_convert -s 0.28 -i fulltest.exr -o fulltest.tiff
	time hdr_convert -s 0.28 -i fulltest-denoised.exr -o fulltest-denoised.tiff
	time convert -average fulltest.tiff fulltest-denoised.tiff fulltest-average.tiff
	time convert fulltest-average.tiff -level 3% -modulate 102,180,106 -gamma 1.14 -sharpen 0x1 fulltest.tiff

