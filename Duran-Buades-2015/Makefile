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

sample: 151018_180328.RAF
	dcraw -v -w -W -d -s all -4 -T 151018_180328.RAF
	convert 151018_180328_0.tiff -crop 252x276+1116+1416 +repage fragment-0.tiff
	convert 151018_180328_1.tiff -crop 252x276+1116+1416 +repage fragment-1.tiff
	convert 151018_180328.JPG -crop 276x252+1576+1106 fragment.jpeg

test:
	./duran-buades fragment-0.tiff out-0.tiff 1
	./duran-buades fragment-1.tiff out-1.tiff 1
	hdr_create -4 -o test.exr out-*tiff
	hdr_convert -i test.exr -o test.tiff
