COBJ = io_png.o
CXXOBJ = libAuxiliary.o libdemosaic.o duran-buades.o imdiff_ipol.o

BIN = imdiff_ipol duran-buades

hdrdir= -I/usr/local/include/ -I/usr/include/
libdir= -L/usr/local/lib/ -L/usr/lib/

COPT = -g -O3 -fopenmp -funroll-loops -fomit-frame-pointer  -Wall -Wextra
CFLAGS  += $(COPT) $(hdrdir)

LDFLAGS +=  $(libdir) -lgomp -lpng

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

$(BIN) : % : %.o  io_png.o libAuxiliary.o libdemosaic.o
	$(CXX)  -std=c++98  -o $@  $^ $(LDFLAGS)

.PHONY : clean
clean:
	$(RM) $(COBJ) $(CXXOBJ) ; rm -f $(BIN)

sample: 151018_180328.RAF
	#dcraw -v -w -d -s all -6 -b 1.7 151018_180328.RAF
	dcraw -v -W -s all -4 151018_180328.RAF
	convert 151018_180328_0.pgm -crop 252x276+1116+1416 +repage fragment-0.png
	convert 151018_180328_1.pgm -crop 252x276+1116+1416 +repage fragment-1.png
	convert 151018_180328.JPG -crop 276x252+1576+1106 fragment.jpeg

test:
	./duran-buades fragment-0.png out-0.png 1
	./duran-buades fragment-1.png out-1.png 1
