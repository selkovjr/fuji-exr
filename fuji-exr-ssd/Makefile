GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)

OBJ = sdd.o linear.o rotate.o cfa_mask.o io_tiff.o libdemosaic.o  libAuxiliary.o fuji-exr.o progressbar.o
BIN = fuji-exr
LIBBIN=.

CXX=g++

HDRDIR=-I/opt/local/include/ -I/usr/local/include/
LIBDIR=-L/opt/local/lib/ -L/usr/local/lib/


CFLAGS += -g -O3 -D VERSION=\"$(GIT_VERSION)\" \
  -fopenmp -funroll-loops -fomit-frame-pointer  -fno-tree-pre -falign-loops -ffast-math -ftree-vectorize \
  -Weffc++ -pedantic -Wall -Wextra  -Wno-write-strings -Wno-deprecated  $(HDRDIR)

LDFLAGS += -g $(CFLAGS) $(LIBDIR) -ltiff -lncurses -lgomp -lpthread


LIBMX=sdd.o linear.o rotate.o cfa_mask.o io_tiff.o libAuxiliary.o libdemosaic.o progressbar.o

default: $(OBJ) $(BIN)

$(OBJ) : %.o : %.cpp
	$(CXX) -c $(CFLAGS)  $< -o $@

$(BIN) : % : %.o  $(LIBMX)
	$(CXX) -o $(LIBBIN)/$@  $^ $(LDFLAGS)


.PHONY : clean
clean:
	$(RM) $(OBJ); cd $(LIBBIN); rm -f $(BIN)
