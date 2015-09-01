include ../make.config

HDF5_C_INCLUDE=-I /usr/local/hdf5/intel/15/mvapich2/2.1a/1.8.13/include
HDF5_C_LIBS=-L /usr/local/hdf5/intel/15/mvapich2/2.1a/1.8.13/lib -lpthread -lssl -lcrypto -lz -lm

#CC = $(CCSERIAL)
CC = mpicc
CFLAGS = $(CFLAGSSERIAL)
INCLDIRS = $(HDF5_C_INCLUDE) -I /usr/include -I ../
LFLAGS = -lm $(HDF5_C_LIBS) -lhdf5 -lhdf5_hl # -lsvml
SOURCES = driver.c pcsource.c tree.c driver_correlate.c correlate.c bins.c pcsource_correlate.c tree_correlate.c
OBJECTS = $(SOURCES:.c=.o)
HEADERS = precompute.h globals.h correlation.h ../params.h
EXEC = correlate

.c.o:
	$(CC) $(CFLAGS) $(INCLDIRS) -c $<

all: $(EXEC)

$(EXEC): $(OBJECTS)
	$(CC) -o $(EXEC) $(OBJECTS) $(LFLAGS)

$(OBJECTS): $(HEADERS) Makefile ../make.config

clean:
	rm -f $(EXEC) *~ ../include/*~ $(OBJECTS)

