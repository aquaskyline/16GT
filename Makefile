CC = g++
LIBFLAG = -lpthread -lm -lz
CFLAGS = -O3 -funroll-loops -fomit-frame-pointer -maccumulate-outgoing-args -funroll-loops -static-libgcc -mpopcnt -fopenmp -fpermissive -w

all: snapshotSnpcaller bam2snapshot RegionIndexBuilder

%.o : %.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

snapshotSnpcaller: snapshotSnpcaller.o SNP.o SAMhandler.o FisherExactTest.o VariantCaller.o SnapshotHandler.o fisher.o likelihood_cache.o ycsq.o SNP_Caller.o SNPFunctions.o interpreter.o lib/lib.a
	$(CC) $(CFLAGS) $^ $(LIBFLAG) -o $@

RegionIndexBuilder: readIndex.o indexFunction.o buildRegionList.o
	$(CC) $(CFLAGS) $^ -o $@

bam2snapshot: bam2snapshot.o SAMhandler.o SNP.o CounterReader.o SnapshotHandler.o indexFunction.o readIndex.o lib/lib.a
	$(CC) $(CFLAGS) $^ $(LIBFLAG) -o $@

clean:
	rm -fr *.o */*.o
	rm -f snapshotSnpcaller bam2snapshot RegionIndexBuilder
