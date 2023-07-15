BINARIES := hll
all: $(BINARIES)

OBJ = main.o count_min_sketch.o pq_array.o

CFLAGS := $(CFLAGS) -lz -march=native -fopenmp -funroll-loops -fipa-bit-cp -ftree-loop-vectorize -ftree-slp-vectorize -ftracer -fsplit-loops -funswitch-loops -O

clean:
	rm -f *.o $(BINARIES)

tags:
	etags *.h *.c *.cc

%.o: %.cpp
	g++ -c $(CFLAGS) $< -o $@ -lz -march=native -fopenmp -funroll-loops -fipa-bit-cp -ftree-loop-vectorize -ftree-slp-vectorize -ftracer -fsplit-loops -funswitch-loops -O

%.o: %.cpp
	gcc -c $(CFLAGS) $< -o $@ -lz -march=native -fopenmp -funroll-loops -fipa-bit-cp -ftree-loop-vectorize -ftree-slp-vectorize -ftracer -fsplit-loops -funswitch-loops -O

hll: $(OBJ)
	g++ $(CFLAGS) $^ -o $@ -lz -march=native -fopenmp -funroll-loops -fipa-bit-cp -ftree-loop-vectorize -ftree-slp-vectorize -ftracer -fsplit-loops -funswitch-loops -O

clean:
	rm *.o
