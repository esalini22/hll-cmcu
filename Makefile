BINARIES := hll
all: $(BINARIES)

OBJ = main.o HyperCount.o CountMinCU.o pq_array.o

CFLAGS := $(CFLAGS) -lz -lm -pthread -Ofast -funroll-loops -fipa-bit-cp -ftree-loop-vectorize -ftree-slp-vectorize -ftracer -fsplit-loops -funswitch-loops -march=native -flto -fopenmp -D_GLIBCXX_PARALLEL -DNDEBUG 

clean:
	rm -f *.o $(BINARIES)

tags:
	etags *.h *.c *.cc

%.o: %.cpp
	g++ -c $(CFLAGS) $< -o $@ -lz -lm -pthread -Ofast -funroll-loops -fipa-bit-cp -ftree-loop-vectorize -ftree-slp-vectorize -ftracer -fsplit-loops -funswitch-loops -march=native -flto -fopenmp -D_GLIBCXX_PARALLEL -DNDEBUG 

%.o: %.cpp
	gcc -c $(CFLAGS) $< -o $@ -lz -lm -pthread -Ofast -funroll-loops -fipa-bit-cp -ftree-loop-vectorize -ftree-slp-vectorize -ftracer -fsplit-loops -funswitch-loops -march=native -flto -fopenmp -D_GLIBCXX_PARALLEL -DNDEBUG 

hll: $(OBJ)
	g++ $(CFLAGS) $^ -o $@ -lz -lm -pthread -Ofast -funroll-loops -fipa-bit-cp -ftree-loop-vectorize -ftree-slp-vectorize -ftracer -fsplit-loops -funswitch-loops -march=native -flto -fopenmp -D_GLIBCXX_PARALLEL -DNDEBUG 

clean:
	rm *.o