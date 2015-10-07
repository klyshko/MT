GXX = g++ -O3 -Wall -Wno-unused-result

all: disc 3d22d

disc: pdbio.cpp dcdio.cpp disc.cpp
	$(GXX) $^ -o $@

3d22d: pdbio.cpp dcdio.cpp 3d22d.cpp
	$(GXX) $^ -o $@ 

.PHONY: all
