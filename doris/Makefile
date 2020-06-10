NVCC ?= nvcc #-g -G 
NVFLAGS ?= -std=c++11 -O3 -Xcompiler -fopenmp -lineinfo -gencode arch=compute_70,code=sm_70 --ptxas-options=-v

# Compile cuda source codes to objects                      
# use -dc when device function is at a different file      

%.o: %.cu
	$(NVCC) $(NVFLAGS) -dc -c -o $@ $<

# Compile objects to executables                                                         
 
vol: voltermMain.o VlasovSer2x3vP1.o VlasovSer2x3vP1_v1.o
	$(NVCC) $(NVFLAGS) -o $@ voltermMain.o VlasovSer2x3vP1.o VlasovSer2x3vP1_v1.o

volInline: voltermMain.cu
	$(NVCC) $(NVFLAGS) --maxrregcount=255 -DINLINE -o $@ $< -I./

surf: surfTermMain.o VlasovSurfElcMagSer2x3vP1.o VlasovSurfElcMagSer2x3vP1_v1.o VlasovSurfElcMagSer2x3vP1_v2.o
	$(NVCC) $(NVFLAGS) -o $@ surfTermMain.o VlasovSurfElcMagSer2x3vP1.o VlasovSurfElcMagSer2x3vP1_v1.o VlasovSurfElcMagSer2x3vP1_v2.o

# Compile all executables                                                                                            
 
all: vol volInline surf

# Clean up                                                                                                  
 
.PHONY: clean
clean:
	rm -f *.o vol volInline surf
 
