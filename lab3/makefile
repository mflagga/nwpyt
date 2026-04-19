ifeq ($(OS),Windows_NT)
	PC = python
	RM = del /q
	CFLAGS = -arch=native -O3
	LIBS = 
else
	PC = python
	RM = rm -f
	CFLAGS = -arch=native -Xcompiler -Wall,-Wextra -Werror all-warnings -O3
	LIBS = -lm
endif

.SUFFIXES:

CXX = nvcc
EXEC = main.out
SRCS = main.cu
HDRS = header.cuh
SCRIPT = plot.py
DATA = S1a.csv S1misc.csv S1rL.csv S3a.csv S3misc.csv S3r.csv S2a.csv S2misc.csv S2rV.csv
MEDIA = S1mask.png S1rho.png S1V.png S1rL.png S3mask.png S3rho.png S3V.png S3Vr.png S2mask.png S2rho.png S2V.png S2rV.png

all: $(MEDIA)

$(EXEC): $(SRCS) $(HDRS)
	$(CXX) $(SRCS) $(CFLAGS) $(LIBS) -o $@

$(DATA): $(EXEC)
	./$<

$(MEDIA): $(SCRIPT) $(DATA)
	$(PC) $(SCRIPT)

clean:
	$(RM) $(EXEC)
	$(RM) $(DATA)
	$(RM) $(MEDIA)

.PHONY: all clean
