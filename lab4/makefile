ifeq ($(OS),Windows_NT)
	PC = python
	RM = del /q
	VR = echo "nie ma valgrinda na windowsie"
	AVE = ../nwenv/Scripts/activate
else
	ifdef VIRTUAL_ENV
		PC = python
	else
		PC = python3
	endif
	RM = rm -f
	VR = valgrind ./$(EXEC)
	AVE = source ../nwenv/bin/activate
endif

.SUFFIXES:

CXX = gcc
LIBS = -lm
CFLAGS = -Wall -Wextra -Werror -Wpedantic -O3
SRCS = main.c
HDRS = header.h
EXEC = main.out
SCRIPT = plot.py
DATA = p1.csv p1misc.csv fpsfile.csv p2.csv p2gamma.csv p3.csv p3gamma.csv
FRAMES = p1frames/framesdone
MEDIA = p1anim.mp4 p2anim.mp4
FPS = $(shell cat fpsfile.csv)

all: $(MEDIA)

$(EXEC): $(SRCS) $(HDRS)
	$(CXX) $(SRCS) $(CFLAGS) $(LIBS) -o $@

$(DATA): $(EXEC)
	./$<

$(FRAMES): $(SCRIPT) $(DATA)
	$(RM) p1frames/*.png p1frames/framesdone
	$(RM) p2frames/*.png
	$(PC) $<
	touch $@

p1anim.mp4: $(FRAMES)
	ffmpeg -framerate $(FPS) -i p1frames/frame_%04d.png -y -loglevel quiet -crf 18 p1anim.mp4

p2anim.mp4: $(FRAMES)
	ffmpeg -framerate $(FPS) -i p2frames/frame_%04d.png -y -loglevel quiet -crf 18 p2anim.mp4

clean:
	$(RM) $(EXEC)
	$(RM) $(DATA)
	$(RM) $(MEDIA)
	$(RM) $(FRAMES) p1frames/*.png p2frames/*.png p3frames/*.png

safe: $(EXEC)
	$(VR)

ve: 
	@echo Run: $(AVE)

gotowe: $(EXEC)
	./$(EXEC)

.PHONY: all clean safe gotowe ve
