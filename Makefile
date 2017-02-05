CC = g++
CFLAGS = -lm -pthread -Wall -O2 --std=c++0x
BUILDDIR := build
SRCDIR := sources
LIBS = -lopencv_core -lopencv_features2d -lopencv_flann -lopencv_imgcodecs -lopencv_imgproc -lopencv_photo
INCLUDE = -I/usr/local/opt/opencv3/include -L/usr/local/opt/opencv3/lib

all: dir main

dir :
	mkdir -p $(BUILDDIR)
main : $(SRCDIR)/main.cpp $(SRCDIR)/detector.cpp
	$(CC) $(SRCDIR)/main.cpp $(SRCDIR)/detector.cpp -o $(BUILDDIR)/detector $(CFLAGS) $(LIBS) $(INCLUDE)