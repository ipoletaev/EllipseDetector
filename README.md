# EllipseDetector
C++ solution for ellipse detection on images.
The key idea is based on this [article](http://ieeexplore.ieee.org/document/1048464/).

### Requirements:
* OpenCV: 2.x (or 3.x) installed on your PC
* g++ compiler: 4.2 or higher

### Build:
```
cd ../EllipseDetetor
make
```

### Examples of usage
```
cd build/
./detector ../data/example.png ../data/example_out.png
```
All of the advanced options you can tune in the `../sources/detector.cpp` file. 