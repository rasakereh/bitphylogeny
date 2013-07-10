#!/bin/bash

INCDIRS="-I/pkgs/pylab-02/Python-2.5.4/include/python2.5 -I/pkgs/pylab-02/Python-2.5.4/lib/python2.5/site-packages/numpy/core/include"
#INCDIRS="-I/Library/Frameworks/Python.framework/Versions/Current/include/python2.6 -I/Library/Frameworks/Python.framework/Versions/Current//lib/python2.6/site-packages/numpy/core/include"

cython lda.pyx
gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing  $INCDIRS -o lda.so lda.c
