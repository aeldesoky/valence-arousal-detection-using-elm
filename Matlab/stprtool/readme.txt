==================================================================
Statistical Pattern Recognition Toolbox for Matlab  
(C) 1999-2009, Version 2.13, 09-Jan-2016, 
Written by Vojtech Franc, Vaclav Hlavac,
{xfrancv,hlavac}@cmp.felk.cvut.cz

Czech Technical University Prague (http://www.cvut.cz)
Faculty of Electrical Engineering (http://www.feld.cvut.cz)
Center for Machine Perception (http://cmp.felk.cvut.cz)
==================================================================

0. Contents
===========
0. Contents
1. Introduction
2. Requirements
3. Installation
4. How to start ?

1. Introduction
================
This toolbox implements a selection of statistical pattern recognition 
methods described in the monograph Schlesinger M.I., Hlavac V.: Ten 
lectures on statistical and structural pattern recognition, Kluwer 
Academic Publishers, 2002. The basic idea of the toolbox was sketched 
by M.I. Schlesinger and V. Hlavac in June 1999. The design of the toolbox 
and implementation is by V. Franc who did it as his master thesis at the 
Czech Technical University Praha, Czech Republic (defended February 2000). 
The toolbox is still being developed and new implemented methods (e.g., 
SVM and other Kernel machines) go beyond the contents of the monograph. 

This software can be used freely for academic, non-profit purposes. 
If you intend to use it for commercial development, please, contact us.

2. Requirements
================
- Matlab, version 5.3 and higher.
- Several algorithms require the Optimization Toolbox and Images toolbox
  by MathWorks. However, most algorithms work with a plain Matlab.

3. Installation
================

- Create a directory of your choice and copy the toolbox there.
- The path to toolbox directories has to be set before using
  the toolbox. This can be done by running
    stprpath
  in the toolbox root. You can add the 'stprpath' command to 
  your 'startup.m' file. 
- The toolbox contains several algorithms implemented in C language
  which must be compiled by the MEX compiler. This algorithms have 
  already been compiled for Windows Matlab 5.3 and Linux Matlab 6.
  If necessary compile all C functions by running
    compilemex
  from the toolbox root.

4. How to start ?
==================

To get started, type 'help stprtool' or 'help stprtool/demos'
then you can run one of the listed demo programs.

An user's guide can be found in 'doc/stprtool.pdf'.

To see HTML documentation open 'doc/manual/index.html' (in the toolbox
root) in your Internet browser. 

If you have any comments or suggestions then, please, send me 
an email at xfrancv@cmp.felk.cvut.cz .
