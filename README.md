# BTsim
MT Bundle transport simulator

A stochastic particle-based simulator of Early Endosome (EE) transport on a Microtubule (MT) bundle. The bundles consist of multiple protofilaments on which EE run, and have variable geometry.

Code by Congping Lin and Peter Ashwin (2018)
Based on research with Gero Steinberg and Martin Schuster 

Available under GNU GPL v3.0

# Installing the code
This graphical user interface requires installation of Qt Creator which can be found on https://www.qt.io/download. Once installed, you can create an application project with Qt Widgets Applications. To run this program, you need to copy all the files into folder of newly created project. The program has been tested on QT version 5.5.1.

Please note that the program uses qcustomplot, and one needs to add a line (below) to the .pro file
QT       += printsupport
