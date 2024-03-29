######################################################################
# Automatically generated by qmake (2.01a) Sun Apr 24 21:30:41 2016
######################################################################

TEMPLATE = app
TARGET = simulator
DEPENDPATH += . libs
INCLUDEPATH += . libs

# Input
HEADERS += CellState.hpp \
           ModelParameters.hpp \
           simparams.h \
           StochSimulator.hpp \
           libs/common.h \
           libs/ludcmp.h \
           libs/nr3.h \
           libs/odeint.h \
           libs/ran.h \
           libs/sort.h \
           libs/sparse.h \
           libs/stepper.h \
           libs/stepperdopr5.h \
           libs/stepperdopr853.h \
           libs/steppersie.h
SOURCES += CellState.cpp \
           main.cpp \
           ModelParameters.cpp \
           simparams.cpp \
           StochSimulator.cpp \
           libs/common.cpp \
           libs/ludcmp.cpp \
           libs/sparse.cpp \
           libs/stepperdopr853.cpp
