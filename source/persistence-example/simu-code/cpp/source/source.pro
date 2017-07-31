######################################################################
# Automatically generated by qmake (3.0) Sat Jul 22 11:46:17 2017
######################################################################

TEMPLATE = app
TARGET = simulator
INCLUDEPATH += .

# Input
HEADERS += alphamodelfunctions.hpp \
           CellState.hpp \
           gemodelfunctions.hpp \
           ModelParameters.hpp \
           partitioningfunctions.hpp \
           simparams.h \
           sizemodelfunctions.hpp \
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
