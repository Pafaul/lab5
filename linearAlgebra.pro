TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    custom.cpp \
    integrator.cpp \
    LA.cpp \
    model.cpp \

HEADERS += \
    custom.h \
    integrator.h \
    LA.h \
    model.h \
