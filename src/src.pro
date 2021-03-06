include(../defaults.pri)


TEMPLATE = lib
TARGET = ../lib/bomd


SOURCES += \
    fileManager/filemanager.cpp \
    generator/generator.cpp \
    molecularSystem/molecularsystem.cpp \
    modifier/modifier.cpp \
    modifier/velocityrescaling.cpp

HEADERS += \
    fileManager/filemanager.h \
    generator/generator.h \
    molecularSystem/molecularsystem.h \
    modifier/modifier.h \
    modifier/velocityrescaling.h

OTHER_FILES += ../include/bomd.h ../install/include/bomd.h

!equals(PWD, $${OUT_PWD}) {
    QMAKE_POST_LINK += $(COPY_DIR) $$OUT_PWD/../lib $$PWD/../
}
