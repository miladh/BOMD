include(../defaults.pri)


TEMPLATE = lib
TARGET = ../lib/bomd


SOURCES += \
    solver/bomd.cpp \
    fileManager/filemanager.cpp \
    generator/generator.cpp

HEADERS += \
    solver/bomd.h \
    fileManager/filemanager.h \
    generator/generator.h

OTHER_FILES += ../include/bomd.h ../install/include/bomd.h

!equals(PWD, $${OUT_PWD}) {
    QMAKE_POST_LINK += $(COPY_DIR) $$OUT_PWD/../lib $$PWD/../
}
