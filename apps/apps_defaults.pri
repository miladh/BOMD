TEMPLATE = app

LIBS += -lhdf5 -lhdf5_cpp
LIBS += -lunittest++
LIBS += -lboost_filesystem -lboost_system -lboost_mpi -lboost_serialization
LIBS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile)

LIBS += -L$$TOP_OUT_PWD/lib -lbomd
INCLUDEPATH += $$TOP_PWD/include

