# $Id: Makefile.writedb_unit_test.app 294689 2011-05-25 20:37:36Z camacho $

APP = writedb_unit_test
SRC = writedb_unit_test

CPPFLAGS = $(ORIG_CPPFLAGS) $(BOOST_INCLUDE)
CXXFLAGS = $(FAST_CXXFLAGS)
LDFLAGS = $(FAST_LDFLAGS)

LIB_ = test_boost writedb seqdb xobjread xobjutil creaders blastdb \
       $(SOBJMGR_LIBS)
LIB = $(LIB_:%=%$(STATIC))
LIBS = $(CMPRS_LIBS) $(NETWORK_LIBS) $(DL_LIBS) $(ORIG_LIBS)

CHECK_REQUIRES = in-house-resources
CHECK_CMD = writedb_unit_test
CHECK_COPY = writedb_unit_test.ini data

WATCHERS = maning madden camacho