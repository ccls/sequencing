# $Id: Makefile.winmasker.app 339995 2011-10-04 18:52:07Z vakatov $

WATCHERS = morgulis camacho mozese2

REQUIRES = objects algo

ASN_DEP = seq

APP = windowmasker
SRC = main win_mask_app win_mask_sdust_masker

LIB = xalgowinmask xalgodustmask blast composition_adjustment seqdb blastdb \
	seqmasks_io tables xobjread creaders xobjutil \
	$(OBJMGR_LIBS:%=%$(STATIC))

LIBS = $(CMPRS_LIBS) $(NETWORK_LIBS) $(DL_LIBS) $(ORIG_LIBS)

CXXFLAGS = $(FAST_CXXFLAGS)
LDFLAGS  = $(FAST_LDFLAGS)

PROJ_TAG = gbench