AUTOMAKE_OPTIONS = foreign

INCLUDES = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(OFFLINE_MAIN)/geant4/include \
  -I`root-config --incdir`

libCentauroJets_la_LDFLAGS = \
  -nodefaultlibs \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --libs`

lib_LTLIBRARIES = \
  libCentauroJets.la

noinst_HEADERS = \
  CentauroJetsLinkDef.h 

include_HEADERS = \
  CentauroJets.h

libCentauroJets_la_SOURCES = \
  CentauroJets.cc \
  CentauroJetsDict.C

libCentauroJets_la_LIBADD = \
  -lg4detectors \
  -lphg4hit \
  -lg4eval  \
  -lphool   \
  -lfastjet \
  -lfastjettools \
  -lCentauro


install-exec-hook:
	@[ $$OSTYPE = darwin ] && ln -sf $(DESTDIR)$(libdir)/libCentauroJets.dylib $(DESTDIR)$(libdir)/libCentauroJets.so || true

BUILT_SOURCES = \
  testexternals.C

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libCentauroJets.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

# Rule for generating CINT dictionaries from class headers.
%Dict.C: %.h %LinkDef.h
	rootcint -f $@ -c $(CINTFLAGS) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

testexternals_SOURCES = testexternals.C
