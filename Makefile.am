ACLOCAL_AMFLAGS= -I ./make/m4

SUBDIRS=method_lib

EXTRA_DIST=CREATED

CLEANFILES=CREATED

dist-hook:
	chmod u+w $(distdir)/CREATED
	echo "Distribution date: `date`" > CREATED

# This is here to prevent installation of LOKI during a distcheck
DISTCHECK_CONFIGURE_FLAGS=--disable-loki

AM_CPPFLAGS=$(BOOST_CPPFLAGS) -I$(top_srcdir)/method_lib $(R_CFLAGS)
AM_LDFLAGS=$(BOOST_LDFLAGS) 

bin_PROGRAMS=plato

plato_SOURCES=Chrom.cc \
Chrom.h \
Controller.h \
dcdflib.cpp \
ExampleModule.cc \
ExampleModule.h \
Finalize.cc \
Finalize.h \
ipmpar.c \
PercentByFamily.cc \
PercentByFamily.h \
Percent.cc \
Percent.h \
ProcessAlleleFrequency.cc \
ProcessAlleleFrequency.h \
ProcessBEAGLEOutput.cc \
ProcessBEAGLEOutput.h \
ProcessBINOutput.cc \
ProcessBINOutput.h \
ProcessCaConChisq.cc \
ProcessCaConChisq.h \
Process.cc \
ProcessClusterMissing.cc \
ProcessClusterMissing.h \
ProcessCMH.cc \
ProcessCMH.h \
ProcessConcordance.cc \
ProcessConcordance.h \
ProcessDeletions.cc \
ProcessDeletions.h \
ProcessEarth.cc \
ProcessEarth.h \
ProcessEigenstratOutput.cc \
ProcessEigenstratOutput.h \
ProcessEpistasis.cc \
ProcessEpistasis.h \
ProcessFBATOutput.cc \
ProcessFBATOutput.h \
ProcessFilterProcess.cc \
ProcessFilterProcess.h \
ProcessFst.cc \
ProcessFst.h \
ProcessGenderCheck.cc \
ProcessGenderCheck.h \
ProcessGRROutput.cc \
ProcessGRROutput.h \
Process.h \
ProcessHomozygous.cc \
ProcessHomozygous.h \
ProcessHWEquilibrium.cc \
ProcessHWEquilibrium.h \
ProcessIBS.cc \
ProcessIBS.h \
ProcessImputeOutput.cc \
ProcessImputeOutput.h \
ProcessKinship.cc \
ProcessKinship.h \
ProcessLAPISOutput.cc \
ProcessLAPISOutput.h \
ProcessLD.cc \
ProcessLD.h \
ProcessLinearReg.cc \
ProcessLinearReg.h \
ProcessLogReg.cc \
ProcessLogReg.h \
ProcessMarkerGenoEff.cc \
ProcessMarkerGenoEff.h \
ProcessMDR.cc \
ProcessMDR.h \
ProcessMDROutput.cc \
ProcessMDROutput.h \
ProcessMDRPDT.cc \
ProcessMDRPDT.h \
ProcessMendelianErrors.cc \
ProcessMendelianErrors.h \
ProcessMitoCheck.cc \
ProcessMitoCheck.h \
ProcessPDT2Output.cc \
ProcessPDT2Output.h \
ProcessPEDOutput.cc \
ProcessPEDOutput.h \
ProcessPHASEOutput.cc \
ProcessPHASEOutput.h \
ProcessPowerMarkerOutput.cc \
ProcessPowerMarkerOutput.h \
ProcessQTDTOutput.cc \
ProcessQTDTOutput.h \
ProcessRunTDT.cc \
ProcessRunTDT.h \
ProcessSampleGenoEff.cc \
ProcessSampleGenoEff.h \
ProcessSTRUCTOutput.cc \
ProcessSTRUCTOutput.h \
ProcessSuperlinkOutput.cc \
ProcessSuperlinkOutput.h \
ProcessTPEDOutput.cc \
ProcessTPEDOutput.h \
sockets.cpp \
sockets.h \
Step.cc \
Step.h \
Vars.h \
wasp.cc \
wasp.h \
DataSetObject.h



plato_LDADD=method_lib/libmethods.la $(BOOST_THREAD_LIB) $(LIBR_LIBS)
