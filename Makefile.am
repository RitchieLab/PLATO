ACLOCAL_AMFLAGS= -I ./make/m4

SUBDIRS=method_lib

EXTRA_DIST=CREATED CHANGELOG README.txt

CLEANFILES=CREATED CHANGELOG README.txt

dist-hook:
	chmod u+w $(distdir)/CREATED
	echo "Distribution date: `date`" > CREATED

# This is here to prevent installation of LOKI during a distcheck
DISTCHECK_CONFIGURE_FLAGS=--disable-loki

AM_CPPFLAGS=$(BOOST_CPPFLAGS) -I$(top_srcdir)/method_lib 
AM_LDFLAGS=$(BOOST_LDFLAGS) 

bin_PROGRAMS=plato

plato_SOURCES= Process.cpp \
Process.h \
main.cpp \
main.h \
ProcessFactory.cpp \
ProcessFactory.h \
InputManager.h \
InputManager.cpp

plato_LDADD=method_lib/libmethods.la $(BOOST_THREAD_LIB) $(BOOST_PROGRAM_OPTIONS_LIB)


orig_plato_src = cdflib.h \
dcdflib.cpp \
ipmpar.c \
Finalize.cc \
Finalize.h \
process/PercentByFamily.cc \
process/PercentByFamily.h \
process/ProcessAlleleFrequency.cc \
process/ProcessAlleleFrequency.h \
process/ProcessBEAGLEOutput.cc \
process/ProcessBEAGLEOutput.h \
process/ProcessBINOutput.cc \
process/ProcessBINOutput.h \
process/ProcessCaConChisq.cc \
process/ProcessCaConChisq.h \
Process.cc \
process/ProcessClusterMissing.cc \
process/ProcessClusterMissing.h \
process/ProcessCMH.cc \
process/ProcessCMH.h \
process/ProcessConcordance.cc \
process/ProcessConcordance.h \
process/ProcessDeletions.cc \
process/ProcessDeletions.h \
process/ProcessEigenstratOutput.cc \
process/ProcessEigenstratOutput.h \
process/ProcessEpistasis.cc \
process/ProcessEpistasis.h \
process/ProcessFBATOutput.cc \
process/ProcessFBATOutput.h \
process/ProcessFilterProcess.cc \
process/ProcessFilterProcess.h \
process/ProcessFst.cc \
process/ProcessFst.h \
process/ProcessGenderCheck.cc \
process/ProcessGenderCheck.h \
process/ProcessGRROutput.cc \
process/ProcessGRROutput.h \
Process.h \
process/ProcessHomozygous.cc \
process/ProcessHomozygous.h \
process/ProcessHWEquilibrium.cc \
process/ProcessHWEquilibrium.h \
process/ProcessIBS.cc \
process/ProcessIBS.h \
process/ProcessImputeOutput.cc \
process/ProcessImputeOutput.h \
process/ProcessKinship.cc \
process/ProcessKinship.h \
process/ProcessLAPISOutput.cc \
process/ProcessLAPISOutput.h \
process/ProcessLD.cc \
process/ProcessLD.h \
process/ProcessLinearReg.cc \
process/ProcessLinearReg.h \
process/ProcessLogReg.cc \
process/ProcessLogReg.h \
process/ProcessMarkerGenoEff.cc \
process/ProcessMarkerGenoEff.h \
process/ProcessMDR.cc \
process/ProcessMDR.h \
process/ProcessMDROutput.cc \
process/ProcessMDROutput.h \
process/ProcessMDRPDT.cc \
process/ProcessMDRPDT.h \
process/ProcessMendelianErrors.cc \
process/ProcessMendelianErrors.h \
process/ProcessMitoCheck.cc \
process/ProcessMitoCheck.h \
process/ProcessPDT2Output.cc \
process/ProcessPDT2Output.h \
process/ProcessPEDOutput.cc \
process/ProcessPEDOutput.h \
process/ProcessPHASEOutput.cc \
process/ProcessPHASEOutput.h \
process/ProcessPowerMarkerOutput.cc \
process/ProcessPowerMarkerOutput.h \
process/ProcessQTDTOutput.cc \
process/ProcessQTDTOutput.h \
process/ProcessRunTDT.cc \
process/ProcessRunTDT.h \
process/ProcessSampleGenoEff.cc \
process/ProcessSampleGenoEff.h \
process/ProcessSTRUCTOutput.cc \
process/ProcessSTRUCTOutput.h \
process/ProcessSuperlinkOutput.cc \
process/ProcessSuperlinkOutput.h \
process/ProcessTPEDOutput.cc \
process/ProcessTPEDOutput.h \
wasp.cc \
wasp.h \
process/ProcessInteraction.h \
process/ProcessInteraction.cpp \
ProcessFactory.cpp \
ProcessFactory.h

