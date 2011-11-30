# The list of executables we will compile.
# #
PROGRAMS=methods plato
LIBPLATODIR=lib/
LIBPLATO=$(LIBPLATODIR)libplato.a
PLATO_AS_LIB=sqllib methods $(LIBPLATO)
METHODDIR=method_lib
LIBDIR=-Llib -L/opt/local/lib #-L/home/cozartc/boost/stage/lib -L/home/cozartc/sqlitewrapped/lib
LIB=-lm -lmethods -lboost_thread-mt#mgw44-mt-1_43 -lsqlite3 -lsqlitewrapped#-lreadline -lintl -lglib-2.0
INCLUDEDIR=-I. -I$(METHODDIR) -I/opt/local/include #-I/home/cozartc/boost -I/home/cozartc/sqlitewrapped/lib#-I/usr/local/include
SYS=MAC
#DB=USE_DB
COMPASLIB=PLATOLIB
#R=USE_R

#
# # Default target builds all programs: (make)
# #
#all: $(PROGRAMS)
#
# # Cleanup the current directory: (make clean)
# #
#clean:
#	rm -f $(PROGRAMS) *.o
#
#
#   # Include the Oracle compiler settings.  Typically overkill, but
# platform-independent.
# #
OTHEROPTS=
#BIT=-m32
COMP=g++ #/opt/local/bin/g++-mp-4.4
CC=$(COMP) $(BIT) -O3  -Wall -Wno-deprecated -g -D_FILE_OFFSET_BITS=64 $(INCLUDEDIR) $(OTHEROPTS)
#CC=g++ -O3 -D_FILE_OFFSET_BITS=64 -mno-cygwin
#CC=mpicxx -m64

ifeq ($(SYS),UNIX)
  CC += -DUNIX -static
endif
ifeq ($(SYS),WIN)
  CC += -DWIN -static
  LIB += -lwsock32
endif
ifeq ($(SYS),MAC)
  CC += -DMAC
endif

OBJECTS = ProcessKinship.o ProcessFst.o Step.o Process.o Percent.o Chrom.o ProcessMarkerGenoEff.o ProcessSampleGenoEff.o PercentByFamily.o ProcessAlleleFrequency.o \
		  ProcessMendelianErrors.o ProcessHWEquilibrium.o ProcessGenderCheck.o ProcessRunTDT.o ProcessGRROutput.o dcdflib.o ProcessPEDOutput.o ProcessBINOutput.o \
		  ProcessCaConChisq.o ipmpar.o ProcessSTRUCTOutput.o ProcessPHASEOutput.o ProcessEigenstratOutput.o \
		  ProcessBEAGLEOutput.o ProcessLAPISOutput.o ProcessMDROutput.o ProcessHomozygous.o ProcessLD.o Finalize.o ProcessPowerMarkerOutput.o \
		  ExampleModule.o ProcessDeletions.o ProcessMitoCheck.o ProcessFBATOutput.o ProcessQTDTOutput.o ProcessPDT2Output.o ProcessConcordance.o \
		  sockets.o ProcessSuperlinkOutput.o ProcessTPEDOutput.o ProcessLogReg.o ProcessCMH.o ProcessLinearReg.o ProcessIBS.o ProcessFilterProcess.o \
		  ProcessMDR.o ProcessClusterMissing.o ProcessMDRPDT.o ProcessEpistasis.o ProcessImputeOutput.o

INTERFACES = ProcessKinship.h ProcessFst.h Step.h Process.h Percent.h Chrom.h ProcessMarkerGenoEff.h ProcessSampleGenoEff.h PercentByFamily.h ProcessAlleleFrequency.h \
		  ProcessMendelianErrors.h ProcessHWEquilibrium.h ProcessGenderCheck.h ProcessRunTDT.h ProcessGRROutput.h dcdflib.h ProcessPEDOutput.h ProcessBINOutput.h \
		  ProcessCaConChisq.h ipmpar.h ProcessSTRUCTOutput.h ProcessPHASEOutput.h ProcessEigenstratOutput.h \
		  ProcessBEAGLEOutput.h ProcessLAPISOutput.h ProcessMDROutput.h ProcessHomozygous.h ProcessLD.h Finalize.h ProcessPowerMarkerOutput.h \
		  ExampleModule.h ProcessDeletions.h ProcessMitoCheck.h ProcessFBATOutput.h ProcessQTDTOutput.h ProcessPDT2Output.h ProcessConcordance.h \
		  sockets.h ProcessSuperlinkOutput.h ProcessTPEDOutput.h ProcessLogReg.h ProcessCMH.h ProcessLinearReg.h ProcessIBS.h ProcessFilterProcess.h \
		  ProcessMDR.h ProcessClusterMissing.h ProcessMDRPDT.h ProcessEpistasis.h ProcessImputeOutput.h Controller.h Vars.h

LIB_OBJECTS := Controller.o Vars.o $(OBJECTS)

OBJECTS += wasp.o

ifeq ($(R),USE_R)
	CC += -DUSE_R
	OBJECTS += ProcessEarth.o
	LIB += -lR
endif
ifeq ($(COMPASLIB), PLATOLIB)
	CC += -DPLATOLIB
	CC += -DUSE_DB
	CC += -DNOSYS
	LIB += -lsqlitewrapped
	INCLUDEDIR := $(INCLUDEDIR) -I./sqlitewrapped
endif

#
# # Generic .cc -> .o production.
# #
#.SUFFIXES: .o .cc
#.cc.o:
#	g++296 $(CFLAGSCC) -c $<

PRELINK = 
AR = ar rv

all: $(PROGRAMS)

clean:
	rm *.o plato wasp; cd method_lib; make clean;

clean_lib:
	rm *.o plato 

sqllib:
	cd sqlitewrapped; make

methods:
	cd method_lib; make

plato:	$(OBJECTS)
	$(CC) -o plato $(OBJECTS) $(LIB) $(SHARED_OCCILIBS) $(CFLAGSCC) $(LIBDIR) 

plato_lib: $(PLATO_AS_LIB)

$(LIBPLATO):	$(LIB_OBJECTS)
	$(CC) $(LIB_OBJECTS) -c $(LIB) $(SHARED_OCCILIBS) $(CFLAGSCC) $(LIBDIR)
	$(PRELINK)
	$(AR) $(LIBPLATO) $?
	@echo $(LIBPLATO) is now up-to-date

wasp.o: wasp.cc wasp.h
	$(CC) wasp.cc -c $(CFLAGSCC) $(LIBDIR) 

ProcessEarth.o: ProcessEarth.cc ProcessEarth.h
	$(CC) ProcessEarth.cc -c $(CFLAGSCC) $(LIBDIR) 

ProcessImputeOutput.o: ProcessImputeOutput.cc ProcessImputeOutput.h
	$(CC) ProcessImputeOutput.cc -c $(CFLAGSCC) $(LIBDIR) 

ProcessKinship.o: ProcessKinship.cc ProcessKinship.h
	$(CC) ProcessKinship.cc -c $(CFLAGSCC) $(LIBDIR) 

ProcessFst.o: ProcessFst.cc ProcessFst.h
	$(CC) ProcessFst.cc -c $(CFLAGSCC) $(LIBDIR) 

sockets.o: sockets.cpp sockets.h sisocks.h
	$(CC) sockets.cpp -c $(CFLAGSCC) $(LIBDIR)

ProcessMendelianErrors.o: ProcessMendelianErrors.cc ProcessMendelianErrors.h
	$(CC) ProcessMendelianErrors.cc -c $(CFLAGSCC) $(LIBDIR)

#Family.o: Family.cc Family.h
#	$(CC) Family.cc -c $(CFLAGSCC) $(LIBDIR) 

#Sample.o: Sample.cc Sample.h
#	$(CC) Sample.cc -c $(CFLAGSCC) $(LIBDIR)

#Marker.o: Marker.cc Marker.h
#	$(CC) Marker.cc -c $(CFLAGSCC) $(LIBDIR)

Process.o: Process.cc Process.h
	$(CC) Process.cc -c $(CFLAGSCC) $(LIBDIR)

Step.o: Step.cc Step.h
	$(CC) Step.cc -c $(CFLAGSCC) $(LIBDIR)

Chrom.o: Chrom.cc Chrom.h
	$(CC) Chrom.cc -c $(CFLAGSCC) $(LIBDIR)

Percent.o: Percent.cc Percent.h
	$(CC) Percent.cc -c $(CFLAGSCC) $(LIBDIR)

ProcessFilterProcess.o: ProcessFilterProcess.cc ProcessFilterProcess.h
	$(CC) ProcessFilterProcess.cc -c $(CFLAGSCC) $(LIBDIR)

#Options.o: Options.cc Options.h
#	$(CC) Options.cc -c $(CFLAGSCC) $(LIBDIR)

#General.o: General.cc General.h
#	$(CC) General.cc -c $(CFLAGSCC) $(LIBDIR)

ProcessMarkerGenoEff.o: ProcessMarkerGenoEff.cc ProcessMarkerGenoEff.h
	$(CC) ProcessMarkerGenoEff.cc -c $(CFLAGSCC) $(LIBDIR)

ProcessAlleleFrequency.o: ProcessAlleleFrequency.cc ProcessAlleleFrequency.h
	$(CC) ProcessAlleleFrequency.cc -c $(CFLAGSCC) $(LIBDIR)

ProcessMDR.o: ProcessMDR.cc ProcessMDR.h
	$(CC) ProcessMDR.cc -c $(CFLAGSCC) $(LIBDIR)

ProcessMDRPDT.o: ProcessMDRPDT.cc ProcessMDRPDT.h
	$(CC) ProcessMDRPDT.cc -c $(CFLAGSCC) $(LIBDIR)

ProcessClusterMissing.o: ProcessClusterMissing.cc ProcessClusterMissing.h
	$(CC) ProcessClusterMissing.cc -c $(CFLAGSCC) $(LIBDIR)

ProcessSampleGenoEff.o: ProcessSampleGenoEff.cc ProcessSampleGenoEff.h
	$(CC) ProcessSampleGenoEff.cc -c $(CFLAGSCC) $(LIBDIR)

ProcessHWEquilibrium.o: ProcessHWEquilibrium.cc ProcessHWEquilibrium.h
	$(CC) ProcessHWEquilibrium.cc -c $(CFLAGSCC) $(LIBDIR)

ProcessGenderCheck.o: ProcessGenderCheck.cc ProcessGenderCheck.h
	$(CC) ProcessGenderCheck.cc -c $(CFLAGSCC) $(LIBDIR)

ProcessRunTDT.o: ProcessRunTDT.cc ProcessRunTDT.h
	$(CC) ProcessRunTDT.cc -c $(CFLAGSCC) $(LIBDIR)
ProcessCMH.o: ProcessCMH.cc ProcessCMH.h
	$(CC) ProcessCMH.cc -c $(CFLAGSCC) $(LIBDIR)

#QualityScore.o: QualityScore.cc QualityScore.h
#	$(CC) QualityScore.cc -c $(CFLAGSCC) $(LIBDIR)

PercentByFamily.o: PercentByFamily.cc PercentByFamily.h
	$(CC) PercentByFamily.cc -c $(CFLAGSCC) $(LIBDIR)

ProcessGRROutput.o: ProcessGRROutput.cc ProcessGRROutput.h
	$(CC) ProcessGRROutput.cc -c $(CFLAGSCC) $(LIBDIR)
ProcessIBS.o: ProcessIBS.cc ProcessIBS.h
	$(CC) ProcessIBS.cc -c $(CFLAGSCC) $(LIBDIR)

ProcessSTRUCTOutput.o: ProcessSTRUCTOutput.cc ProcessSTRUCTOutput.h
	$(CC) ProcessSTRUCTOutput.cc -c $(CFLAGSCC) $(LIBDIR)

ProcessEigenstratOutput.o: ProcessEigenstratOutput.cc ProcessEigenstratOutput.h
	$(CC) ProcessEigenstratOutput.cc -c $(CFLAGSCC) $(LIBDIR)
ProcessPEDOutput.o: ProcessPEDOutput.cc ProcessPEDOutput.h
	$(CC) ProcessPEDOutput.cc -c $(CFLAGSCC) $(LIBDIR)
ProcessBINOutput.o: ProcessBINOutput.cc ProcessBINOutput.h
	$(CC) ProcessBINOutput.cc -c $(CFLAGSCC) $(LIBDIR)
ProcessEpistasis.o: ProcessEpistasis.cc ProcessEpistasis.h
	$(CC) ProcessEpistasis.cc -c $(CFLAGSCC) $(LIBDIR)
ProcessTPEDOutput.o: ProcessTPEDOutput.cc ProcessTPEDOutput.h
	$(CC) ProcessTPEDOutput.cc -c $(CFLAGSCC) $(LIBDIR)
#QSOutput.o: QSOutput.cc QSOutput.h
#	$(CC) QSOutput.cc -c $(CFLAGSCC) $(LIBDIR)
ProcessPHASEOutput.o: ProcessPHASEOutput.cc ProcessPHASEOutput.h
	$(CC) ProcessPHASEOutput.cc -c $(CFLAGSCC) $(LIBDIR)
#PartialOutput.o: PartialOutput.cc PartialOutput.h
#	$(CC) PartialOutput.cc -c $(CFLAGSCC) $(LIBDIR)
#ChiSquareAllelic.o: ChiSquareAllelic.cpp ChiSquareAllelic.h
#	$(CC) ChiSquareAllelic.cpp -c $(CFLAGSCC) $(LIBDIR)
#ChiSquareArmitage.o: ChiSquareArmitage.cpp ChiSquareArmitage.h
#	$(CC) ChiSquareArmitage.cpp -c $(CFLAGSCC) $(LIBDIR)
ProcessCaConChisq.o: ProcessCaConChisq.cc ProcessCaConChisq.h
	$(CC) ProcessCaConChisq.cc -c $(CFLAGSCC) $(LIBDIR)
#FisherExact.o: FisherExact.cpp FisherExact.h 
#	$(CC) FisherExact.cpp -c $(CFLAGSCC) $(LIBDIR)
#StepOptions.o: StepOptions.cc StepOptions.h
#	$(CC) StepOptions.cc -c $(CFLAGSCC) $(LIBDIR)
ProcessBEAGLEOutput.o: ProcessBEAGLEOutput.cc ProcessBEAGLEOutput.h
	$(CC) ProcessBEAGLEOutput.cc -c $(CFLAGSCC) $(LIBDIR)
ProcessLAPISOutput.o: ProcessLAPISOutput.cc ProcessLAPISOutput.h
	$(CC) ProcessLAPISOutput.cc -c $(CFLAGSCC) $(LIBDIR)
ProcessMDROutput.o: ProcessMDROutput.cc ProcessMDROutput.h
	$(CC) ProcessMDROutput.cc -c $(CFLAGSCC) $(LIBDIR)
ProcessHomozygous.o: ProcessHomozygous.cc ProcessHomozygous.h
	$(CC) ProcessHomozygous.cc -c $(CFLAGSCC) $(LIBDIR)
ProcessLD.o: ProcessLD.cc ProcessLD.h
	$(CC) ProcessLD.cc -c $(CFLAGSCC) $(LIBDIR)
ProcessLogReg.o: ProcessLogReg.cc ProcessLogReg.h
	$(CC) ProcessLogReg.cc -c $(CFLAGSCC) $(LIBDIR)
ProcessLinearReg.o: ProcessLinearReg.cc ProcessLinearReg.h
	$(CC) ProcessLinearReg.cc -c $(CFLAGSCC) $(LIBDIR)
Finalize.o: Finalize.cc Finalize.h
	$(CC) Finalize.cc -c $(CFLAGSCC) $(LIBDIR)
ProcessPowerMarkerOutput.o: ProcessPowerMarkerOutput.cc ProcessPowerMarkerOutput.h
	$(CC) ProcessPowerMarkerOutput.cc -c $(CFLAGSCC) $(LIBDIR)
ipmpar.o: ipmpar.c max.h
	$(CC) ipmpar.c -c $(CFLAGSCC) $(LIBDIR)
ProcessDeletions.o: ProcessDeletions.cc ProcessDeletions.h
	$(CC) ProcessDeletions.cc -c $(CFLAGSCC) $(LIBDIR)
ProcessMitoCheck.o: ProcessMitoCheck.cc ProcessMitoCheck.h
	$(CC) ProcessMitoCheck.cc -c $(CFLAGSCC) $(LIBDIR)
ProcessFBATOutput.o: ProcessFBATOutput.cc ProcessFBATOutput.h
	$(CC) ProcessFBATOutput.cc -c $(CFLAGSCC) $(LIBDIR)
ProcessQTDTOutput.o: ProcessQTDTOutput.cc ProcessQTDTOutput.h
	$(CC) ProcessQTDTOutput.cc -c $(CFLAGSCC) $(LIBDIR)
ExampleModule.o: ExampleModule.cc ExampleModule.h
	$(CC) ExampleModule.cc -c $(CFLAGSCC) $(LIBDIR)
ProcessPDT2Output.o: ProcessPDT2Output.cc ProcessPDT2Output.h
	$(CC) ProcessPDT2Output.cc -c $(CFLAGSCC) $(LIBDIR)
ProcessSuperlinkOutput.o: ProcessSuperlinkOutput.cc ProcessSuperlinkOutput.h
	$(CC) ProcessSuperlinkOutput.cc -c $(CFLAGSCC) $(LIBDIR)
ProcessConcordance.o: ProcessConcordance.cc ProcessConcordance.h
	$(CC) ProcessConcordance.cc -c $(CFLAGSCC) $(LIBDIR)
Controller.o: Controller.cpp Controller.h
	$(CC) Controller.cpp -c $(CFLAGSCC) $(LIBDIR)
Vars.o: Vars.cpp Vars.h
	$(CC) Vars.cpp -c $(CFLAGSCC) $(LIBDIR)
