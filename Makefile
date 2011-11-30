# The list of executables we will compile.
# #
PROGRAMS=plato
METHODDIR=../../development/libraries/method_lib
LIBDIR=-L../../development/libraries/lib
LIB=-lm -lmethods
INCLUDEDIR=-I. -I$(METHODDIR)
SYS=MAC

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
CC=g++ -O3 -g -D_FILE_OFFSET_BITS=64 $(INCLUDEDIR)
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

OBJECTS = Step.o wasp.o Process.o Percent.o Chrom.o ProcessMarkerGenoEff.o ProcessSampleGenoEff.o PercentByFamily.o ProcessAlleleFrequency.o \
		  ProcessMendelianErrors.o ProcessHWEquilibrium.o ProcessGenderCheck.o ProcessRunTDT.o ProcessGRROutput.o dcdflib.o ProcessPEDOutput.o \
		  ProcessCaConChisq.o ipmpar.o ProcessSTRUCTOutput.o ProcessPHASEOutput.o \
		  ProcessBEAGLEOutput.o ProcessLAPISOutput.o ProcessMDROutput.o ProcessHomozygous.o ProcessLD.o Finalize.o ProcessPowerMarkerOutput.o \
		  ExampleModule.o ProcessDeletions.o ProcessMitoCheck.o ProcessFBATOutput.o ProcessQTDTOutput.o ProcessPDT2Output.o ProcessConcordance.o \
		  sockets.o ProcessSuperlinkOutput.o ProcessTPEDOutput.o ProcessLogReg.o ProcessCMH.o

#
# # Generic .cc -> .o production.
# #
#.SUFFIXES: .o .cc
#.cc.o:
#	g++296 $(CFLAGSCC) -c $<

all: $(PROGRAMS)

clean:
	rm *.o plato

plato:	$(OBJECTS)
	$(CC) -o plato $(OBJECTS) $(LIB) $(SHARED_OCCILIBS) $(CFLAGSCC) $(LIBDIR)

wasp.o: wasp.cc wasp.h
	$(CC) wasp.cc -c $(CFLAGSCC) $(LIBDIR) 

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

#Options.o: Options.cc Options.h
#	$(CC) Options.cc -c $(CFLAGSCC) $(LIBDIR)

#General.o: General.cc General.h
#	$(CC) General.cc -c $(CFLAGSCC) $(LIBDIR)

ProcessMarkerGenoEff.o: ProcessMarkerGenoEff.cc ProcessMarkerGenoEff.h
	$(CC) ProcessMarkerGenoEff.cc -c $(CFLAGSCC) $(LIBDIR)

ProcessAlleleFrequency.o: ProcessAlleleFrequency.cc ProcessAlleleFrequency.h
	$(CC) ProcessAlleleFrequency.cc -c $(CFLAGSCC) $(LIBDIR)

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

ProcessSTRUCTOutput.o: ProcessSTRUCTOutput.cc ProcessSTRUCTOutput.h
	$(CC) ProcessSTRUCTOutput.cc -c $(CFLAGSCC) $(LIBDIR)

ProcessPEDOutput.o: ProcessPEDOutput.cc ProcessPEDOutput.h
	$(CC) ProcessPEDOutput.cc -c $(CFLAGSCC) $(LIBDIR)
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
	
