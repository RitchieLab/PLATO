# Makefile for datasim library
CC = g++
LIBDIR = ../lib/
OS = MAC
RLIB = -lR
RLIBS = $(LIBDIR)libR.a #libR.dylib $(LIBDIR)libRblas.dylib
LIB= $(LIBDIR)libmethods.a
CFLAGS=  -Wall -Wno-deprecated
CFLAGS =  -Wall -Wno-deprecated
CFLAGSCC=  -Wall -Wno-deprecated
MDRPDTDIR = mdrpdt
FILTDIR = filters
MYSQLDEF = 

ifeq ($(OS),UNIX)
	CFLAGS += -DUNIX
	CFLAGSCC += -DUNIX
endif
ifeq ($(OS),WIN)
	CFLAGS += -DWIN
	CFLAGSCC += -DWIN
endif
ifeq ($(OS),MAC)
	CFLAGS += -g -DMAC
	CFLAGSCC += -g -DMAC
endif

#CFLAGS= -g -Wall -Wno-deprecated
#CFLAGSCC = -g

OBJS = Kinship.o earth.o LogisticRegression.o Marker.o Sample.o ComboGenerator.o \
	DataSet.o FlatIndex.o Family.o General.o Options.o \
	StepOptions.o AlleleFrequency.o MarkerGenoEff.o HWEquilibrium.o CaConChisq.o \
	FisherExact.o ChiSquareArmitage.o dcdflib.o ipmpar.o ChiSquareAllelic.o \
	CondLogRegression.o conditional_lr.o RunTDT.o MendelianErrors.o \
	MitoCheck.o MDR.o ContingencyTable.o UncertaintyCoefficient.o \
	LikelihoodRatio.o OddsRatio.o NMI.o Concordance.o Deletions.o \
	Homozygous.o GenderCheck.o ChiSquare.o SampleGenoEff.o \
	BEAGLEOutput.o FBATOutput.o GRROutput.o LAPISOutput.o MDROutput.o \
	PEDOutput.o PDT2Output.o PHASEOutput.o PowerMarkerOutput.o \
	QTDTOutput.o STRUCTOutput.o SuperlinkOutput.o LD.o InputFilter.o TPEDOutput.o \
	CMH.o Model.o LinearRegression.o MultComparison.o LogTransform.o \
	SquareRootTransform.o BoxCoxTransform.o EigenstratOutput.o IBS.o \
	CovTrait.o Covariate.o Trait.o ClusterMissing.o Fst.o BINOutput.o

INTERFACES =  Kinship.h earth.h LogisticRegression.h Marker.h Sample.h ComboGenerator.h \
	DataSet.h FlatIndex.h Family.h General.h Globals.h Options.h \
	Helper.h StepOptions.h AlleleFrequency.h MarkerGenoEff.h HWEquilibrium.h CaConChisq.h cdflib.h \
	FisherExact.h ChiSquareArmitage.h max.h ChiSquareAllelic.h CondLogRegression.h \
	conditional_lr.h RunTDT.h MendelianErrors.h MitoCheck.h MDR.h ContingencyTable.h \
	UncertaintyCoefficient.h  LikelihoodRatio.h OddsRatio.h NMI.h Concordance.h \
	Deletions.h Homozygous.h GenderCheck.h ChiSquare.h SampleGenoEff.h \
	BEAGLEOutput.h FBATOutput.h GRROutput.h LAPISOutput.h MDROutput.h \
	PEDOutput.h PDT2Output.h PHASEOutput.h PowerMarkerOutput.h \
	QTDTOutput.h STRUCTOutput.h SuperlinkOutput.h LD.h InputFilter.h TPEDOutput.h \
	CMH.h Model.h LinearRegression.h MultComparison.h LogTransform.h \
	SquareRootTransform.h BoxCoxTransform.h EigenstratOutput.h IBS.h \
	CovTrait.h Covariate.h Trait.h Fst.h BINOutput.h

MDRPDTINTERFACES = $(INTERFACES) $(MDRPDTDIR)/evaluationmethod.h \
	$(MDRPDTDIR)/familyrepository.h $(MDRPDTDIR)/foldproduction.h $(MDRPDTDIR)/foldstatistic.h \
	$(MDRPDTDIR)/genotypecounter.h $(MDRPDTDIR)/genotyperepository.h $(MDRPDTDIR)/individual.h \
	$(MDRPDTDIR)/lineparser.h $(MDRPDTDIR)/matchedoddsratio.h $(MDRPDTDIR)/MDRPDTRepository.h \
	$(MDRPDTDIR)/orderreport.h $(MDRPDTDIR)/pdtmodel.h $(MDRPDTDIR)/pedigree.h \
	$(MDRPDTDIR)/ptestdistribution.h $(MDRPDTDIR)/random.h $(MDRPDTDIR)/rbtree.h \
	$(MDRPDTDIR)/sibship.h $(MDRPDTDIR)/strings.h $(MDRPDTDIR)/tfinalreport.h \
	$(MDRPDTDIR)/tstatistic.h $(MDRPDTDIR)/types.h $(MDRPDTDIR)/utility.h $(MDRPDTDIR)/randomc.h

MDRPDTOBJS = $(MDRPDTDIR)/evaluationmethod.o \
	$(MDRPDTDIR)/familyrepository.o $(MDRPDTDIR)/foldproduction.o $(MDRPDTDIR)/foldstatistic.o \
	$(MDRPDTDIR)/genotyperepository.o $(MDRPDTDIR)/individual.o $(MDRPDTDIR)/mersenne.o \
	$(MDRPDTDIR)/orderreport.o $(MDRPDTDIR)/pdtmodel.o $(MDRPDTDIR)/pedigree.o \
	$(MDRPDTDIR)/ptestdistribution.o $(MDRPDTDIR)/random.o \
	$(MDRPDTDIR)/sibship.o $(MDRPDTDIR)/strings.o $(MDRPDTDIR)/tfinalreport.o \
	$(MDRPDTDIR)/tstatistic.o  $(MDRPDTDIR)/lineparser.o $(MDRPDTDIR)/matchedoddsratio.o \
	$(MDRPDTDIR)/MDRPDTRepository.o MDRPDT.o

FILTERINTERFACES = $(INTERFACES) $(FILTDIR)/AlleleFrequencyFilter.h $(FILTDIR)/Analysis.h \
	$(FILTDIR)/MarsFilter.h $(FILTDIR) \
	$(FILTDIR)/AnalysisSet.h $(FILTDIR)/ConfigReader.h $(FILTDIR)/DataDefinitions.h \
	$(FILTDIR)/ArmitageFilter.h  $(FILTDIR)/BioBinaryReader.h  \
	$(FILTDIR)/BioFilter.h $(FILTDIR)/ChiSquareFilter.h \
	$(FILTDIR)/ConditionalLRFilter.h $(FILTDIR)/ConfigFileReader.h \
	$(FILTDIR)/ContingencyFilter.h \
	$(FILTDIR)/DataExcept.h $(FILTDIR)/FileResults.h  $(FILTDIR)/Filter.h \
	$(FILTDIR)/FilterExcept.h $(FILTDIR)/LikelihoodRatioFilter.h $(FILTDIR)/LogInfo.h \
	$(FILTDIR)/LogRegFilter.h $(FILTDIR)/MarkerGenoFilter.h $(FILTDIR)/MaxScoreFilter.h \
	$(FILTDIR)/MDRFilter.h $(FILTDIR)/MDRPDTFilter.h $(FILTDIR)/OddsRatioFilter.h \
	$(FILTDIR)/OutputWriterFactory.h $(FILTDIR)/PermScoreFilter.h $(FILTDIR)/Permutation.h \
	$(FILTDIR)/ResultsManager.h  $(FILTDIR)/Stringmanip.h \
	$(FILTDIR)/UncertaintyCoeffFilter.h $(FILTDIR)/Filters.h \
	$(FILTDIR)/PlatoDefinitions.h $(FILTDIR)/PlatoExcept.h $(FILTDIR)/ResultSet.h \
	$(FILTDIR)/FilterFactory.h $(FILTDIR)/NMIFilter.h $(FILTDIR)/SnpFileReader.h \
	$(FILTDIR)/BinSnpFileReader.h $(FILTDIR)/TextSnpFileReader.h

FILTEROBJS = $(FILTDIR)/AlleleFrequencyFilter.o $(FILTDIR)/Analysis.o $(FILTDIR)/AnalysisSet.o \
	$(FILTDIR)/MarsFilter.o \
	$(FILTDIR)/ArmitageFilter.o  $(FILTDIR)/BioBinaryReader.o \
	$(FILTDIR)/BioFilter.o $(FILTDIR)/ChiSquareFilter.o \
	$(FILTDIR)/ConditionalLRFilter.o $(FILTDIR)/ConfigFileReader.o  $(FILTDIR)/ContingencyFilter.o \
	$(FILTDIR)/DataExcept.o $(FILTDIR)/FileResults.o  $(FILTDIR)/Filter.o \
	$(FILTDIR)/FilterExcept.o $(FILTDIR)/LikelihoodRatioFilter.o $(FILTDIR)/LogInfo.o \
	$(FILTDIR)/LogRegFilter.o $(FILTDIR)/MarkerGenoFilter.o $(FILTDIR)/MaxScoreFilter.o \
	$(FILTDIR)/MDRFilter.o $(FILTDIR)/MDRPDTFilter.o $(FILTDIR)/OddsRatioFilter.o \
	$(FILTDIR)/OutputWriterFactory.o $(FILTDIR)/PermScoreFilter.o $(FILTDIR)/Permutation.o \
	$(FILTDIR)/ResultsManager.o  $(FILTDIR)/Stringmanip.o \
	$(FILTDIR)/UncertaintyCoeffFilter.o FilterProcess.o \
	$(FILTDIR)/FilterFactory.o $(FILTDIR)/NMIFilter.o $(FILTDIR)/BinSnpFileReader.o \
	$(FILTDIR)/TextSnpFileReader.o

ifeq ($(MYSQLDEF),-DMYSQLDEFINED)
FILTEROBJS +=  $(FILTDIR)/DB.o $(FILTDIR)/BioRecordResults.o $(FILTDIR)/BioFilterCombinations.o
FILTERINTERFACES += $(FILTDIR)/DB.h $(FILTDIR)/BioRecordResults.h $(FILTDIR)/BioFilterCombinations.h
endif

# GNU compiler flags (works with g++ and egcs)
#BIT		    = -m32
CC          = g++ -O2 $(BIT) #/opt/local/bin/g++-mp-4.4 -O2
PRELINK     =
LD          = g++ -w  $(BIT) #/opt/local/bin/g++-mp-4.4 -w
AR          = ar rv
INSTALL     = install -c
RM          = rm
INC_DIRS    = -I. 

all: $(LIB)

MDRPDT.o: MDRPDT.cpp $(MDRPDTINTERFACES) $(INTERFACES)
	$(CC) $(CFLAGS) $(INC_DIRS) -c $< -o $@
FilterProcess.o: FilterProcess.cpp $(INTERFACES) $(FILTERINTERFACES)
	$(CC) $(CFLAGS) $(INC_DIRS) -c $< -o $@
earth.o: earth.cc 
	$(CC) $(CFLAGS) $(INC_DIRS) -c $< -o $@
%.o: %.cpp $(INTERFACES)
	$(CC) $(CFLAGS) $(INC_DIRS) -c $< -o $@
%.o: %.cc $(INTERFACES)
	$(CC) $(CFLAGSCC) $(INC_DIRS) -c $< -o $@
$(LIB): $(OBJS) $(MDRPDTOBJS) $(FILTEROBJS)
	$(PRELINK)
	$(AR) $(LIB) $?
	@echo $(LIB) is now up-to-date
.cpp.a:;
clean:
	rm -f *.o $(MDRPDTOBJS) $(FILTEROBJS) *~ core $(LIB) 
