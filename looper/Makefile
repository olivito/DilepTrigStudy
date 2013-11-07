CC = g++
#INCLUDE = -I../ -I./ -I/opt/root_v5.27.06/include/
#ROOFITINCLUDE = $(shell scramv1 tool info roofitcore | grep INCLUDE | sed 's/^INCLUDE=/-I/')
INCLUDE = -I../ -I./ $(ROOFITINCLUDE)
CFLAGS = -Wall -g -fPIC $(shell root-config --cflags) $(INCLUDE) $(EXTRACFLAGS) -DTOOLSLIB
LINKER = g++

LINKERFLAGS = $(shell root-config --ldflags) $(shell root-config --libs) -lMinuit
ifeq ($(shell root-config --platform),macosx)
	LINKERFLAGS = -dynamiclib -undefined dynamic_lookup -Wl,-x -O -Xlinker -bind_at_load -flat_namespace $(shell root-config --libs) -lEG -lGenVector
endif

SOURCES = dilepStudyLooper.cc 
OBJECTS = $(SOURCES:.cc=.o) LinkDef_out.o
LIB = libdilepStudyLooper.so

CORESOURCES = ../CORE/CMS2.cc ../CORE/utilities.cc ../CORE/ssSelections.cc ../CORE/electronSelections.cc ../CORE/electronSelectionsParameters.cc ../CORE/MITConversionUtilities.cc ../CORE/muonSelections.cc ../CORE/eventSelections.cc ../CORE/trackSelections.cc ../CORE/metSelections.cc ../CORE/jetSelections.cc ../CORE/photonSelections.cc ../CORE/triggerUtils.cc ../CORE/triggerSuperModel.cc ../CORE/mcSelections.cc ../CORE/susySelections.cc ../CORE/mcSUSYkfactor.cc ../CORE/SimpleFakeRate.cc ../Tools/goodrun.cc ../Tools/vtxreweight.cc ../Tools/msugraCrossSection.cc ../Tools/pfjetMVAtools.cc ../CORE/jetsmear/JetSmearer.cc ../CORE/jetsmear/JetResolution.cc ../CORE/jetsmear/SigInputObj.cc ../CORE/jetSmearingTools.cc ../CORE/QuarkGluonTagger/QGLikelihoodCalculator.cc ../CORE/QuarkGluonTagger/QuarkGluonTagger.cc 

COREOBJECTS = $(CORESOURCES:.cc=.o) 
CORELIB = libdilepStudyCORE.so

#FWLIB = libMiniFWLite.so
FWLIB = ../Tools/MiniFWLite/libMiniFWLite.so

LIBS = $(LIB) $(CORELIB) $(FWLIB) 

.PHONY: all help compile clean cms2env

libs:	$(LIBS)

$(LIB):	$(OBJECTS) 
	echo "Linking $(LIB)"; \
	$(LINKER) $(LINKERFLAGS) -shared $(OBJECTS) -o $@

$(CORELIB):	$(COREOBJECTS) 
	echo "Linking $(CORELIB)"; \
	$(LINKER) $(LINKERFLAGS) -shared $(COREOBJECTS) -o $@

$(FWLIB):
	echo "making MiniFWLite"; \
        cd ../Tools/MiniFWLite; \
        $(MAKE) -f Makefile; cd -; \
	cp $(FWLIB) .;

LinkDef_out.cxx: LinkDef.h dilepStudyLooper.h
	rootcint -f $@ -c $(INCLUDE) dilepStudyLooper.h $<

# General rule for making object files
%.d:	%.cc
	$(CC) -MM -MT $@ -MT ${@:.d=.o} $(CFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@
%.d:	%.cxx
	$(CC) -MM -MT $@ -MT ${@:.d=.o} $(CFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@

%.o: 	%.cc 
	$(CC) $(CFLAGS) $< -c -o $@

%.o: 	%.cxx
	$(CC) $(CFLAGS) $< -c -o $@

.PHONY: all
all:	$(LIBS)  

.PHONY: clean
clean:  
	rm -v -f \
        ../CORE/*.o \
        ../CORE/*.d \
        ../CORE/*.so \
        ../Tools/*.o \
	../Tools/*.d \
        ../Tools/MiniFWLite/*.o \
	../Tools/MiniFWLite/*.d \
        $(CORELIB) $(LIB) $(FWLIB) \
        LinkDef_out.cxx  LinkDef_out.h \
        *.o *.d *.so; echo "Done"

test: 
	@echo $(INCLUDE)

-include $(CORESOURCES:.cc=.d) $(CORE:.cc=.d)
-include $(LIBDIR)/LinkDef_out.d

