##################################################################################
#
# General Makefile - IMPA - daniel@impa.br
#
##################################################################################
#
#	User customized definitions
#
##################################################################################

#
#	Define executable and/or library name
PROGNAME = rpn
LIBNAME =

#
#	Compilers and linker
CC = gcc
CXX = g++
FC = gcc
LD = $(CXX)

#
#	Include, C-preprocessor, C, C++ and Fortran compilation flags
INCFLAGS = -I ../../../../wave/ode -I ../../../../rpnumerics/methods/Coincidence -I ../../../../wave/util/MultiArray
CPPFLAGS =
CXXFLAGS =
CFLAGS =
FFLAGS =

#
#	Libraries to link and special flags. Ex: LDLIBS = -L/lib -lm, LDFLAGS = -fPIC
LDLIBS = -lblas -llapack -lgfortran -lfltk -lpng -lX11 -ldl
LDFLAGS = -fopenmp -g -p

#
#	Flags used to create a library. Only needed if LIBNAME was defined.
LIBEXT = so
LIBFLAGS = -shared  -Wl,-soname,$(LIBNAME).$(LIBEXT)

#
#	Source code files lists

#
#	Define if the list of source files will be set internally (see bellow)
INTERNAL_SOURCE_DEF = no

#
#	Define if the list of source files will be set by external file (use 'make sourcedef')
EXTERNAL_SOURCE_DEF = yes
EXTERNAL_SOURCE_FILE = source_files.mk

#
#	Define if the list of source files will be automatically searched
AUTO_SOURCE_DEF = no
AUTO_SOURCE_DIR = ../../../..
AUTO_SOURCE_EXT = c cc cpp F


##################################################################################
#
#	Internal source code files definition. Relative or full path required.
#
##################################################################################

ifeq ($(INTERNAL_SOURCE_DEF),yes)
#	C source files (.c)
cSRC =


#	C++ source files (.cc)
ccSRC = 


#	C++ source files (.cpp)
cppSRC = 


#	Fortran source files (.F)
FSRC =
endif






##################################################################################
#
#
#
#	WARNING: change the code below at youe ow: risk!
#
#
#
##################################################################################



##################################################################################
#
#	External source code files definition
#
##################################################################################

ifeq ($(EXTERNAL_SOURCE_DEF),yes)
EXISTS = $(shell if [ -f $(EXTERNAL_SOURCE_FILE) ]; then echo yes; fi )

ifeq ($(EXISTS),yes)
-include $(EXTERNAL_SOURCE_FILE)
else
abort_1:
	@echo Error: No $(EXTERNAL_SOURCE_FILE) found!
	@echo You can use \'make sourcedef\' to create one.
endif
endif



##################################################################################
#
#	If AUTO_SOURCE_DEF then automatically search for source files
#
##################################################################################

ifeq ($(AUTO_SOURCE_DEF),yes)

EXISTS = $(shell if [ -d $(AUTO_SOURCE_DIR) ]; then echo yes; fi )

ifeq ($(EXISTS),yes)
#	C source files (.c)
cSRC += $(shell for src in `find $(AUTO_SOURCE_DIR) -name \*.c`; do cSRC="$$cSRC "$$src; done; echo $$cSRC)

#	C++ source files (.cc)
ccSRC += $(shell for src in `find $(AUTO_SOURCE_DIR) -name \*.cc`; do ccSRC="$$ccSRC "$$src; done; echo $$ccSRC)

#	C++ source files (.cpp)
cppSRC += $(shell for src in `find $(AUTO_SOURCE_DIR) -name \*.cpp`; do cppSRC="$$cppSRC "$$src; done; echo $$cppSRC)

#	Fortran source files (.F)
FSRC += $(shell for src in `find $(AUTO_SOURCE_DIR) -name \*.F`; do FSRC="$$FSRC "$$src; done; echo $$FSRC)
else
abort_2:
	@echo Error: $(AUTO_SOURCE_DIR) is not a directory!
endif
endif



##################################################################################
#
#	General Definitions
#
##################################################################################

#
#	Source and include file paths
#
VPATH := $(shell echo $(dir $(cSRC) $(ccSRC) $(cppSRC) $(FSRC)) | sed "s, ,$'\n,g" | sort | uniq )
INCLUDE = $(patsubst %,-I %,$(subst :, ,$(VPATH))) $(INCFLAGS)

#
#	List of objects and dependency files
#
cSRC_ = $(notdir $(cSRC))
ccSRC_ = $(notdir $(ccSRC))
cppSRC_ = $(notdir $(cppSRC))
FSRC_ = $(notdir $(FSRC))

OBJECTS = $(cSRC_:%.c=$(OBJDIR)/%.o) \
	  $(ccSRC_:%.cc=$(OBJDIR)/%.o) \
          $(cppSRC_:%.cpp=$(OBJDIR)/%.o) \
          $(FSRC_:%.F=$(OBJDIR)/%.o)

DEPENDS = $(cSRC_:%.c=$(OBJDIR)/%.d) \
	  $(ccSRC_:%.cc=$(OBJDIR)/%.d) \
	  $(cppSRC_:%.cpp=$(OBJDIR)/%.d) \
          $(FSRC_:%.F=$(OBJDIR)/%.d)

ifeq ($(strip $(OBJECTS)),)
abort_3:
	@echo Error: No source code defined!
endif

#
#	Definitions for cleaning
#
CLEAN = rm -f
CLEANSRC = $(OBJDIR)/*.[ado] $(PROGNAME) $(LIBNAME)

#
#	Definitions for printing
#
PRINT = lpr -p
PRINTSRC = makefile *.[hcF]*

#
#	Tar definitions
#
TAR = tar
TARDIR = .
TARFLAGS = cf
TARSRC = makefile* *.[hcfF]*

#
#	Archiver definitions
#
AR = ar
ARFLAGS = r
RANLIB = ranlib

#
#	Miscellaneous definitions
#
MAKE = make
RMDIR = rmdir --ignore-fail-on-non-empty

#
#	Location for object, archive, and dependency files
#
RPHOSTTYPE = $(shell uname -s | tr '[:upper:]' '[:lower:]')_$(shell uname -m)
OBJDIR = $(RPHOSTTYPE)

.PRECIOUS: 



##################################################################################
#
#	Main Targets
#
##################################################################################

.PHONY: all

all: $(PROGNAME) $(LIBNAME)

#
#	make a program
#
$(PROGNAME): $(OBJECTS)
	$(LINK) -o $@ $^  $(LDLIBS) $(LDFLAGS)

#
#	make a TEST1 program
#
subphysics: $(OBJECTS)
	$(LINK)  $(CPPFLAGS) $(INCLUDE) ../../../../../main_subphysics.cpp -o $@ $^  $(LDLIBS) $(LDFLAGS)

#
#	Make a library
#
$(LIBNAME): $(OBJECTS)
	$(LINK) $(LIBFLAGS) -o $@ $^  $(LDLIBS) $(LDFLAGS)



##################################################################################
#
#	Compilation Rules
#
##################################################################################

SUFFIXES = .o .d .c .cc .cpp .F
.SUFFIXES:
.SUFFIXES: $(SUFFIXES)

CPPFLAGS := $(CPPFLAGS) $(INCLUDE)
DEPEND = $(CC) -M $(CPPFLAGS)
LINK = $(LD)

#
#	C definitions & rules
#
COMPILE.c = $(CC) $(CFLAGS) $(CPPFLAGS) -c

$(OBJDIR)/%.o: %.c
	$(COMPILE.c) -o $@ $<

$(OBJDIR)/%.d: %.c
	@set -e; rm -f $@; $(DEPEND) $< | sed 's,\($*\)\.o[ :]*,$(OBJDIR)/\1.o $@ : ,g' > $@;


#
#	C++ definitions & rules
#
COMPILE.cc = $(CXX) $(CXXFLAGS) $(CPPFLAGS) -c

$(OBJDIR)/%.o: %.cc
	$(COMPILE.cc) -o $@ $<

$(OBJDIR)/%.o: %.cpp
	$(COMPILE.cc) -o $@ $<

$(OBJDIR)/%.d: %.cc
	@set -e; rm -f $@; $(DEPEND) $< | sed 's,\($*\)\.o[ :]*,$(OBJDIR)/\1.o $@ : ,g' > $@;

$(OBJDIR)/%.d: %.cpp
	@set -e; rm -f $@; $(DEPEND) $< | sed 's,\($*\)\.o[ :]*,$(OBJDIR)/\1.o $@ : ,g' > $@;


#
#	Fortran definitions & rules
#
COMPILE.F = $(FC) $(FFLAGS) $(CPPFLAGS) -c

$(OBJDIR)/%.o: %.F
	$(COMPILE.F) -o $@ $<

$(OBJDIR)/%.d: %.F
	@set -e; rm -f $@; $(DEPEND) $< | sed 's,\($*\)\.o[ :]*,$(OBJDIR)/\1.o $@ : ,g' > $@;




##################################################################################
#
#	General Rules
#
##################################################################################

$(OBJECTS) $(DEPENDS): | $(OBJDIR)
     
$(OBJDIR):
	@echo Creating object repository directory \"$(OBJDIR)\".
	@mkdir $(OBJDIR)


clean: FORCE
	$(CLEAN) $(CLEANSRC)


distclean: FORCE clean
	@$(RMDIR) $(OBJDIR)
	$(CLEAN) $(EXTERNAL_SOURCE_FILE)


print: FORCE
	$(PRINT) $(PRINTSRC)


$(TAGSFILE): $(TAGSSRC)
	$(TAGS) $(TAGSFLAGS) $(TAGSSRC) > $(TAGSFILE)


tar: FORCE
	pwd=`pwd`; name=`basename $$pwd`.tar;\
	$(TAR) $(TARFLAGS) $(TARDIR)/$$name $(TARSRC)


sourcedef: FORCE
	@set -e; > $(EXTERNAL_SOURCE_FILE) ; \
	for ext in $(AUTO_SOURCE_EXT); do \
	    for source in `find $(AUTO_SOURCE_DIR)  -name \*.$$ext | sed -e 's,^\./,,g'`; do  \
	        src="$$src \\\\\n\\t$$source"; \
	    done; \
	    echo "\n"$$ext"SRC += $$src\n" >> $(EXTERNAL_SOURCE_FILE); src=; \
	done


debug:
	$(info Debug Makefile:)
	$(info )
	$(info SOURCEFILES = $(cSRC) $(ccSRC) $(cppSRC) $(FSRC))
	$(info OBJECTS     = $(OBJECTS))
	$(info VPATH       = $(VPATH))
	$(info INCLUDE     = $(INCLUDE))
	$(info DEPENDS     = $(DEPENDS))
	@echo
	
FORCE:


-include $(DEPENDS)



##################################################################################

