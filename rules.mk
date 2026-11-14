#Directories
TESTDIR ?= test
SRCDIR ?= src
INCDIR ?= inc
#install dirs
PREFIX ?= /usr
INSTALL_INC ?= $(PREFIX)/include/${PROJECT_NAME}
INSTALL_LIB ?= $(PREFIX)/lib/${PROJECT_NAME}
HEADERS = $(INCDIR)/*/*.h

#Toolchains
CXX ?= g++-13
CXX_SUFFIX ?= cc
CXX_STD ?=-std=c++23
LD = $(CXX)
OPT ?= -O3
MARCH ?= -march=native

#Compile Options
CXXFLAGS = $(CXX_STD) $(OPT) $(MARCH) -MMD -MP -I$(INCDIR) -I /usr/include/matadata
LDLIBS ?=
LDFLAGS ?= 
TEST_LDFLAGS = $(LDFLAGS)
TEST_LDLIBS = $(LDLIBS)

#Shell type
SHELL := /bin/bash

#build rulesS
OBJ = $(SRC:%.$(CXX_SUFFIX)=$(SRCDIR)/%.o)
TEST_SOURCES = $(TEST_TARGETS:%=%.$(CXX_SUFFIX))
TEST_OBJ = $(TEST_SRC:%.$(CXX_SUFFIX)=$(TESTDIR)/%.o)
#objects a test links against
TEST_DEPS ?= $(OBJ)
#header dependencies written by -MMD
DEPS = $(OBJ:.o=.d) $(TEST_OBJ:.o=.d)

define add_section
	$(shell sed -i 's|$(1)[ ]*=.*|$(1)=$(2)|' $(3))
endef

all: $(TARGETS)

tests:$(TEST_TARGETS)

$(TARGETS): $(OBJ)
	$(LD) -o $@ $(OBJ) $(LDFLAGS) $(LDLIBS)

$(TEST_TARGETS): %: $(TESTDIR)/%.o $(TEST_DEPS)
	$(LD) -o $@ $(TESTDIR)/$@.o $(TEST_DEPS) $(TEST_LDFLAGS) $(TEST_LDLIBS)

uninstall:
	$(RM) -r  $(INSTALL_INC) $(INSTALL_LIB)

clean:
	$(RM) $(OBJ) $(TEST_OBJ) $(DEPS)

cleanall: clean
	$(RM) $(TARGETS) $(TEST_TARGETS) Doxyfile *.png *.csv
	$(RM) -r $(OUTPUT_DIRECTORY)

doxygen: Doxyfile
	$(call add_section,PROJECT_NAME,$(PROJECT_NAME),Doxyfile)
	$(call add_section,INPUT,$(INCDIR) $(SRCDIR) $(DOXYGEN_INPUT),Doxyfile)
	$(call add_section,RECURSIVE,YES,Doxyfile)
	$(call add_section,IMAGE_PATH,$(DOXYGEN_IMAGE_PATH),Doxyfile)
	$(call add_section,USE_MDFILE_AS_MAINPAGE,$(DOXYGEN_MAINPAGE),Doxyfile)
	$(call add_section,PROJECT_BRIEF,$(PROJECT_BRIEF),Doxyfile)
	$(call add_section,OUTPUT_DIRECTORY,$(OUTPUT_DIRECTORY),Doxyfile)
	$(call add_section,EXTENSION_MAPPING,h=c++,Doxyfile)
	$(call add_section,BUILTIN_STL_SUPPORT,YES,Doxyfile)
	$(call add_section,EXTRACT_ALL,YES,Doxyfile)
	$(call add_section,EXTRACT_PRIVATE,YES,Doxyfile)
	$(call add_section,INTERNAL_DOCS,YES,Doxyfile)
	$(call add_section,QUIET,YES,Doxyfile)
	$(call add_section,FILE_PATTERNS,*.h *.$(CXX_SUFFIX) *.md,Doxyfile)
	$(call add_section,SOURCE_BROWSER,YES,Doxyfile)
	$(call add_section,GENERATE_LATEX,NO,Doxyfile)
	$(call add_section,HAVE_DOT,YES,Doxyfile)
	$(call add_section,INLINE_INHERITED_MEMBERS,YES,Doxyfile)
	$(call add_section,REPEAT_BRIEF,YES,Doxyfile)
	$(call add_section,WARN_IF_UNDOCUMENTED,NO,Doxyfile)
	$(call add_section,WARN_IF_DOC_ERROR,NO,Doxyfile)
	@doxygen Doxyfile

Doxyfile:
	@doxygen -g $@  > /dev/null

-include $(DEPS)

.PHONY: all clean help

.DEFAULT_GOAL := all
