ESCAPE_HOME := ../

TARGETS :=  test_cc test_chic

OBJECTS := $(TARGETS:%=%.o)

include ../common.mk

INCLUDES += -I..
LDFLAGS += -L..

test_% : test_%.o ../libescape.a
	$(COMPILE_AND_LINK)
