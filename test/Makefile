# SPDX-License-Identifier: GPL-3.0-only
#
#  Copyright (C) 2020 Daniel Hu <daddy.of.qq@gmail.com>
#

CXX = g++
CC = gcc

ROOTDIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))
OBJDIR := $(CURDIR)/build

SRCS = main.cpp hash.cpp sha3/sha3.c $(wildcard mbedtls/*.c)

OBJS := $(patsubst %.cpp, $(OBJDIR)/%.o, $(patsubst %.c, $(OBJDIR)/%.o, $(SRCS)))

CFLAGS := $(CFLAGS) -I. -I../ -I -Imbedtls/ -O3 -Werror -DHAVE_LOG -DDEBUG_CAST
CPPFLAGS := $(CFLAGS) -std=c++17

all: $(OBJDIR)/tester

$(OBJDIR)/tester: $(OBJS)
	@echo "LD $@"
	@mkdir -p $(dir $@)
	$(CXX) $(LDFLAGS) -o $@ $^

$(OBJDIR)/%.o: $(ROOTDIR)/%.cpp $(OBJDIR)/%.d
	@echo "CXX $<"
	@mkdir -p $(dir $@)
	@$(CXX) -c $(DEPFLAGS) $(CPPFLAGS) -o $@ -c $< $(INCLUDES)

$(OBJDIR)/%.o: $(ROOTDIR)/%.c $(OBJDIR)/%.d
	@echo "CC $(LTO) $<"
	@mkdir -p $(dir $@)
	@$(CC) -c $(DEPFLAGS) $(CFLAGS) -o $@ -c $< $(INCLUDES)

$(OBJDIR):
	@mkdir -p $@

$(OBJDIR)/%.d: ;
.PRECIOUS: $(OBJDIR)/%.d

-include $(patsubst %.o, %.d, $(OBJS))

clean:
	rm -rf build
