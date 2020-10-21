CXX = g++
CC = gcc

ROOTDIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))
OBJDIR := $(CURDIR)/build

SRCS = demo.cpp test/hash.cpp test/sha3/sha3.c $(wildcard test/mbedtls/*.c)

OBJS := $(patsubst %.cpp, $(OBJDIR)/%.o, $(patsubst %.c, $(OBJDIR)/%.o, $(SRCS)))

CFLAGS := $(CFLAGS) -I. -I./test -Itest/mbedtls/ -O3 -Werror 
CPPFLAGS := $(CFLAGS) -std=c++17

all: $(OBJDIR)/demo

$(OBJDIR)/demo: $(OBJS)
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
