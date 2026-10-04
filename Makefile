# ir-core-lab -- build, test and simulation
#
# Deliberately free of subdirectory rules: this Makefile must also work with
# GNU Make 3.80 as shipped with some Windows toolchains, where $(shell mkdir)
# tricks are unavailable.  Objects are not written to an intermediate
# directory; each target compiles straight to an executable in the root.
#
# CC is set with `=` (not `?=`) so `make CC=clang` still overrides it, while
# the default stays a plain `gcc`.  On Windows/MSYS2 use, for example:
#     make CC=/d/msys2/ucrt64/bin/gcc.exe
#
# -ffp-contract=off is mandatory: without it clang (and newer gcc) may fuse
# a*b+c into one FMA, which changes the last bits of the NUC and TEC results
# and makes the floating point assertions fail on macOS only.

CC      = gcc
CFLAGS  = -std=c99 -O2 -Wall -Wextra -Werror -ffp-contract=off
LDLIBS  = -lm

# Windows needs the .exe suffix in both the -o argument and the target name,
# otherwise make rebuilds every time (gcc creates test_ir.exe while the rule
# claims to produce test_ir).
ifeq ($(OS),Windows_NT)
  EXE       = .exe
  RUNPREFIX =
else
  EXE       =
  RUNPREFIX = ./
endif

BIN_TEST = test_ir$(EXE)
BIN_SIM  = sim_ir$(EXE)

SRC = src/ir_lab.c \
      src/tec_ctrl.c \
      src/nuc.c \
      src/badpix.c \
      src/agc.c \
      src/frame.c \
      src/vcs_pwr.c

HDR = src/ir_lab.h \
      src/tec_ctrl.h \
      src/nuc.h \
      src/badpix.h \
      src/agc.h \
      src/frame.h \
      src/vcs_pwr.h

all: $(BIN_TEST) $(BIN_SIM)

$(BIN_TEST): test/test_ir.c $(SRC) $(HDR)
	$(CC) $(CFLAGS) -Isrc -o $(BIN_TEST) test/test_ir.c $(SRC) $(LDLIBS)

$(BIN_SIM): src/main_sim.c $(SRC) $(HDR)
	$(CC) $(CFLAGS) -Isrc -o $(BIN_SIM) src/main_sim.c $(SRC) $(LDLIBS)

# Run the self tests.  results/ is part of the repository, so nothing here
# needs to create a directory.
test: $(BIN_TEST)
	@$(RUNPREFIX)$(BIN_TEST)

# Run the end to end simulation; it writes results/*.csv and metrics.txt.
# Do NOT pipe this into `head` in CI: the broken pipe raises SIGPIPE and the
# job reports a non-zero exit even though the program itself succeeded.
sim: $(BIN_SIM)
	@$(RUNPREFIX)$(BIN_SIM)

clean:
	rm -f test_ir test_ir.exe sim_ir sim_ir.exe
	rm -f _test_write.tmp

.PHONY: all test sim clean
