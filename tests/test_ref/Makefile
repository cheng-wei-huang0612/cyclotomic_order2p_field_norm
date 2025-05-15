# Makefile — uses reference source in ../../ref
# -------------------------------------------------

# Adjust these paths if your directory layout changes
REF_DIR := ../../ref
REF_SRC := $(REF_DIR)/norm_ref.c    # reference implementation (compute_norm)

CC      := gcc
CFLAGS  := -O0 -g -std=c11
LDFLAGS := -lflint -lmpfr -lgmp -lm

# Local sources
LOCAL_SRC := poly_sampler.c test_ref.c

SRC     := $(REF_SRC) $(LOCAL_SRC)
OBJ     := $(patsubst %.c,%.o,$(SRC))

all: test_ref

test_ref: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -I$(REF_DIR) -c $< -o $@

run: test_ref
	./test_ref

clean:
	rm -f $(OBJ) test_ref

.PHONY: all run clean
