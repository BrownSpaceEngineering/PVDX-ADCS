# Simple Makefile for PVDX-ADCS
# Builds all .c files in src/ and links into out/adcs executable

CC = gcc
CFLAGS = -Iinclude -Wall -Wextra -O2 -std=c11
SRCS = $(wildcard src/*.c)
OBJS := $(patsubst src/%.c,out/output_files/%.o,$(SRCS))
TARGET = out/adcs

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@

out/output_files/%.o: src/%.c
	@mkdir -p out/
	@mkdir -p out/output_files
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
