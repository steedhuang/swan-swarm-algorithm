CXXFLAGS = `pkg-config --cflags playerc++` -g -std=gnu++11
CFLAGS = `pkg-config --cflags playerc` -g
LDLIBS = `pkg-config --libs playerc++`
CC = g++

all: main7

