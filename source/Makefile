# $Id: gfx-config.in 343 2008-09-13 18:34:59Z garland $

CXX = g++
CXXFLAGS = -g -O2 -Wall -Wno-sign-compare -Iinclude -DHAVE_CONFIG_H 
OBJS = Solver.o Particle.o shader.o TinkerToy.o RodConstraint.o SpringForce.o CircularWireConstraint.o imageio.o

project1: $(OBJS)
	$(CXX) -o $@ $^ -lGL -lGLU -lglut -lpng -lglew 
clean:
	rm $(OBJS) project1
