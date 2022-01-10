#.SILENT:

CPP_CC = g++
C_CC   = gcc

AR = ar
RM = rm

CDEBUG = -g

INC_DIR = ../include
OBJ_DIR = ../obj
LIB_DIR = ../lib

all: dirs $(LIB_DIR)/libCMath.a

dirs:
	@if [ ! -e ../obj ]; then mkdir ../obj; fi
	@if [ ! -e ../lib ]; then mkdir ../lib; fi
	@if [ ! -e ../bin ]; then mkdir ../bin; fi

CPP_SRC = \
CArcToBezier.cpp \
CBezierToLine.cpp \
CDelaunay.cpp \
CGiftWrap.cpp \
CGrahamHull.cpp \
CGrahamHullI.cpp \
CHull3D.cpp \
CInvNorm.cpp \
CMathGen.cpp \
CMathGeom2D.cpp \
CMathGeom3D.cpp \
CMathPoly2D.cpp \
CMathRand.cpp \
CMathRound.cpp \
COffsetBezier.cpp \
CParticle2D.cpp \
CPrime.cpp \
CRigidBody.cpp \
CTriangulate2D.cpp \

CPP_OBJS = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(CPP_SRC))

CPPFLAGS = \
-std=c++14 \
$(CDEBUG) \
-I$(INC_DIR) \
-I../../CMath/include \
-I../../COS/include \
-I../../CUtil/include \
-I. \

clean:
	$(RM) -f $(OBJ_DIR)/*.o
	$(RM) -f $(LIB_DIR)/libCMath.a

.SUFFIXES: .cpp

$(CPP_OBJS): $(OBJ_DIR)/%.o: %.cpp
	$(CPP_CC) -c $< -o $(OBJ_DIR)/$*.o $(CPPFLAGS)

$(LIB_DIR)/libCMath.a: $(CPP_OBJS) $(C_OBJS)
	$(AR) crv $(LIB_DIR)/libCMath.a $(CPP_OBJS) $(C_OBJS)
