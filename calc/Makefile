
######## Use compiler flags for resutils/Transistor2D #######

include ${HOME}/lib/resutils/Make.common

CXX = ${CXX_TRAN2D}
CXXFLAGS = ${CXXFLAGS_TRAN2D} -I include

######## Object files to be generated #######

SRCS_WOMAIN = ${wildcard src/*.cpp}
OBJS_WOMAIN = ${patsubst %.cpp, %.o, ${SRCS_WOMAIN}}

SRCS_MAIN = ${wildcard bin/*.cpp}
TARGETS = ${patsubst bin/%.cpp, %, ${SRCS_MAIN}}


######## Build rules #######

.PHONY: clean ${TARGETS}

%.o : %.cpp
	${CXX} ${CXXFLAGS} -c $< -o $@

main:  ${OBJS_WOMAIN} bin/main.o
	${CXX} -o bin/$@ bin/$@.o ${OBJS_WOMAIN} ${CXXLDFLAGS}

mainTimeDep:  ${OBJS_WOMAIN} bin/mainTimeDep.o
	${CXX} -o bin/$@ bin/$@.o ${OBJS_WOMAIN} ${CXXLDFLAGS}

mainAlphaDep:  ${OBJS_WOMAIN} bin/mainAlphaDep.o
	${CXX} -o bin/$@ bin/$@.o ${OBJS_WOMAIN} ${CXXLDFLAGS}

mainVgLgDep:  ${OBJS_WOMAIN} bin/mainVgLgDep.o
	${CXX} -o bin/$@ bin/$@.o ${OBJS_WOMAIN} ${CXXLDFLAGS}

mainVgWgDep:  ${OBJS_WOMAIN} bin/mainVgWgDep.o
	${CXX} -o bin/$@ bin/$@.o ${OBJS_WOMAIN} ${CXXLDFLAGS}

mainVgDep-grounded:  ${OBJS_WOMAIN} bin/mainVgDep-grounded.o
	${CXX} -o bin/$@ bin/$@.o ${OBJS_WOMAIN} ${CXXLDFLAGS}

mainVgLgDep-appr:  ${OBJS_WOMAIN} bin/mainVgLgDep-appr.o
	${CXX} -o bin/$@ bin/$@.o ${OBJS_WOMAIN} ${CXXLDFLAGS}

mainGain-VgLgDep:  ${OBJS_WOMAIN} bin/mainGain-VgLgDep.o
	${CXX} -o bin/$@ bin/$@.o ${OBJS_WOMAIN} ${CXXLDFLAGS}

mainGain-VgWgDep:  ${OBJS_WOMAIN} bin/mainGain-VgWgDep.o
	${CXX} -o bin/$@ bin/$@.o ${OBJS_WOMAIN} ${CXXLDFLAGS}

testTR:  ${OBJS_WOMAIN} bin/testTR.o
	${CXX} -o bin/$@ bin/$@.o ${OBJS_WOMAIN} ${CXXLDFLAGS}

clean:
	rm -f ${OBJS_WOMAIN} ${OBJS_MAIN} *~ src/*~ bin/*~ include/*~ bin/*.o
