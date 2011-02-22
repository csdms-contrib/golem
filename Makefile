
GCC=gcc
CFLAGS=
LIBS=-lm
PROG=golem
VERSION=5.14
EXTRA_DIST= data/example.in \
            data/example.strat \
            data/example2.in \
            data/flat130x130.base
SOURCES= golem.c
OBJS=${SOURCES:.c=.o}

all: ${PROG}

${PROG}: ${OBJS}
	${GCC} ${CFLAGS} ${LIBS} -o ${PROG} ${SOURCES}

dist:
	@mkdir -p ${PROG}-${VERSION}
	@cp ${SOURCES} ${PROG}-${VERSION}
	@for f in ${EXTRA_DIST}; do \
		mkdir -p ${PROG}-${VERSION}/`dirname $$f` && \
		cp -pr $$f ${PROG}-${VERSION}/`dirname $$f`; \
	done
	@cp Makefile ${PROG}-${VERSION}
	@tar cvfz ${PROG}-${VERSION}.tar.gz ${PROG}-${VERSION} 
	@rm -rf ${PROG}-${VERSION}

clean:
	@rm -f ${PROG} ${OBJS} core

.c.o:
	${GCC} ${CFLAGS} -c $<

