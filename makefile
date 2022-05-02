target := swm.x

iflags :=
iflags += -I${GEOLYTICAL}/include
iflags += -I${HYWALL}/include
iflags += -I${PTL}/include

lflags :=
lflags += -L${GEOLYTICAL}/lib -lgeolytical
lflags += -L${HYWALL}/lib -lHyWall
lflags += -L${PTL}/lib -lPTL

main:
	mpicxx -O3 ${iflags} main.cc -o ${target} ${lflags}

run:
	./${target}