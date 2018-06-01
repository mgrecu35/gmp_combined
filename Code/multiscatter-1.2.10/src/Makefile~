include ../Makefile_include

UTILITY_OBJECTS = variance.o config.o distribution.o anisotropic.o \
		dump_input.o range_spacing_is_regular.o print_algorithms.o
CORE_OBJECTS = singlescatter.o multiscatter.o small_angle.o \
		wide_angle.o tdts.o fast_small_angle.o explicit.o \
		wide_angle_regrid.o
ADJOINT_OBJECTS = singlescatter_AD.o multiscatter_AD.o \
		wide_angle_AD.o tdts_AD.o fast_small_angle_AD.o \
		anisotropic_AD.o variance_AD.o
EXTRA_OBJECTS = fast_small_angle_lag.o jacobian.o jacobian_linear.o \
		fortran_interface.o

OBJECTS = $(UTILITY_OBJECTS) $(CORE_OBJECTS) $(ADJOINT_OBJECTS) $(EXTRA_OBJECTS)

LIBMS = ../lib/libmultiscatter.a

PROGRAM = ../bin/multiscatter
PROGRAMF = ../bin/multiscatterf
PROGRAM_FDEMO = ../bin/multiscatter_fortran_demo

all: $(LIBMS) $(PROGRAM) $(PROGRAMF) $(PROGRAM_FDEMO)

$(LIBMS): $(OBJECTS)
	ar -r $(LIBMS) $(OBJECTS)


$(PROGRAM): $(LIBMS) multiscatter_ascii.o
	$(CC) $(CFLAGS) $(INCLUDES) multiscatter_ascii.o $(LIBS) -o $(PROGRAM)


$(PROGRAMF): $(LIBMS) multiscatterf_ascii.f90 ../include/multiscatter.inc
	$(FC) $(FFLAGS) $(INCLUDES) multiscatterf_ascii.f90 $(LIBS) \
	   -o $(PROGRAMF)

$(PROGRAM_FDEMO): $(LIBMS) multiscatter_fortran_demo.f90 ../include/multiscatter.inc
	$(FC) $(FFLAGS) $(INCLUDES) multiscatter_fortran_demo.f90 $(LIBS) \
	   -o $(PROGRAM_FDEMO)

%.o: %.c ms.h ../include/multiscatter.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

clean:
	rm -f $(OBJECTS) $(LIBMS) $(PROGRAM) $(PROGRAMF) \
	   ${PROGRAM_FDEMO} multiscatter_ascii.o

clean_autosaves:
	rm -f *~

*.f90 *.c: ms.h ../include/multiscatter.h
*.f90: ../include/multiscatter.inc
