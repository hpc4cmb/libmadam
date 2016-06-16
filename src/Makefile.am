
# Initialize compile and linking flags

AM_CPPFLAGS =

AM_FCFLAGS =
AM_LDFLAGS =
AM_LIBS = 

# Append CFITSIO, FFTW, BLAS and LAPACK to compile and link flags
AM_FCFLAGS += $(CFITSIO_CPPFLAGS) $(FFTW_CPPFLAGS)
AM_LIBS += $(CFITSIO) $(FFTW) $(BLAS_LIBS) $(LAPACK_LIBS)

# Always need libm
AM_LIBS += -lm

AM_CPPFLAGS += -DMPIWRITECOVMAT

lib_LTLIBRARIES = libmadam.la

libmadam_la_SOURCES = \
	coordinate_conversion.f90 \
	mpi_wrappers.f90 \
        timing.f90 \
        matrix.f90 \
        fourier_fftw_2003.f90 \
	planck_config.F90 \
	planck_types.f90 \
	ls_misc_utils.f90 \
	linebufmod.f90 \
	fitsmod2.f90 \
        healpix_routines.f90 \
        satellite.f90 \
	simulation.f90 \
        commonparam.f90 \
        inputparam.f90 \
        map_routines.f90 \
        submap_transfer.f90 \
	tod_storage.f90 \
        noise_routines.f90 \
        parameter_control.f90 \
        pointing.f90 \
	maps_and_baselines.f90 \
        memory_and_time.f90 \
        maptod_transfer.f90 \
        read_data.f90 \
        output.f90 \
        smadam_routines.f90 \
        covmat_util.f90 \
        covmat.f90 \
        smadam.f90

libmadam_la_LIBADD = $(AM_LIBS)

#noinst_PROGRAMS = test_libmadam test_libcmadam

#test_libmadam_SOURCES = test_libmadam.cpp
#test_libmadam_LDADD = libmadam.la $(AM_LIBS)

#test_libcmadam_SOURCES = test_libcmadam.c
#test_libcmadam_LDADD = libmadam.la $(AM_LIBS)

# force use of the C++ linker...
#nodist_EXTRA_test_libcmadam_SOURCES = dummy.cpp


#TESTS = test_libmadam test_libcmadam test_libmadam.py

clean-local:
	rm -f *.mod
