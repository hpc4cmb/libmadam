
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

include_HEADERS = madam.h

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
        covmat_util.F90 \
        covmat.F90 \
        smadam.F90

# Automake does not recognize Fortran dependencies so we list them
# here explicitly

planck_config.lo : planck_config.F90

planck_types.lo : planck_types.f90 planck_config.lo

ls_misc_utils.lo : ls_misc_utils.f90 planck_config.lo

fitsmod2.lo : fitsmod2.f90 planck_config.lo planck_types.lo ls_misc_utils.lo \
	linebufmod.lo

commonparam.lo : commonparam.f90 planck_config.lo mpi_wrappers.lo

inputparam.lo : inputparam.f90 commonparam.lo mpi_wrappers.lo

map_routines.lo : map_routines.f90 commonparam.lo matrix.lo mpi_wrappers.lo \
	timing.lo

submap_transfer.lo : submap_transfer.f90 commonparam.lo mpi_wrappers.lo

parameter_control.lo : parameter_control.f90 commonparam.lo \
	noise_routines.lo pointing.lo mpi_wrappers.lo maps_and_baselines.lo \
	memory_and_time.lo

pointing.lo : pointing.f90 commonparam.lo mpi_wrappers.lo memory_and_time.lo

matrix.lo : matrix.f90 commonparam.lo

maps_and_baselines.lo : maps_and_baselines.f90 commonparam.lo mpi_wrappers.lo \
	memory_and_time.lo

tod_storage.lo : tod_storage.f90 commonparam.lo mpi_wrappers.lo memory_and_time.lo

noise_routines.lo : noise_routines.f90 commonparam.lo fourier_fftw_2003.lo \
	mpi_wrappers.lo timing.lo tod_storage.lo memory_and_time.lo

memory_and_time.lo : memory_and_time.f90 commonparam.lo mpi_wrappers.lo

output.lo : output.f90 commonparam.lo submap_transfer.lo \
	mpi_wrappers.lo timing.lo maps_and_baselines.lo tod_storage.lo \
	fitsmod2.lo

maptod_transfer.lo : maptod_transfer.f90 commonparam.lo mpi_wrappers.lo \
	memory_and_time.lo

smadam_routines.lo : smadam_routines.f90 commonparam.lo mpi_wrappers.lo \
	maptod_transfer.lo noise_routines.lo map_routines.lo \
	pointing.lo timing.lo

read_data.lo : read_data.f90 commonparam.lo submap_transfer.lo \
	fitsmod2.lo pointing.lo maps_and_baselines.lo coordinate_conversion.lo \
	healpix_routines.lo tod_storage.lo timing.lo \
	mpi_wrappers.lo maptod_transfer.lo

covmat_util.lo : covmat_util.F90 planck_config.lo mpi_wrappers.lo

covmat.lo : covmat.F90 commonparam.lo fourier_fftw_2003.lo pointing.lo covmat_util.lo \
	noise_routines.lo mpi_wrappers.lo

smadam.lo : smadam.F90 commonparam.lo inputparam.lo \
	parameter_control.lo pointing.lo maptod_transfer.lo submap_transfer.lo \
	fourier_fftw_2003.lo noise_routines.lo smadam_routines.lo \
	map_routines.lo read_data.lo output.lo mpi_wrappers.lo \
	maps_and_baselines.lo tod_storage.lo memory_and_time.lo timing.lo \
	covmat.lo covmat_util.lo


libmadam_la_LIBADD = $(AM_LIBS)

#libmadam_la_LDFLAGS = -Wl,--no-undefined

noinst_PROGRAMS = test_libmadam

test_libmadam_SOURCES = test_libmadam.c
test_libmadam_LDADD = libmadam.la $(AM_LIBS)

AM_TESTS_ENVIRONMENT = \
	mkdir -p maps;

TESTS = test_libmadam

clean-local:
	rm -rf *.mod maps/*
