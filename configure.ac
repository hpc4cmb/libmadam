dnl
dnl This is the Autoconf file for libmadam
dnl
dnl +------------------------
dnl | Initialize package info
dnl +------------------------
AC_INIT([libmadam], [0.2.5], [github:hpc4cmb/libmadam])
AC_CONFIG_SRCDIR([Makefile.am])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_HEADERS(config.h)
AC_CONFIG_MACRO_DIR([m4])

dnl +-------------------------
dnl | If the library source code has changed at all since the last update, 
dnl | then increment revision (c:r:a becomes c:r+1:a).  If any interfaces 
dnl | have been added, removed, or changed since the last update, increment
dnl | current, and set revision to 0.  If any interfaces have been added 
dnl | since the last public release, then increment age.  If any interfaces 
dnl | have been removed since the last public release, then set age to 0.
dnl +-------------------------
LIBMADAM_LT_VERSION="1:3:0"
AC_SUBST(LIBMADAM_LT_VERSION)

dnl +---------------------
dnl | Compiler features
dnl +---------------------
dnl Set default optimization to '-O3' instead of the default '-g -O2'
if test x"${CFLAGS}" = x; then
   CFLAGS="-O3"
fi
if test x"${FCFLAGS}" = x; then
   FCFLAGS="-O3"
fi
AC_CANONICAL_HOST
AC_PROG_INSTALL

dnl +-------------------------
dnl | Fortran support -- use MPI compiler throughout
dnl +-------------------------
dnl AC_PROG_FC
AX_PROG_FC_MPI([],[],[
  AC_MSG_FAILURE([MPI Fortran compiler is required but not found.])
])

dnl Use ".f90" as the Fortran source code suffix
AC_FC_SRCEXT(f90)

dnl +-------------------------
dnl | C support -- use MPI compiler throughout
dnl +-------------------------
dnl AC_PROG_CC
AX_PROG_CC_MPI([],[],[
  AC_MSG_FAILURE([MPI CC compiler is required but not found.])
])

dnl libtool must only be called after the F90/c compilers are overridden
dnl by MPI versions using AX_PROG_CC_MPI and AX_PROG_FC_MPI

LT_INIT

dnl +------------------------------------------------
dnl | This is an F2003 project...
dnl +------------------------------------------------
AC_LANG([Fortran])

dnl +------------------------------------------------
dnl | ... that is OpenMP parallel
dnl +------------------------------------------------
AC_OPENMP
AC_LANG_PUSH([C])
AC_OPENMP
AC_LANG_POP([C])

AC_SUBST([FCFLAGS],["$FCFLAGS $OPENMP_FCFLAGS"])
AC_SUBST([LDFLAGS],["$LDFLAGS $OPENMP_FCFLAGS"])
AC_SUBST([CFLAGS],["$CFLAGS $OPENMP_CFLAGS"])

dnl +------------------------------------------------
dnl | Check for BLAS and LAPACK
dnl +------------------------------------------------
AC_LANG_PUSH([C])
AX_CHECK_BLAS(,AC_MSG_ERROR([Could not find a BLAS library!]))
AX_CHECK_LAPACK(,AC_MSG_ERROR([Could not find a LAPACK library!]))
AC_LANG_POP([C])

dnl +------------------------------------------------
dnl | Check for CFITSIO
dnl +------------------------------------------------
AC_LANG_PUSH([C])
ACX_CFITSIO([], [AC_MSG_ERROR([Could not find the CFITSIO library!])])
AC_LANG_POP([C])

dnl +------------------------------------------------
dnl | Check for FFTW
dnl +------------------------------------------------
AC_LANG_PUSH([C])
ACX_FFTW([],[
  AC_MSG_FAILURE([FFTW is required.])
])
AC_LANG_POP([C])

FFTW_INC=""
if test x"${with_fftw}" != x; then
    FFTW_INC="$with_fftw/include"
fi
AX_F90_HEADER([fftw3.f03],[fftw3.f03],[
    use, intrinsic :: iso_c_binding
    include 'fftw3.f03'
],[${FFTW_INC}],[],[
  AC_MSG_FAILURE([fftw3.f03 is required.])
])
AC_SUBST([FFTW_CPPFLAGS],["$FFTW_CPPFLAGS $F90_HEADER_fftw3_f03"])

dnl +------------------------------------------------
dnl | Set outputs
dnl +------------------------------------------------
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])

dnl +------------------------------------------------
dnl | Generate Makefiles
dnl +------------------------------------------------
AC_OUTPUT

dnl +------------------------------------------------
dnl | Print out detected build options
dnl +------------------------------------------------
AC_MSG_NOTICE([                                           ])
AC_MSG_NOTICE([=========== Build Configuration ===========])
AC_MSG_NOTICE([                                           ])
AC_MSG_NOTICE([  Fortran Compiler          : $FC])
AC_MSG_NOTICE([  Fortran Compile Flags     : $FCFLAGS])
AC_MSG_NOTICE([  Fortran Linking Flags     : $LDFLAGS])
AC_MSG_NOTICE([  C Compiler                : $CC])
AC_MSG_NOTICE([  C Compile Flags           : $CFLAGS])
AC_MSG_NOTICE([  BLAS Link                 : $BLAS_LIBS])
AC_MSG_NOTICE([  LAPACK Link               : $LAPACK_LIBS])
AC_MSG_NOTICE([  FFTW_CPPFLAGS             : $FFTW_CPPFLAGS])
AC_MSG_NOTICE([  FFTW                      : $FFTW])
AC_MSG_NOTICE([  CFITSIO Compile           : $CFITSIO_CPPFLAGS])
AC_MSG_NOTICE([  CFITSIO Link              : $CFITSIO])
AC_MSG_NOTICE([                                           ])
