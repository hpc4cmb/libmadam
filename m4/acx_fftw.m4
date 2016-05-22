#
# SYNOPSIS
#
#   ACX_FFTW([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro looks for a version of FFTW3.  The FFTW_CPPFLAGS and FFTW
#   output variables hold the compile and link flags.
#
#   To link an application with FFTW, you should link with:
#
#   	$FFTW
#
#   The user may use:
# 
#       --with-fftw-cpp=<flags> --with-fftw-libs=<flags>
#
#   to manually specify the include and linking flags.
#
#   ACTION-IF-FOUND is a list of shell commands to run if the FFTW library is
#   found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it is
#   not found. If ACTION-IF-FOUND is not specified, the default action will
#   define HAVE_FFTW and set output variables above.
#
#   This macro requires autoconf 2.50 or later.
#
# LAST MODIFICATION
#
#   2011-03-31
#
# COPYLEFT
#
#   Copyright (c) 2009-2011 Theodore Kisner <tskisner@lbl.gov>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Macro Archive. When you make and
#   distribute a modified version of the Autoconf Macro, you may extend this
#   special exception to the GPL to apply to your modified version as well.
#
#  Revisions:
#    2016-05-21 Adapted to Fortran

AC_DEFUN([ACX_FFTW], [
AC_PREREQ(2.50)

acx_fftw_ok=no
acx_fftw_threads=no

FFTW_CPPFLAGS=""
FFTW=""

AC_ARG_WITH(fftw-cpp, [AC_HELP_STRING([--with-fftw-cpp=<flags>], [use FFTW preprocessing flags <flags>.  Set to "no" to disable.])])

AC_ARG_WITH(fftw-libs, [AC_HELP_STRING([--with-fftw-libs=<flags>], [use FFTW linking flags <flags>.  Set to "no" to disable.])])

if test x"$with_fftw_cpp" != x; then
   if test x"$with_fftw_cpp" != xno; then
      FFTW_CPPFLAGS="$with_fftw_cpp"
   else
      acx_fftw_ok=disable
   fi
fi

if test x"$with_fftw_libs" != x; then
   if test x"$with_fftw_libs" != xno; then
      FFTW="$with_fftw_libs"
   else
      acx_fftw_ok=disable
   fi
fi

if test $acx_fftw_ok = disable; then
   echo "**** FFTW explicitly disabled by configure."
else

   # Save environment

   acx_fftw_save_CPPFLAGS="$CPPFLAGS"
   acx_fftw_save_LIBS="$LIBS"

   # Test serial compile and linking

   # First, try user-specified location or linked-by-default (for example using
   # a compiler wrapper script)

   CPPFLAGS="$CPPFLAGS $FFTW_CPPFLAGS"
   LIBS="$FFTW $acx_fftw_save_LIBS -lm"

   AC_CHECK_HEADERS([fftw3.h])

   AC_MSG_CHECKING([for fftw_malloc in user specified location])
   AC_TRY_LINK_FUNC(fftw_malloc, [acx_fftw_ok=yes;AC_DEFINE(HAVE_FFTW,1,[Define if you have the FFTW library.])], [])
   AC_MSG_RESULT($acx_fftw_ok)

   if test $acx_fftw_ok = yes; then
      acx_fftw_threads_ok=no
      AC_MSG_CHECKING([for fftw_plan_with_nthreads in user specified location])
      AC_TRY_LINK_FUNC(fftw_plan_with_nthreads, [acx_fftw_threads_ok=yes;AC_DEFINE(HAVE_FFTW_THREADS,1,[Define if you have the FFTW threads library.])], [])
      AC_MSG_RESULT($acx_fftw_threads_ok)
   fi
   
   if test $acx_fftw_ok = no; then
      # try pkg-config
      if test `pkg-config --exists fftw`; then
         FFTW_CPPFLAGS=`pkg-config --cflags fftw`
         FFTW=`pkg-config --libs fftw`
         CPPFLAGS="$acx_fftw_save_CPPFLAGS $FFTW_CPPFLAGS"
         LIBS="$FFTW $acx_fftw_save_LIBS -lm"
	
         AC_CHECK_HEADERS([fftw3.h])

         AC_MSG_CHECKING([for fftw_malloc in pkg-config location])
         AC_TRY_LINK_FUNC(fftw_malloc, [acx_fftw_ok=yes;AC_DEFINE(HAVE_FFTW,1,[Define if you have the FFTW library.])], [])
         AC_MSG_RESULT($acx_fftw_ok)

         if test $acx_fftw_ok = yes; then
            acx_fftw_threads_ok=no
            LIBS="-lfftw3_threads $FFTW $acx_fftw_save_LIBS -lm"
            AC_MSG_CHECKING([for fftw_plan_with_nthreads in pkg-config location])
            AC_TRY_LINK_FUNC(fftw_plan_with_nthreads, [acx_fftw_threads_ok=yes;AC_DEFINE(HAVE_FFTW_THREADS,1,[Define if you have the FFTW threads library.])], [])
            AC_MSG_RESULT($acx_fftw_threads_ok)

            if test $acx_fftw_threads_ok = yes; then
               FFTW="-lfftw3_threads $FFTW"
            fi
         fi

      fi

      if test $acx_fftw_ok = no; then
         # try just default linking flags
         FFTW_CPPFLAGS=
         FFTW="-lfftw3"

         CPPFLAGS="$acx_fftw_save_CPPFLAGS $FFTW_CPPFLAGS"
         LIBS="$FFTW $acx_fftw_save_LIBS -lm"
	
         AC_CHECK_HEADERS([fftw3.h])

         AC_MSG_CHECKING([for fftw_malloc in pkg-config location])
         AC_TRY_LINK_FUNC(fftw_malloc, [acx_fftw_ok=yes;AC_DEFINE(HAVE_FFTW,1,[Define if you have the FFTW library.])], [])
         AC_MSG_RESULT($acx_fftw_ok)

         if test $acx_fftw_ok = yes; then
            acx_fftw_threads_ok=no
            LIBS="-lfftw3_threads $FFTW $acx_fftw_save_LIBS -lm"
            AC_MSG_CHECKING([for fftw_plan_with_nthreads in pkg-config location])
            AC_TRY_LINK_FUNC(fftw_plan_with_nthreads, [acx_fftw_threads_ok=yes;AC_DEFINE(HAVE_FFTW_THREADS,1,[Define if you have the FFTW threads library.])], [])
            AC_MSG_RESULT($acx_fftw_threads_ok)

            if test $acx_fftw_threads_ok = yes; then
               FFTW="-lfftw3_threads $FFTW"
            fi
         fi
      fi
   fi

   if test $acx_fftw_ok = no; then
      FFTW_CPPFLAGS=""
      FFTW=""
   fi

   # Restore environment

   LIBS="$acx_fftw_save_LIBS"
   CPPFLAGS="$acx_fftw_save_CPPFLAGS"

fi

# Define exported variables

AC_SUBST(FFTW_CPPFLAGS)
AC_SUBST(FFTW)

# Execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:

if test x"$acx_fftw_ok" = xyes; then
   ifelse([$1],,[echo "**** Enabling FFTW support."],[$1])
else
   ifelse([$2],,[echo "**** FFTW not found - disabling support."],[$2])
fi

])
