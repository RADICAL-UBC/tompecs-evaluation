###############################################################################
##  ransampl: library for random number sampling                             ##
##  configure.ac: used by 'autoreconf -i' to prepare for 'configure'         ##
##  (C) Joachim Wuttke 2013-                                                 ##
###############################################################################

AC_INIT([ransampl], [1:1], [http://apps.jcns.fz-juelich.de/ransampl])

# standard configuration for shared libraries
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([foreign]) # don't insert GNU blala files
LT_INIT([disable-static])

AC_PROG_CC_C99
AC_SUBST(AM_CFLAGS,"-pedantic -Wall -Werror")
AC_CONFIG_HEADERS([config.h]) # to avoid endless -D options
AC_HEADER_STDC

# consistency check: is source code present?
AC_CONFIG_SRCDIR([lib/ransampl.c])

# dependency checks
AC_CHECK_HEADERS([gsl/gsl_rng.h],,
	[AC_MSG_ERROR(GSL header gsl_rng.h not found)])

  # -lgsl will only work if we also provide -lgslcblas
AC_SEARCH_LIBS([cblas_dgemm], [gslcblas],,
			 [AC_MSG_ERROR(libgslcblas not found or out of sync)])
AC_SEARCH_LIBS([gsl_rng_uniform], [gsl],,
        [AC_MSG_ERROR(libgsl not found or corrupted: missing gsl_rng_uniform)])

# make these Makefiles
AC_CONFIG_FILES([Makefile lib/Makefile man/Makefile demo/Makefile ransampl.pc])
AC_OUTPUT
