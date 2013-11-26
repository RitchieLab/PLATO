# Ritchie Lab configure.ac template

AC_PREREQ(2.59)
AC_INIT(PLATO, 2.0.0, software@ritchielab.psu.edu)

# Checks for --enable-debug added
AX_CHECK_ENABLE_DEBUG

AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX

RL_CONFIGURE([1.36],[HAVE_BOOST HAVE_BOOST_THREAD HAVE_GSL])

AM_INIT_AUTOMAKE([1.8 foreign])

AC_CONFIG_SRCDIR([main.cpp])
AC_CONFIG_HEADER([config.h])


AC_CONFIG_FILES([Makefile
                 method_lib/Makefile
                ])
AC_OUTPUT
