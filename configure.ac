# Ritchie Lab configure.ac template

AC_PREREQ(2.59)
AC_INIT(PLATO, 1.2.0, software@ritchielab.psu.edu)

AM_INIT_AUTOMAKE([1.8 foreign])

AC_CONFIG_SRCDIR([wasp.cc])
AC_CONFIG_HEADER([config.h])

RL_CONFIGURE([1.40],[HAVE_BOOST HAVE_BOOST_THREAD])

AC_CONFIG_FILES([Makefile
                 method_lib/Makefile
                ])
AC_OUTPUT
