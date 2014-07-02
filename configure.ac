# Ritchie Lab configure.ac template

AC_PREREQ(2.59)
AC_INIT(PLATO, 2.0.0, software@ritchielab.psu.edu)

# Checks for --enable-debug added
AX_CHECK_ENABLE_DEBUG
#AC_PROG_CC
#AC_PROG_CPP
#AC_PROG_CXX

AX_PROG_CC_MPI(,[AC_DEFINE([HAVE_CC_MPI],[],[C MPI Found])],[AC_DEFINE([NO_CC_MPI],[],[No C MPI Found])])
AX_PROG_CXX_MPI(,[AC_DEFINE([HAVE_CXX_MPI],[],[CXX MPI Found])],[AC_DEFINE([NO_CXX_MPI],[],[No CXX MPI Found])])

RL_CONFIGURE([1.42],[HAVE_BOOST HAVE_BOOST_THREAD HAVE_BOOST_PROGRAM_OPTIONS \
                     HAVE_BOOST_FILESYSTEM HAVE_BOOST_REGEX HAVE_BOOST_IOSTREAMS \ 
                     HAVE_GSL])

AM_INIT_AUTOMAKE([1.8 foreign])

AC_CONFIG_SRCDIR([main.cpp])
AC_CONFIG_HEADER([config.h])


AC_CONFIG_FILES([Makefile
                 method_lib/Makefile
                ])
AC_OUTPUT
