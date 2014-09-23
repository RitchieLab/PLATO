# Ritchie Lab configure.ac template

AC_PREREQ(2.59)
AC_INIT(PLATO, 2.1.0, software@ritchielab.psu.edu)

# Checks for --enable-debug added
AX_CHECK_ENABLE_DEBUG
#AC_PROG_CC
#AC_PROG_CPP
#AC_PROG_CXX


AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
[compile with MPI (parallelization) support. If none is found, MPI is not used. Default: auto]) ],
,[with_mpi=auto])

AX_PROG_CXX_MPI([test x"$with_mpi" != xno],[AC_DEFINE([HAVE_CXX_MPI],[],[CXX MPI Found])],[
   AC_DEFINE([NO_CXX_MPI],[],[No C MPI Found])
   if test x"$with_mpi" = xyes; then
     AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
   else
     AC_MSG_WARN([No MPI compiler found, won't use MPI.])
   fi
])

AX_PROG_CC_MPI([test x"$with_mpi" != xno],[AC_DEFINE([HAVE_CC_MPI],[],[C MPI Found])],[
   AC_DEFINE([NO_CC_MPI],[],[No C MPI Found])
   if test x"$with_mpi" = xyes; then
     AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
   else
     AC_MSG_WARN([No MPI compiler found, won't use MPI.])
   fi
])


RL_CONFIGURE([1.42],[HAVE_BOOST HAVE_BOOST_THREAD HAVE_BOOST_PROGRAM_OPTIONS \
                     HAVE_BOOST_FILESYSTEM HAVE_BOOST_REGEX HAVE_BOOST_IOSTREAMS \ 
                     HAVE_BOOST_SERIALIZATION HAVE_GSL])

AM_INIT_AUTOMAKE([1.8 foreign])

AC_CONFIG_SRCDIR([main.cpp])
AC_CONFIG_HEADER([config.h])


AC_CONFIG_FILES([Makefile
                 method_lib/Makefile
                ])
AC_OUTPUT
