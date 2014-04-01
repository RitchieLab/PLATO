Building PLATO
====================

To build PLATO, from this directory, type:

./configure
make
make install

You will need to have root permissions to execute "make install" with the 
default configure options.  Additional helpful options that can be passed
to configure are:

--prefix=<directory>
	This defines the directory to install PLATO into.  PLATO will be copied 
	into $prefix/bin, and by default, the prefix is /usr/local.
	
--enable-debug
	This option enables all debugging symbols and warnings and turns off
	any optimizations.  This is helpful if you encounter a bug and would like
	to fix it, but running in debug mode will likely result in MUCH slower 
	execution times.

NOTE: Due to a bug in Red Hat 6.4 (https://bugzilla.redhat.com/attachment.cgi?id=676407),
PLATO will not configure cleanly with stock BOOST libraries.  To fix this, you
can apply the following patch: https://bugzilla.redhat.com/attachment.cgi?id=676407.
