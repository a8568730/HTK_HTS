# ----------------------------------------------------------- 
#                                                             
#                          ___                                
#                       |_| | |_/   SPEECH                    
#                       | | | | \   RECOGNITION               
#                       =========   SOFTWARE                  
#                                                             
#                                                             
# ----------------------------------------------------------- 
#         Copyright: Cambridge University
#          1995-2005 Engineering Department
#                    http://htk.eng.cam.ac.uk
#                    http://mi.eng.cam.ac.uk
#                 
#   Use of this software is governed by a License Agreement   
#    ** See the file License for the Conditions of Use  **    
#    **     This banner notice must not be removed      **    
#                                                             
# ----------------------------------------------------------- 
# File: configure.ac source file for generating configure script
# ----------------------------------------------------------- 

#  *** THIS IS A MODIFIED VERSION OF HTK ***                        #
#  ---------------------------------------------------------------  #
#           The HMM-Based Speech Synthesis System (HTS)             #
#                       HTS Working Group                           #
#                                                                   #
#                  Department of Computer Science                   #
#                  Nagoya Institute of Technology                   #
#                               and                                 #
#   Interdisciplinary Graduate School of Science and Engineering    #
#                  Tokyo Institute of Technology                    #
#                                                                   #
#                     Copyright (c) 2001-2007                       #
#                       All Rights Reserved.                        #
#                                                                   #
#  Permission is hereby granted, free of charge, to use and         #
#  distribute this software in the form of patch code to HTK and    #
#  its documentation without restriction, including without         #
#  limitation the rights to use, copy, modify, merge, publish,      #
#  distribute, sublicense, and/or sell copies of this work, and to  #
#  permit persons to whom this work is furnished to do so, subject  #
#  to the following conditions:                                     #
#                                                                   #
#    1. Once you apply the HTS patch to HTK, you must obey the      #
#       license of HTK.                                             #
#                                                                   #
#    2. The source code must retain the above copyright notice,     #
#       this list of conditions and the following disclaimer.       #
#                                                                   #
#    3. Any modifications to the source code must be clearly        #
#       marked as such.                                             #
#                                                                   #
#  NAGOYA INSTITUTE OF TECHNOLOGY, TOKYO INSTITUTE OF TECHNOLOGY,   #
#  HTS WORKING GROUP, AND THE CONTRIBUTORS TO THIS WORK DISCLAIM    #
#  ALL WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL       #
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT   #
#  SHALL NAGOYA INSTITUTE OF TECHNOLOGY, TOKYO INSTITUTE OF         #
#  TECHNOLOGY, HTS WORKING GROUP, NOR THE CONTRIBUTORS BE LIABLE    #
#  FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY        #
#  DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,  #
#  WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTUOUS   #
#  ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR          #
#  PERFORMANCE OF THIS SOFTWARE.                                    #
#                                                                   #
#  ---------------------------------------------------------------  #

# After editing this file, run "autoconf" to regenerate "configure".
AC_INIT(HTS,2.0.1, hts-users@sp.nitech.ac.jp)
AC_COPYRIGHT(Copyright 2003-2005 University of Cambridge)
AC_COPYRIGHT(Copyright 2001-2007 Nagoya Institute of Technology and Tokyo Institute of Technology)

AC_PREFIX_DEFAULT(/usr/local/HTS-2.0.1)

# Checks for programs.
AC_PROG_MAKE_SET
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_RANLIB

# Check for ar
AN_MAKEVAR([AR], [AC_PROG_AR])
AN_PROGRAM([ar], [AC_PROG_AR])
AC_DEFUN([AC_PROG_AR], [AC_CHECK_TOOL(AR, ar, :)])
AC_PROG_AR

# Checks for libraries.
AC_CHECK_LIB([X11], [main])
AC_CHECK_LIB([m], [main])

# Checks for header files.
AC_CONFIG_HEADERS
AC_PATH_X
AC_HEADER_STDC
AC_CHECK_HEADERS(errno.h fcntl.h float.h limits.h malloc.h memory.h stdlib.h string.h sys/ioctl.h sys/socket.h sys/time.h unistd.h)

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_VOLATILE
AC_TYPE_SIZE_T
AC_HEADER_TIME
AC_STRUCT_TM

# Checks for library functions.
AC_PROG_GCC_TRADITIONAL
AC_FUNC_MEMCMP
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_TYPE_SIGNAL
AC_FUNC_VPRINTF
AC_HEADER_STDBOOL
AC_CHECK_FUNCS(floor gettimeofday memmove memset modf pow socket sqrt strchr strcspn strrchr strspn strstr strtol strdup)


AC_CANONICAL_HOST
AC_C_BIGENDIAN
AC_CHECK_SIZEOF(int *)


# Traditional HTK build?
# if enabled this will install in $bindir.$host_cpu
AC_ARG_ENABLE(trad-htk,
		AS_HELP_STRING([--enable-trad-htk],
		[HTK installed into $prefix/bin.$cpu and installed when compiled]))

# Build HDecode
AC_ARG_ENABLE(hdecode,
		AS_HELP_STRING([--enable-hdecode],
		[build large vocabulary recognition system (HDecode)]))

# Build LM Tools
AC_ARG_ENABLE(hlmtools,
	AS_HELP_STRING([--enable-hlmtools],
	[build Language Modelling tools]))

# Disable building for X
# i.e. use HGraf.null, and don't build HSLab
AC_ARG_ENABLE(hslab,
		AS_HELP_STRING([--disable-hslab],
		[don't build HSLab]))

# Build LV Recognition
AC_ARG_ENABLE(htkbook,
		AS_HELP_STRING([--enable-htkbook],
		[build HTK book]))

# Disable building hts_engine
AC_ARG_ENABLE(hts_engine,
	AS_HELP_STRING([--disable-hts_engine],
	[don't build hts_engine]))

# Use -Wall if using gcc
case "$CC" in
	 gcc*)
		CFLAGS="-Wall -Wno-switch $CFLAGS"
		;;		
esac

# customise for htk comilpation
trad_bin_dir=$host
case "$host_cpu" in
	 powerpc)
	 	host=darwin
		trad_bin_dir=$host
		;;
esac

case "$host" in
		*x86_64*linux*)
		CFLAGS="-ansi -D_SVID_SOURCE -DOSS_AUDIO -D'ARCH=\"$host_cpu\"' $CFLAGS"
		AC_MSG_CHECKING(whether to compile in 32-bit or 64-bit mode)
		if test x$ac_cv_sizeof_int_p = x"4"; then
			AC_MSG_RESULT(32-bit)
				LDFLAGS="-L/usr/X11R6/lib $LDFLAGS"
		elif test x$ac_cv_sizeof_int_p = x"8"; then
			AC_MSG_RESULT(64-bit)
			LDFLAGS="-L/usr/X11R6/lib64 $LDFLAGS"
		fi
				ARCH=linux
				trad_bin_dir=linux
		;;
		*linux*)
				CFLAGS="-ansi -D_SVID_SOURCE -DOSS_AUDIO -D'ARCH=\"$host_cpu\"' $CFLAGS"
				LDFLAGS="-L/usr/X11R6/lib $LDFLAGS"
				ARCH=linux
				trad_bin_dir=linux
		;;
		*pc-solaris*)
				CFLAGS="-ansi -I/usr/openwin/include -O2 -DSOLARIS16_AUDIO -DBSD_COMP -D'ARCH=\"$host_cpu\"' $CFLAGS"
				LDFLAGS="-L/usr/openwin/lib -ldl -lsocket -lnsl -lintl $LDFLAGS"
		;;
		*sun*)
				CFLAGS="-ansi -I/usr/openwin/include -O2 -DSOLARIS16_AUDIO -DBSD_COMP -D'ARCH=\"$host_cpu\"' $CFLAGS"
				LDFLAGS="-L/usr/openwin/lib -ldl -lsocket -lnsl -lintl $LDFLAGS"
		;;
		*sgi*)
				CFLAGS="-O -DSGI_AUDIO  -D'ARCH=\"$host_cpu\"' $CFLAGS"
				LDFLAGS="-laudio -lmalloc $LDFLAGS"
		;;
		*hp*)
				CFLAGS="+e -O -I/usr/include/X11R5 -I/usr/include/audio -I/opt/audio/include -D_HPUX_SOURCE -Aa -DHPUX_AUDIO -DUSE_TIMER  -D'ARCH=\"$host_cpu\"' $CFLAGS"
				LDFLAGS="-L/opt/audio/lib -lAlib $LDFLAGS"
		;;
		*cygwin)
				CFLAGS="-ansi -DCYGWIN -D'ARCH=\"ASCII\"' $CFLAGS"
				LDFLAGS="-L/usr/X11R6/lib $LDFLAGS"
				ARCH=ASCII
		CPU=cygwin
				BINARY_EXTENSION=.exe
		;;				
	*darwin*)
		CFLAGS="-ansi -g -O2 -DNO_AUDIO -D'ARCH=\"darwin\"' -I/usr/include/malloc $CFLAGS"
				LDFLAGS="-L/usr/X11R6/lib $LDFLAGS"
				ARCH=darwin
				Objcopy=echo
				PRILF="-x"
				CPU=darwin
				SHRLF="-shared"
				LIBEXT=dylib
		;;				
esac

# Process customised options
make_all=htktools
make_install=install-htktools
build_notes="Build notes: "
# Tweak bindir and libdir for backwards compatibility 
# with previous versions of HTK
case "$enable_trad_htk" in
     yes)
	bindir=${bindir}.${trad_bin_dir}
	libdir=${libdir}.${trad_bin_dir}
	;;
esac
case "$enable_hlmtools" in
	yes)
	make_all="$make_all hlmtools"
	make_install="$make_install install-hlmtools"
	build_notes="$build_notes Language Modelling tools will be built."
	;;
	*)
		build_notes="$build_notes The Language Modelling tools will not be built. You can build them manually later by running 'make hlmtools install-hlmtools'."
	;;
esac
case "$enable_hdecode" in
     yes)
	make_all="$make_all hdecode"
	make_install="$make_install install-hdecode"
	build_notes="$build_notes Large Vocabulary recogniser HDecode will be built."
	;;
     *)
		build_notes="$build_notes HDecode will not be built. You can build it manually later by running 'make hdecode install-hdecode'."
	;;
esac
case "$enable_htkbook" in
	 yes)
	 make_all="$make_all book"
	 make_install="$make_install install-book"
	 build_notes="$build_notes The HTK Book will be built."
	 ;;
esac
case "$enable_hslab" in
     no)
	HSLAB=
	HGRAF=HGraf.null.o
	build_notes="$build_notes HSLab will not be built.  If you want to change this then run 'configure' again with the '--enable-hslab' option."
	;;
     *)
	HSLAB=HSLab
	HGRAF=HGraf.o
	build_notest="$build_notes GUI tool HSLab will be built."
	;;
esac
case "$enable_hts_engine" in
	no)
		build_notes="$build_notes hts_engine will not be built. You can build it manually later by running 'make hts_engine install-hts_engine'."
	;;
	*)
		make_all="$make_all hts_engine"
		make_install="$make_install install-hts_engine"
		build_notes="$build_notes hts_engine will be built."
	;;
esac

AC_SUBST(TRADHTK, $enable_trad_htk)
AC_SUBST(TRADHTKBIN, $trad_bin_dir)
AC_SUBST(make_all, $make_all)
AC_SUBST(make_install, $make_install)
AC_SUBST(HSLAB, $HSLAB)
AC_SUBST(HGRAF, $HGRAF)
AC_SUBST(BINARY_EXTENSION, $BINARY_EXTENSION)

AC_OUTPUT(HTKLib/Makefile HTKTools/Makefile HLMLib/Makefile HLMTools/Makefile HTKLVRec/Makefile hts_engine/Makefile Makefile)

case "$enable_trad_htk" in
     yes)
cat<<EOF
**************************************************

HTK is now ready to be built.

Type "make all" to build and install the HTK
libraries and tools.

The tools will be installed in ${prefix}/`basename ${bindir}`

`echo $build_notes | fmt -w 50`

**************************************************
EOF
	;;
	*)
cat<<EOF
**************************************************

HTK is now ready to be built.

Type "make all" to build the HTK libraries
and tools.

Then "make install" to install them.

The tools will be installed in ${prefix}/`basename ${bindir}`

`echo $build_notes | fmt -w 50`

**************************************************
EOF
	;;
esac
