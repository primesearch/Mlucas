#                                               -*- Autoconf -*-
# configure.ac - template used by autoconf to generate the `configure' script
# Copyright (C) 2015  Alex Vong
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([Mlucas], [14.1], [ewmayer@aol.com], [mlucas],
[http://hogranch.com/mayer/README.html])
AM_INIT_AUTOMAKE([-Wall -Werror gnu])
AM_SILENT_RULES([yes])
AC_CONFIG_SRCDIR([src/rng_isaac.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
: ${CFLAGS=""}
AC_PROG_CC

# Checks for libraries.
AC_SEARCH_LIBS([ceil, log, pow, sqrt, sincos, floor, lrint, atan], [m], [],
[
AC_MSG_ERROR([cannot find ceil(), log(), pow(), sqrt(), sincos(), floor(), dnl
lrint() or atan() in libm])
]
)

# Checks for header files.
AC_CHECK_HEADERS([fenv.h limits.h mach/mach.h stddef.h stdlib.h string.h dnl
sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([clock_gettime gethrtime gettimeofday memset pow sqrt dnl
strerror strstr strtoul])

# Enable Large File Support, useful on 32-bit systems.
AC_CHECK_DECL([_LARGEFILE_SOURCE],
[AC_SUBST([LFSCPPFLAGS], [['-D_FILE_OFFSET_BITS=64']])],
[AC_SUBST([LFSCPPFLAGS], [['-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64']])]
)

# Notify the user the existence of custom options.
AC_ARG_ENABLE([foo],
AS_HELP_STRING([[
Custom Options:
(please note that all options mentioned below does NOT override one another)]],
[]))

# Adding options for user specified NORMAL-CFLAGS and TRICKY-CFLAGS.
AC_ARG_ENABLE([NORMAL-CFLAGS],
AS_HELP_STRING([[--enable-NORMAL-CFLAGS='FOO' or --disable-NORMAL-CFLAGS]],
[set CFLAGS for normal source files to FOO or nothing ['-Os']]))
AC_ARG_ENABLE([TRICKY-CFLAGS],
AS_HELP_STRING([[--enable-TRICKY-CFLAGS='FOO' or --disable-TRICKY-CFLAGS]],
[set CFLAGS for optimization-sensitive source files to FOO or] dnl
[nothing ['-O1']]))

# Check if the user set NORMAL-CFLAGS and TRICKY-CFLAGS.
AS_IF([test "x$enable_NORMAL_CFLAGS" = "x"],
[AC_SUBST([NORMALCFLAGS], [['-Os']])],
[test "x$enable_NORMAL_CFLAGS" = "xno"],
[AC_SUBST([NORMALCFLAGS], [[]])],
[AC_SUBST([NORMALCFLAGS], [[$enable_NORMAL_CFLAGS]])]
)
AS_IF([test "x$enable_TRICKY_CFLAGS" = "x"],
[AC_SUBST([TRICKYCFLAGS], [['-O1']])],
[test "x$enable_NORMAL_CFLAGS" = "xno"],
[AC_SUBST([NORMALCFLAGS], [[]])],
[AC_SUBST([TRICKYCFLAGS], [[$enable_TRICKY_CFLAGS]])]
)

# Notify the user CFLAGS, CPPFLAGS and LDFLAGS are honored.
AC_ARG_ENABLE([foo],
AS_HELP_STRING([[CFLAGS='FOO']],
[set CFLAGS for all source files to FOO ['']]))
AC_ARG_ENABLE([foo],
AS_HELP_STRING([[CPPFLAGS='FOO']],
[set CPPFLAGS for all source files to FOO ['']]))
AC_ARG_ENABLE([foo],
AS_HELP_STRING([[LDFLAGS='FOO']],
[set LDFLAGS for all linking to FOO ['']]))

# Notify the user --host is honored.
AC_ARG_ENABLE([foo],
AS_HELP_STRING([[--host='FOO']],
[skip autodetection and force detecting host as FOO ]
[(set --host=none to disable host-specific compilation flags)]))

# Add appropriate flag if host cpu is i386, i486, i586, i686 or amd64.
AC_CANONICAL_HOST
AS_IF([test "x$host_cpu" = "xi386" ||
test "x$host_cpu" = "xi486" ||
test "x$host_cpu" = "xi586" ||
test "x$host_cpu" = "xi686"],
[
AC_SUBST([ARCHCFLAGS], [['-m32']])
AC_SUBST([ARCHCPPFLAGS], [['-DCPU_IS_X86']])
]
)
AS_IF([test "x$host_cpu" = "xx86_64"],
[AC_SUBST([ARCHCFLAGS], [['-m64']])]
)

# Adding options for building with extended instruction sets support.
AC_ARG_ENABLE([instruction],
AS_HELP_STRING([[[--enable-instruction=[SSE2|AVX|AVX2] or]]] dnl
[[[--disable-instruction]]],
[build with SSE2, AVX, AVX2 or without extended instruction set]))

# Checks for extended instruction sets.
AS_IF([test "x$enable_instruction" = "xsse2" ||
test "x$enable_instruction" = "xSSE2"],
[
AC_SUBST([INSTRUCTIONCFLAGS], [["-msse2"]])
AC_SUBST([INSTRUCTIONCPPFLAGS], [["-DUSE_SSE2"]])
],
[test "x$enable_instruction" = "xavx" ||
test "x$enable_instruction" = "xAVX"],
[
AC_SUBST([INSTRUCTIONCFLAGS], [["-mavx"]])
AC_SUBST([INSTRUCTIONCPPFLAGS], [["-DUSE_AVX"]])
],
[test "x$enable_instruction" = "xavx2" ||
test "x$enable_instruction" = "xAVX2"],
[
AC_SUBST([INSTRUCTIONCFLAGS], [["-mavx2"]])
AC_SUBST([INSTRUCTIONCPPFLAGS], [["-DUSE_AVX2"]])
],
[test "x$host_cpu" = "xx86_64"],
[amd64_instruction_unspecified=true]
)

# Set AMD64_INSTRUCTION_UNSPECIFIED if amd64_instruction_unspecified is true.
AM_CONDITIONAL([AMD64_INSTRUCTION_UNSPECIFIED],
[test "x$amd64_instruction_unspecified" = "xtrue"])

# Adding options for disabling multithread-support.
AC_ARG_ENABLE([threads],
AS_HELP_STRING([[--disable-threads]],
[force building without multithread-support]))

# Set USE_THREADS if the user wants thread-support.
AM_CONDITIONAL([USE_THREADS],
[test "x$enable_threads" != "xno"])

# Adding options for setting MLUCAS_DEFAULT_PATH to "$HOME/.mlucas.d/".
AC_ARG_ENABLE([MLUCAS-DEFAULT-PATH],
AS_HELP_STRING([[--enable-MLUCAS-DEFAULT-PATH]],
[set MLUCAS_DEFAULT_PATH (unset by default) to "$HOME/.mlucas.d/"]))

# Check if the user set MLUCAS-DEFAULT-PATH.
AS_IF([test "x$enable_MLUCAS_DEFAULT_PATH" = "xyes"],
[AC_SUBST([PATHCPPFLAGS],
[[-DMLUCAS_DEFAULT_PATH="'"'"'"'""'"'$$'"'"HOME/.mlucas.d/"'"'"'"'"]])]
)

# Notify the user --disable-silent-rules is honored.
AC_ARG_ENABLE([foo],
AS_HELP_STRING([[--disable-silent-rules]],
[increase make verbosity]))

# Adding options for enabaling verbose compiler.
AC_ARG_ENABLE([verbose-compiler],
AS_HELP_STRING([[--enable-verbose-compiler]],
[enable verbose compiler]))

# Check if the user wants the compiler to be verbose.
AS_IF([test "x$enable_verbose_compiler" != "xyes"],
[AC_SUBST([REDIRECTION], ['2> /dev/null'])]
)

# Notify the user the end of custom options.
AC_ARG_ENABLE([foo],
AS_HELP_STRING([[
End of Custom Options]],
[]))

# End of configuration.
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
