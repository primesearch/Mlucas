#!/usr/bin/env perl
#                                      Hey, EMACS: -*- cperl -*-
# manpage.pl - script to generate mlucas(1) man page from embedded Pod
# Copyright (C) 2015 Alex Vong <alexvong1995@gmail.com>,
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# First parameter, NAME, should be all caps
# Second parameter, SECTION, should be 1-8, maybe w/ subsection
# other parameters are allowed: see man(7), man(1)
#
# Some roff macros, for reference:
# .nh        disable hyphenation
# .hy        enable hyphenation
# .ad l      left justify
# .ad b      justify to both left and right margins
# .nf        disable filling
# .fi        enable filling
# .br        insert line break
# .sp <n>    insert n+1 empty lines
# for manpage-specific macros, see man(7)
# for reference for writing man pages, see man-page(7).
# for Perl Pod format, see perlpod(1).

use strict;
use warnings;

use Pod::Man;

# my global variables
my $ERR_USAGE = 67;
my $PROG_NAME = "MLUCAS";
my $PROG_VER = "Mlucas 14.1";
my $SECTION_NAME = "User Commands";
my $SECTION_NUM = "1";
my $LAST_UPDATE = "2015-08-04";
my $ERR_HANDLING = "die";
chomp(my $USAGE = <<EOF);
Usage: $0 [-h | --help | --m4 | --pmp]

  -h, --help  show this message
  --m4        use m4 preprocessor (default)
  --pmp       use poor man's preprocessor (use if m4 is not found)
EOF
chomp(my $GREET = <<'EOF');
.\"                                      Hey, EMACS: -*- nroff -*-
EOF
chomp(my $INTRO = <<'EOF');
.\" mlucas.1 - man page of mlucas written for Debian GNU/Linux
.\" Copyright (C) 2015 Alex Vong <alexvong1995@gmail.com>,
.\"
.\" This program is free software; you can redistribute it and/or modify
.\" it under the terms of the GNU General Public License as published by
.\" the Free Software Foundation; either version 2 of the License, or
.\" (at your option) any later version.
.\"
.\" This program is distributed in the hope that it will be useful,
.\" but WITHOUT ANY WARRANTY; without even the implied warranty of
.\" MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
.\" GNU General Public License for more details.
.\"
.\" You should have received a copy of the GNU General Public License along
.\" with this program; if not, write to the Free Software Foundation, Inc.,
.\" 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
.\"
.\" First parameter, NAME, should be all caps
.\" Second parameter, SECTION, should be 1-8, maybe w/ subsection
.\" other parameters are allowed: see man(7), man(1)
.\"
.\" Some roff macros, for reference:
.\" .nh        disable hyphenation
.\" .hy        enable hyphenation
.\" .ad l      left justify
.\" .ad b      justify to both left and right margins
.\" .nf        disable filling
.\" .fi        enable filling
.\" .br        insert line break
.\" .sp <n>    insert n+1 empty lines
.\" for manpage-specific macros, see man(7)
.\" for reference for writing man pages, see man-page(7).
.\" for Perl Pod format, see perlpod(1).
EOF

# poor man's preprocessor
sub pmp
{
  my $tmp;
  my %def;

  $_ = join "", @_;

  # build hash table of macro definitions
  while (/define\(\s*\[\+(.+?)\+\]\s*,
	  \s*\[\+(.+?)\+\]\s*\)dnl\n
	 /sx)
    {
      s/define\(\s*\[\+(.+?)\+\]\s*,
	\s*\[\+(.+?)\+\]\s*\)dnl\n
       //sx;
      $def{$1} = "$2";
    }

  # remove blank lines starting with dnl
  s/\ndnl\n/\n/gs;

  # implement _INDENT_BLOCK macro
  while(/_INDENT_BLOCK\(\s*\[\+(.+?)\+\]\s*,
	 \s*\[\+(.+?)\+\]\s*\)
	/sx)
    {
      my ($grp1, $grp2) = ("$1", "$2");

      # Use _INDENT_BLOCK_NORMAL if $1 does not contain B<>
      if ($grp1 !~ /B<[^<]+>\([012345678]?\)/)
	{
	  s/_INDENT_BLOCK\(\s*\[\+.+?\+\]\s*,
	    \s*\[\+.+?\+\]\s*\)
	   /$def{_INDENT_BLOCK_NORMAL}/sx;
	}
      # Otherwise, use _INDENT_BLOCK_HACKISH
      else
	{
	  s/_INDENT_BLOCK\(\s*\[\+.+?\+\]\s*,
	    \s*\[\+.+?\+\]\s*\)
	   /$def{_INDENT_BLOCK_HACKISH}/sx;
	}

      s/\$1/$grp1/gs;
      s/\$2/$grp2/gs;
    }

  # implement _START_NO_BLANK_LINES and _END_NO_BLANK_LINES
  s/\n_START_NO_BLANK_LINES\n
   /\n$def{_START_NO_BLANK_LINES}\n/gsx;
  s/\n_END_NO_BLANK_LINES\n
   /\n$def{_END_NO_BLANK_LINES}\n/gsx;

  # implement _JOIN
  while (/_JOIN\(\s*\[\+([^_]+)\+\]\s*\)/s)
    {
      $tmp = "$1";
      $tmp =~ s/\+\]\s*,\s*\[\+//gs;
      s/_JOIN\(([^_]+)\)/$tmp/;
    }

  # remove [++] style quotes
  s/\[\+|\+\]//gs;

  return split "$/";
}

# preprocess embedded Pod in this script using M4
# or pmp if M4 is not available
sub preprocess
{
  $_ = shift;
  if ($_ ne "--pmp")
    {
      return <<`EOF`;
printf "change""quote(\\\`[+', \\\`+]')" | cat - "$0" | /usr/bin/env m4
EOF
    }
  else
    {
      open my $self, "< $0" or die "cannot open < $0 $!";
      my @lines = pmp <$self>;
      close $self;
      return @lines;
    }
}

# convert Pod to Man using Pod::Man
sub pod2man
{
  my $output;
  my $parser = Pod::Man->new
    (
     center => "$SECTION_NAME",
     date => "$LAST_UPDATE",
     errors => "$ERR_HANDLING",
     name => "$PROG_NAME",
     release => "$PROG_VER",
     section => "$SECTION_NUM",
    );

  $parser->output_string(\$output);
  $parser->parse_lines(@_);
  return split "$/", "$output";
}

# add greeting and intorduction in the beginning
sub add_greeting_intro
{
  my $header = shift;
  my $blank_line = shift;

  unshift @_, "$GREET", "$header", "$blank_line", "$INTRO", "$blank_line";
  return @_;
}

# postprocess Man using perlre
sub postprocess
{
  foreach (@_)
    {
      # convert .IP macro to .TP macro
      s/^\.(el|ie\sn)\s\.IP\s"([^"]+)"\s(\d+)
       /.$1 \\{\\PERL_NEWLINE.TP $3PERL_NEWLINE$2\\}/gx;
      s/^\.IP\s"([^"]+)"\s(\d+)
       /.TP $2PERL_NEWLINE$1/gx;

      # convert `\fBfunction\fR()' to `.BR function ()'
      # convert `\fBfunction\fR(section)' to `.BR function (section)'
      s/\\fB([^\\]+)\\fR\((\d?)\)([[:punct:]]?)\s*(\.?[^\s[:punct:]]+)
       /PERL_NEWLINE.BR $1 ($2)$3PERL_NEWLINE$4/gx;
      s/\\fB([^\\]+)\\fR\((\d?)\)([[:punct:]]?)\s*
       /PERL_NEWLINE.BR $1 ($2)$3/gx;

      # remove \& followed by newline
      s/\\&PERL_NEWLINE//g;

      # convert dummy newline to real one in order to avoid a bug
      # where newline `\n' fail to render
      # when being placed between a comma `,' and a dot `.'
      s/PERL_NEWLINE/\n/g;

      print;
    }
}

# main
$\ = "$/";
foreach (@ARGV)
{
  if ($_ eq "-h" or $_ eq "--help")
    {
      print "$USAGE";
      exit;
    }
  elsif ($_ eq "--pmp")
    {
      postprocess(add_greeting_intro(pod2man(preprocess "$_")));
      exit;
    }
  elsif ($_ eq "--m4")
    {
      postprocess(add_greeting_intro(pod2man(preprocess "$_")));
      exit;
    }
  else
    {
      print STDERR "invalid option: `$_'";
      print STDERR "$USAGE";
      exit "$ERR_USAGE";
    }
}
unless (system "/usr/bin/env which m4 > /dev/null")
  {
    postprocess(add_greeting_intro(pod2man(preprocess "--m4")));
    exit;
  }
else
  {
    postprocess(add_greeting_intro(pod2man(preprocess "--pmp")));
    exit;
  }

=pod

=for comment
[+start removing blank lines+]
dnl
define([+_START_NO_BLANK_LINES+],
[+=for man .PD 0

=begin man

.de Sp
.br
..

=end man+])dnl

=for comment
[+stop removing blank lines+]
dnl
define([+_END_NO_BLANK_LINES+],
[+=for man .PD

=begin man

.de Sp
.if t .sp .5v
.if n .sp
..

=end man+])dnl

=for comment
[+front end to invoke _INDENT_BLOCK_NORMAL and _INDENT_BLOCK_HACKISH
_INDENT_BLOCK_NORMAL is used
if no function or manpage references are present+]
dnl
define([+_INDENT_BLOCK+],
[+ifelse(regexp([+$1+],[+B<[^<]+>([012345678]?)+]),
[+-1+],
[+_INDENT_BLOCK_NORMAL([+$1+], [+$2+])+],
[+_INDENT_BLOCK_HACKISH([+$1+], [+$2+])+])+])dnl

=for comment
[+primary way to start indentation of text block
a text block is a piece of text without blank lines+]
dnl
define([+_INDENT_BLOCK_NORMAL+],
[+=over 4

=item [+$1+]

=for man .Sp

_START_NO_BLANK_LINES

[+$2+]

_END_NO_BLANK_LINES

=back+])dnl

=for comment
[+secondary way to start indentation of text block+]
dnl
define([+_INDENT_BLOCK_HACKISH+],
[+=for man .RE

=for man .RS 7

[+$1+]

=for man .RS 4

[+$2+]

=for man .RE+])dnl

=for comment
[+join strings together (e.g., `_JOIN([+a+], [+b+], [+c+])' => `abc')+]
dnl
define([+_JOIN+],
[+ifelse([+$#+],
[+1+],
[+$1+],
[+[+$1+]_JOIN(shift($@))+])+])dnl

=head1 NAME

mlucas - program to perform Lucas-Lehmer test on a Mersenne number,
2 ^ p - 1

=head1 SYNOPSIS

_START_NO_BLANK_LINES

B<mlucas>

B<mlucas -h>

B<mlucas> B<-s> B<tiny> | B<t> | B<small> | B<s> | B<medium> | B<m> |
B<large> | B<l> | B<huge> | B<h> | B<all> | B<a>
[B<-iters> B<100> | B<1000> | B<10000> [B<-nthread> I<threads>]]

B<mlucas> B<-m> I<exponent> | B<-f> I<exponent>
[B<-iters> B<100> | B<1000> | B<10000> [B<-nthread> I<threads>]]

B<mlucas> B<-fftlen> I<fft_length>
[B<-radset> I<radix_set>]
[B<-m> I<exponent> | B<-f> I<exponent>]
B<-iters> B<100> | B<1000> | B<10000>
[B<-nthread> I<threads>]

_END_NO_BLANK_LINES

=head1 DESCRIPTION

This manual page documents briefly the B<mlucas> command.

B<mlucas> is an open-source (and free/libre) program
for performing Lucas-Lehmer test on prime-exponent Mersenne numbers,
that is, integers of the form 2 ^ p - 1, with prime exponent p.
In short, everything you need to search for world-record Mersenne primes!
It has been used in the verification of various Mersenne primes,
including the 45th, 46th and 48th found Mersenne prime.

You may use it to test any suitable number as you wish,
but it is preferable that you do so in a coordinated fashion,
as part of the B<Great Internet Mersenne Prime Search> (B<GIMPS>).
For more information on B<GIMPS>,
see the B<Great Internet Mersenne Prime Search> subsection
within the B<NOTES> section and B<SEE ALSO> section.
Note that B<mlucas> is not (yet) as efficient as the main B<GIMPS> client,
George Woltman's B<Prime95> program
(a.k.a. B<mprime> for the (gnu/)linux version),
but that program is not truly open-source (and free/libre), since
it requires the user to abide by the prize-sharing rules set by its author
(incompatible with I<freedom to run the program as you wish,
for any purpose>),
should a user be lucky enough to find a new prime eligible for
one of the monetary prizes offered by the Electronic Freedom Foundation
(see L<EFF Cooperative Computing Awards|https://www.eff.org/awards/coop>
for details).

B<mlucas> reads the exponents from the F<$MLUCAS_PATH/worktodo.ini> file.
Results are written to the F<$MLUCAS_PATH/results.txt> file
and the exponent-specific F<$MLUCAS_PATH/*.stat> file
(see section B<FILES> for details).
Error messages are written to I<stderr>
and the F<$MLUCAS_PATH/*.stat> file.
Exponents can also be passed as command-line arguments
but this is mainly used for debugging (see section B<OPTIONS> for details).
In addition, B<mlucas> can perform the Pe'pin primality test
on Fermat numbers 2 ^ (2 ^ n) + 1,
using an exponent-optimized fast-transform length
much like that used for testing Mersenne numbers.

New users are urged to jump straight to the B<EXAMPLE> section
and follow the examples and pointers to other sections.
Users with little time for in-depth reading
should at least read the B<NOTES>, B<BUGS> and B<EXAMPLE> sections
for a brief introduction to the B<Great Internet Mersenne Prime Search>,
undesirable restrictions and common usages.
B<FILES> section is also highly recommended
since it describes the B<mlucas> configuration files
used for host-specific optimization and other B<mlucas>-generated files.
Advanced users should also peruse the B<OPTIONS> section
since it introduces less-commonly-used advanced options.
Experienced users who find this manual inadequate
should consult the B<SEE ALSO> section for further information.
Lastly, the F<Mlucas README>, available both online and offline,
is highly recommended
since it is written and maintained by the author of B<mlucas>
and should be considered the final authority.

=head1 OPTIONS

B<mlucas> follows the traditional POSIX (see B<standards>(7) for details)
command line syntax,
with short options starting with one dashes (`I<->').
A summary of options is included below.
A complete description is in the B<SEE ALSO> section.

=over 7

=item B<-h>

Show version of program and summary of options.

=item B<-s t, -s tiny>

Run 100-iteration self-test on a set of 32 Mersenne exponents,
ranging from 173431 to 2455003.
This will take around 1 minute on a fast (pre-2010) CPU.

=item B<-s s, -s small>

Run 100-iteration self-test on a set of 24 Mersenne exponents,
ranging from 173431 to 1245877.
This will take around 10 minutes on a fast (pre-2010) CPU.

=item B<-s m, -s medium>

Run 100-iteration self-test on a set of 24 Mersenne exponents,
ranging from 1327099 to 9530803.
This will take around an hour on a fast (pre-2010) CPU.

=item B<-s l, -s large>

Run 100-iteration self-test on a set of 24 Mersenne exponents,
ranging from 10151971 to 72851621.
This will take around an hour on a fast (pre-2010) CPU.

=item B<-s h, -s huge>

Run 100-iteration self-test on a set of 16 Mersenne exponents,
ranging from 77597293 to 282508657.
This will take a couple of hours on a fast (pre-2010) CPU.

=item B<-s a, -s all>

Run 100-iteration self-test on all Mersenne exponents
and all FFT radix sets.
This will take several hours on a fast (pre-2010) CPU.

=item B<-fftlen> I<fft_length>

This allows the user to specify the length of the fast-transform (FFT)
used to effect the large-integer modular multiply
which is at the heart of all such nonfactorial primality tests.
The length unit here is in terms of the number of double-precision
machine words used in the multiword-integer encoding
of the primality test residue
which is both input and result of each of said multiplies.
Because mlucas is intended for testing numbers with many millions of bits,
we generally speak of these FFT lengths in terms of kilodoubles
(= 2 ^ 10 or 1024 doubles).
If I<fft_length> is one of the available FFT lengths (in kilodoubles),
run all available FFT radices available at that length,
unless the I<-radset> flag is also invoked (see below for details).
If I<-fftlen> is invoked with either the I<-m> or I<-f> flag,
the self-tests will perform the first 100 iterations
of a Lucas-Lehmer test (I<-m>) or Pe'pin test (I<-f>)
on the user-specified Mersenne or Fermat number.
If no user-set exponent is invoked,
do 100 Lucas-Lehmer test iterations using
the default self-test Mersenne or Fermat exponent for that FFT length.
The program uses this to find the optimal radix set for a given FFT length
on your hardware.

=item B<-iters> B<100> | B<1000> | B<10000>

Do I<100>, I<1000> or I<10000> self-test iterations of the type determined
by the modulus-related options
(I<-s> / I<-m> = Lucas-Lehmer test iterations with initial seed 4,
I<-f> = Pe'pin test squarings with initial seed 3).
Default is I<100> iterations.

=item B<-radset> I<radix_set>

Specify index of a set of complex FFT radices to use,
based on the big selection table in the function B<get_fft_radices>().
This requires a supported value of I<-fftlen> to be specified,
meaning (for an FFT length supported by the program)
an index B<0>, B<1>, B<2>, ... and so on.
B<0> is always a valid radix set index;
how high one can go in the enumeration depends on the FFT length.
As soon as the user tries an index out of range of the current FFT length,
the program will error-exit with an informational message to that effect,
which also notes the maximum allowable radix set index for that FFT length.

=item B<-nthread> I<threads>

For multithread-enabled (default) build,
perform the test in parallel mode with this many threads.

=item B<-m> I<exponent>

Perform a Lucas-Lehmer primality test of the Mersenne number
M(I<exponent>) = 2 ^ I<exponent> - 1,
where I<exponent> must be an odd prime.
If I<-iters> is also invoked, this indicates a timing test.
This requires suitable added arguments
(I<-fftlen> and, optionally, I<-radset>) to be supplied.
If the I<-fftlen> option (and optionally I<-radset>) is also invoked
but I<-iters> is not,
the program first checks
the first line of the F<$MLUCAS_PATH/worktodo.ini> file to see
if the assignment specified there is a Lucas-Lehmer test
with the same exponent as specified via the I<-m> argument.
If so, the I<-fftlen> argument is treated as a user override
of the default FFT length for the exponent.
If I<-radset> is also invoked, this is similarly treated as
a user-specified radix set for the user-set FFT length;
otherwise the program will use the F<$MLUCAS_PATH/mlucas.cfg> file
to select the radix set to be used for the user-forced FFT length.
If the F<$MLUCAS_PATH/worktodo.ini> file entry
does not match the I<-m> value,
a set of timing self-tests is run on the user-specified Mersenne number
using all sets of FFT radices available at the specified FFT length.
If the I<-fftlen> option is not invoked,
the tests use all sets of FFT radices
available at that exponent's default FFT length.
Use this to find the optimal radix set
for a single given Mersenne exponent on your hardware,
similarly to the I<-fftlen> option.
Perform 100 iterations, or as many as specified via the I<-iters> flag.

=item B<-f> I<exponent>

Perform a base-3 Pe'pin test on the Fermat number
F(I<exponent>) = 2 ^ (2 ^ I<exponent>) + 1.
If desired this can be invoked together with the I<-fftlen> option
as for the Mersenne-number self-tests (see above notes on the I<-m> flag;
note that not all FFT lengths supported for I<-m> are available for I<-f>:
I<-m> permits FFT lengths of form I<odd> * 2 ^ n
with I<odd> = any of B<1>, B<3>, B<5>, B<7>, B<9>, B<11>, B<13>, B<15>;
I<-f> allows odd = B<1>, B<7>, B<15> and B<63>)
Optimal radix sets and timings
are written to the F<$MLUCAS_PATH/fermat.cfg> file.
Perform 100 iterations, or as many as specified via the I<-iters> flag.

=back

=head1 EXIT STATUS

_INDENT_BLOCK(

[+The list of exit status values is limited.
It is not possible to determine the cause of failure
from the exit status value alone.
However, B<mlucas> make use of I<stderr> to print error messages
as well as saving them to the F<$MLUCAS_PATH/*.stat> file,
where I<*> is in the form+],

[+pI<exponent>+])

_INDENT_BLOCK(

[+for Mersenne number 2 ^ I<exponent> - 1 or+],

[+fI<exponent>+])

for Fermat number 2 ^ (2 ^ I<exponent>) + 1.
(see B<FILES> section for details).

=over 7

=item B<0>

Exit successfully.

=item B<1>

_START_NO_BLANK_LINES

Assertion failure.

Cannot determine the number of CPUs.

Unknown fetal error.

Radix set index not available for given FFT length.

_END_NO_BLANK_LINES

=item B<255>

_START_NO_BLANK_LINES

B<thread_policy_set>() failure.

B<malloc>(3), B<calloc>(3) or B<realloc>(3) failure.

B<pthread_create>(3) or B<pthread_join>(3) failure.

_END_NO_BLANK_LINES

=back

=head1 ENVIRONMENT

B<mlucas> honors the following environment variables, if they exist:

=over 7

=item B<MLUCAS_PATH>

The path to read B<mlucas> configuration files
and to write B<mlucas> generated files (see B<FILES> section for details).
B<MLUCAS_PATH> must end with a slash (e.g., I</home/foolish/bar/>.
If B<MLUCAS_PATH> is not set,
then B<MLUCAS_PATH> defaults to I<$HOME/.mlucas.d/>,
where the environmental variable I<$HOME> will be expanded
in the environment where B<mlucas> is invoked.
B<mlucas> will attept to make the directory with parents
pointed by B<MLUCAS_PATH> using the B<mkdir>(1) command.
The effect is similar to executing I<mkdir -p $MLUCAS_PATH> in the shell
provided that the I<-p> flag is honored.

=back

=head1 FILES

This section details B<mlucas> configuration files
and B<mlucas> generated files.
As noted in the B<ENVIRONMENT> section,
B<$MLUCAS_PATH> defaults to I<$HOME/mlucas.d/> but this can be
overridden at run-time by setting the B<MLUCAS_PATH> environment variable.

=over 7

=item B<$MLUCAS_PATH/*.stat>

_INDENT_BLOCK(

[+The filename-prefix wildcard I<*>
is as described in the EXIT STATUS section;
for the primality test of the Mersenne number 2 ^ I<exponent> - 1
it is of the form+],

[+pI<exponent>+])

_INDENT_BLOCK(

[+All important events, per-10000-iteration residues
(or per-100000-iteration if more than 4 threads are used for the test)
and the final residue during Lucas-Lehmer test of I<exponent>
are recorded in this file.
It can be seen as an I<exponent>-specific detailed
F<$MLUCAS_PATH/results.txt>
(see F<$MLUCAS_PATH/results.txt> below for details).
This file is useful for debugging purposes.
Its format looks like:+],

[+B<INFO>: I<event>

...

B<M>I<exponent>: B<using FFT length> I<fft_length>B<K>
B<=> I<fft_length * 1024> B<8-byte floats.>
Bz<this gives an average> I<bits> B<bits per digit>

B<Using complex FFT radices> I<radix_set>
(product of all elements of radix_set = fft_length / 2)

...

B<[>I<date_and_time>B<]>
B<M>I<exponent>
B<Iter#> B<=> I<iterations>
B<clocks> B<=> I<time_taken_per_10000_iterations>
B<[>I<   time_taken_per_iteration> B<sec/iter>B<]>
B<Res64:> I<residue>B<.>
B<AvgMaxErr> B<=> I<roe_avg>B<.>
B<MaxErr> B<=> I<roe_max>

...

[B<Restarting> B<M>I<exponent> B<at iteration> B<=> I<iteration>B<.>
B<Res64:> I<residue>

B<M>I<exponent>: B<using FFT length> I<fft_length>B<K>
B<=> I<fft_length * 1024> B<8-byte floats.>

B<this gives an average> I<bits> B<bits per digit>

B<Using complex FFT radices> I<radix_set>]
(product of all elements of radix_set = fft_length / 2)

...

B<M>I<exponent> B<is not prime.>
B<Res64:> I<residue>B<.>
B<Program: E14.1>

B<M>I<exponent> B<mod 2^36>     B<=>          I<remainder_1>

B<M>I<exponent> B<mod 2^35 - 1> B<=>          I<remainder_2>

B<M>I<exponent> B<mod 2^36 - 1> B<=>          I<remainder_3>+])

=item B<$MLUCAS_PATH/fermat.cfg>

The format of this file is exactly the same as
the format of F<$MLUCAS_PATH/mlucas.cfg>
(see F<$MLUCAS_PATH/mlucas.cfg> below for details).

=item B<$MLUCAS_PATH/mlucas.cfg>

_INDENT_BLOCK(

[+This file stores the radix set with best timing for each FFT length.
Its format looks like:+],

[+B<14.1>

I<fft_length> B<msec/iter> B<=> I<timing>
B<ROE[avg,max]> B<=> B<[>I<roe_avg>B<,> I<roe_max>B<]>
B<radices> B<=> I<radix_set>

...+])

_END_NO_BLANK_LINES

Normally, the I<timing> entry for each line should be monotonic
from above to below since larger FFT length should take longer to test.
But it is OK for a given I<fft_length> to have a higher I<timing> than
the one after it since B<mlucas> checks the timings listed in this file
for all FFT lengths >= the default FFT length for the number being tested,
and uses the FFT length having the smallest listed timing.
However, if you notice that this file has any entries such that
a given I<fft_length> has a timing 5% or more greater than the next-larger
FFT length, or higher timing than two or more larger FFT lengths,
please contact the author (see B<BUGS> section for details).

=item B<$MLUCAS_PATH/nthreads.ini>

This file sets the number of threads used.
It should only contain a positive integer
since the content of this file is read by
B<sscanf(>I<in_line>, I<"%d">B<,> I<&NTHREADS>B<);>
where the variable I<in_line>
contains the content of the F<$MLUCAS_PATH/nthreads.ini> file.
If this file is not present,
B<mlucas> will use as many threads as the number of CPUs detected.
The number of threads used set by this file
can be overridden by setting I<-nthread> flag at run-time.
This file is for those who want to set the number of threads
to be greater or less than the number of CPUs detected.
This can be useful since some users reported up to 10% performance gain
when using more threads than the number of CPUs detected.

=item B<$MLUCAS_PATH/results.txt>

_INDENT_BLOCK(

[+Important events which occurred during Lucas-Lehmer test
and the final residue obtained are recorded in this file.
This file summarizes important information
in all F<$MLUCAS_PATH/*.stat> files
(see F<$MLUCAS_PATH/*.stat> above for details) into a single file.
This file (more specifically, any results
which were added to it since your last checkin from)
should be submitted to the PrimeNet server (see subsection
B<Great Internet Mersenne Prime Search> in section B<NOTES> for details)
since the Lucas-Lehmer test exponents are obtained from the PrimeNet server
(see F<$MLUCAS_PATH/worktodo.ini> below for details).
Its format looks like:+],

[+B<INFO:> I<event>

...

[B<M>I<exponent> B<Roundoff warning on iteration> I<iteration>B<,>
B<maxerr> B<=> I<roundoff_error>

B< Retrying iteration interval to see if roundoff error is reproducible.>

[B<Retry of iteration interval with fatal roundoff error was successful.>]]

...

B<M>I<exponent>B< is not prime.>
B<Res64:> I<residue>B<.>
B<Program: E14.1>

B<M>I<exponent> B<mod 2^36>     B<=>          I<remainder_1>

B<M>I<exponent> B<mod 2^35 - 1> B<=>          I<remainder_2>

B<M>I<exponent> B<mod 2^36 - 1> B<=>          I<remainder_3>

...+])

=item B<$MLUCAS_PATH/worktodo.ini>

_INDENT_BLOCK(

[+This file contains Lucas-Lehmer test assignments to be tested.
Its format looks like:+],

[+I<assignment>B<=>I<ID>B<,>I<exponent>B<,>I<trial
factored up to>B<,>I<has P-1 factoring>

...+])

The I<assignment> field contains B<Test>
if the assignment is a first-time Lucas-Lehmer test,
or B<DoubleCheck> if the assignment is a double-check Lucas-Lehmer test.
(The program handles both cases the same way.)

_START_NO_BLANK_LINES

I<ID> is a unique 32-digit hex number.

I<exponent> specifies the Mersenne number
(of the form 2 ^ I<exponent> - 1) to be tested.

I<trial factored up to> is the number of bit this Mersenne number
has been trial factored up to without finding a factor.

I<has P-1 factoring> B<=> B<0> if no prior P-1 factoring has been done,
B<=> B<1> if P-1 factoring (without finding a factor) has been done.
Since mlucas currently has no P-1 factoring capability
it simply discards these data,
but users should prefer B<=> B<1> here
since such an assignment is slightly more likely (5-10%) to yield a prime.

_END_NO_BLANK_LINES

To do Lucas-Lehmer test,
you should reserve exponents from the PrimeNet server
and copy lines in the above format into the
F<$MLUCAS_PATH/worktodo.ini> file (see subsection
B<Great Internet Mersenne Prime Search> in section B<NOTES> for details).
You may need to create the F<$MLUCAS_PATH/worktodo.ini> file
if it does not exist.

=item B<Save files in $MLUCAS_PATH>

All files matching the following extended regular expression
(see B<regex>(7) for details)
in I<$MLUCAS_PATH> directory are save files:

    ^[fpq][0123456789]+([.][0123456789]+0M)?$

For both of the supported test types,
duplicate pairs of savefiles are written at each checkpoint,
to guard against corruption of the on-disk savefiles.
Lucas-Lehmer test savefile-pair names start with <p> and <q>, respectively,
while Pe'pin test savefile-pair names start with <f> and <q>, respectively.
They should not be modified but backups may be made by the user.
By default, the program will save a persistent backup of the primary
(B<p> or B<f>) save file every 10 millionth iteration,
for examples upon completion of the Lucas-Lehmer test of M57885161
the user will find the following exponent-associated files
in the I<$MLUCAS_PATH> directory:

    p57885161.stat
    p57885161.10M
    p57885161.20M
    p57885161.30M
    p57885161.40M
    p57885161.50M

=back

=head1 NOTES

=head2 Great Internet Mersenne Prime Search

This subsection needs to be compeleted...

=head1 BUGS

The argument parser is buggy.
The relative position of arguments is relevant to B<mlucas>,
the order of arguments in B<SYNOPSIS>
should be followed to avoid confusing the parser.
Only I<100>, I<1000> and I<10000> are supported for I<-iters> flag.
However, the parser will not reject unsupported arguments.
Using unsupported arguments for I<-iters> flag
may trigger strange behaviour.

Sometimes there is more than one applicable exit status values
(see B<EXIT STATUS> section for details).
In such case, there is no guarantee which will be returned.
For example,
if B<malloc>(3) failure triggers an assertion failure.
It is possible that B<mlucas> returns I<1>
instead of I<255> as exit status value.

For problems regarding the program B<mlucas>,
please contact the author Ernst W. Mayer <ewmayer@aol.com>.
For installation and documentation related problems
regarding the Debian package and this manual,
please use B<reportbug>(1) to contact Alex Vong <alexvong1995@gmail.com>.

=head1 EXAMPLE

There are 3 common cases where you will want to run this program.
Normally, you should do a spot-check first to quick-test your build,
followed by the self-test range for `medium' exponents.
Finally, full-blown Lucas-Lehmer testing
which is the main purpose of this program.

=over 7

=item B<mlucas -fftlen 192 -iters 100 -radset 0 -nthread 2>

Perform spot-check to see if B<mlucas> works
and fill-in a bug report if it does not.
The spot check should produce residues
matching the internal tabulated ones.
If the residues does not match,
B<mlucas> should emit a verbose error message.

=item B<mlucas -s m>

Perform timing self-test for `medium' exponents
to tune code parameters for your platform.
Ordinary users are recommended to do this self-test only.
For best results,
run any self-tests under zero- or constant-load conditions.
The self-tests append
(or create if F<$MLUCAS_PATH/mlucas.cfg> does not exist)
new timing data to the F<$MLUCAS_PATH/mlucas.cfg>
(see B<FILES> section for details).
Before doing any self-tests,
you should first check if there is an existing
F<$MLUCAS_PATH/mlucas.cfg> file
and either delete it or do a backup-via-rename to
to prevent mixing old and new timing data.
F<$MLUCAS_PATH/mlucas.cfg> normally locates at
I<$HOME/.mlucas.d/> directory
although this can be overridden at run-time
by settingthe B<MLUCAS_PATH> environment variable
(see B<ENVIRONMENT> section for details).

=item B<mlucas &>

_INDENT_BLOCK(

[+Perform Lucas-Lehmer test on Mersenne numbers
by running B<mlucas> as a background job
(see B<JOB CONTROL> section in B<bash>(1)
and B<Builtins> subsection in B<dash>(1) for details).
To perform Lucas-Lehmer test on a given Mersenne number,
you must first perform a self-test
for `medium' exponents mentioned above,
or if you only desire to test a single selected Mersenne number,
a self-test for the default FFT length for that number:+],

[+mlucas -m I<exponent> -iters 100+])

In the case of multi-exponent "production testing",
you should reserve exponent from the PrimeNet server
and add them into F<$MLUCAS_PATH/worktodo.ini>
(see the subsection B<Great Internet Mersenne Prime Search>
within the section B<NOTES> and B<FILES> section for details).

=back

=head2 Advanced Usage Tips

To start B<mlucas> in terminal 1,
add the following lines to your login shell initialization file,
such as I<$HOME/.profile>
(see B<INVOCATION> section in B<bash>(1)
and B<Invocation> subsection B<dash>(1) for details).

    # Test if we are in tty1
    if test `tty` = '/dev/tty1'
    then
	# turn on job control
	set -m
	# start mlucas
	nice mlucas > /dev/null 2>&1 &
    fi

=head1 SEE ALSO

B<bash>(1), B<dash>(1), B<reportbug>(1)

L<http://hogranch.com/mayer/README.html>,
F</usr/share/doc/mlucas/html/README.html>

B<mlucas> is documented fully by F<Mlucas README>,
available both online and offline as shown above.

L<F<Great Internet Mersenne Prime Search>|http://www.mersenne.org>

L<F<Mersenne Forum>|http://www.mersenneforum.org>

_JOIN([+L<F<Chris Caldwell's web page on Mersenne numbers>|+],
[+http://www.utm.edu/research/primes/mersenne.shtml>+])

_JOIN([+L<F<Richard Crandall and Barry Fagin, +],
[+Discrete Weighted Transforms and Large-Integer Arithmetic.>|+],
[+http://www.faginfamily.net/barry/Papers/+],
[+Discrete%20Weighted%20Transforms.pdf>+])

_JOIN([+L<F<Richard E. Crandall, Ernst W. Mayer, +],
[+and Jason S. Papadopoulos, +],
[+The Twenty-Fourth Fermat Number is Composite.>|+],
[+http://hogranch.com/mayer/F24.pdf>+])

=cut
