[+ autogen5 template -*- Scheme -*- am +]
[+ #|
 amd64_instruction_unspecified.tpl - template to be included by main.tpl
 Copyright (C) 2015  Alex Vong

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software Foundation,
 Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.  |# +]

[+ #| What get cleaned by `$ make mostlyclean'  |# +]
[+ (if-use-threads
    (list
     "MOSTLYCLEANFILES = mlucas.cfg\n"
     (amd64-rules->sse2-avx-avx2-rules
      "MOSTLYCLEANFILES += amd64/*.o\n")
     (amd64-rules->sse2-avx-avx2-rules
      "MOSTLYCLEANFILES += AMD64_TRICKY_O.tmp AMD64_TRICKY_O.stamp\n"))
    (amd64-rules->sse2-avx-avx2-rules
     "MOSTLYCLEANFILES += AMD64_NORMAL_O-AMD64_THREADS_O.tmp \
AMD64_NORMAL_O-AMD64_THREADS_O.stamp")
    (amd64-rules->sse2-avx-avx2-rules
     "MOSTLYCLEANFILES += AMD64_NORMAL_O.tmp AMD64_NORMAL_O.stamp")) +]

[+ #| Recipe of mlucas (wrapper script)  |# +]
bin_SCRIPTS = mlucas

[+ (generate-script-making-rule
    '("pkglibexecdir")
    (script-name->path-name (get "script_in"))
    "mlucas"
    "mlucas"
    (script-name->path-name (get "script_in"))
    "echo '  GEN      mlucas'"
    "./") +]

[+ #| Recipe of sse2/mlucas, avx/mlucas and avx2/mlucas  |# +]
nobase_pkglibexec_PROGRAMS = sse2/mlucas avx/mlucas avx2/mlucas

[+ (amd64-rules->sse2-avx-avx2-rules "amd64_mlucas_SOURCES=") +]

[+ (amd64-rules->sse2-avx-avx2-rules
    "amd64_mlucas_CPPFLAGS = $(AM_CPPFLAGS) -DUSE_AMD64") +]

[+ (amd64-rules->sse2-avx-avx2-rules
    "amd64_mlucas_CFLAGS = $(AM_CFLAGS) -mamd64") +]

[+ (amd64-rules->sse2-avx-avx2-rules
    (if-use-threads
     (list "amd64_mlucas_LDADD="
	   (object-name->path-name
	    "amd64/"
	    (source-name->object-name (get "normal_c") (get "tricky_c"))))
     (list (object-name->path-name
	    "amd64/"
	    (source-name->object-name (get "threads_c")))
	   " -lpthread -lrt")
     "")) +]

[+ (amd64-rules->sse2-avx-avx2-rules
    (if-use-threads
     (list "amd64_mlucas_DEPENDENCIES="
	   (object-name->path-name
	    "amd64/"
	    (source-name->object-name (get "normal_c") (get "tricky_c"))))
     (object-name->path-name
      "amd64/"
      (source-name->object-name (get "threads_c")))
     "")) +]

[+ #| Recipe of object files  |# +]
[+ (amd64-rules->sse2-avx-avx2-rules
    (if-use-threads
     ""
     (generate-compilation-rule
      (append-prev-dir
       (source-name->path-name (get "normal_c") (get "threads_c")))
      "-DUSE_AMD64 -mamd64 $(NORMALCFLAGS)"
      "AMD64_NORMAL_O-AMD64_THREADS_O"
      (object-name->path-name
       "amd64/"
       (source-name->object-name (get "normal_c") (get "threads_c")))
      (source-name->path-name (get "normal_c") (get "threads_c"))
      "echo '  CC       $$AMD64_NORMAL_O $$AMD64_THREADS_O'"
      "amd64/")
     (generate-compilation-rule
      (append-prev-dir
       (source-name->path-name (get "normal_c")))
      "-DUSE_AMD64 -mamd64 $(NORMALCFLAGS)"
      "AMD64_NORMAL_O"
      (object-name->path-name
       "amd64/"
       (source-name->object-name (get "normal_c")))
      (source-name->path-name (get "normal_c"))
      "echo '  CC       $$AMD64_NORMAL_O'"
      "amd64/"))) +]

[+ (amd64-rules->sse2-avx-avx2-rules
    (generate-compilation-rule
     (append-prev-dir (source-name->path-name (get "tricky_c")))
     "-DUSE_AMD64 -mamd64 $(TRICKYCFLAGS)"
     "AMD64_TRICKY_O"
     (object-name->path-name
      "amd64/"
      (source-name->object-name (get "tricky_c")))
     (source-name->path-name (get "tricky_c"))
     "echo '  CC       $$AMD64_TRICKY_O'"
     "amd64/")) +]
