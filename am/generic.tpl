[+ autogen5 template -*- Scheme -*- am +]
[+ #|
 generic.tpl - template to be included by main.tpl
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
    "MOSTLYCLEANFILES = mlucas.cfg TRICKY_O.tmp TRICKY_O.stamp"
    " NORMAL_O-THREADS_O.tmp NORMAL_O-THREADS_O.stamp"
    " NORMAL_O.tmp NORMAL_O.stamp") +]

[+ #| Recipe of mlucas (binary)  |# +]
bin_PROGRAMS = mlucas

mlucas_SOURCES = [+ #| '()  |# +]
[+ (if-use-threads
    (list "mlucas_LDADD="
	  (source-name->object-name (get "normal_c") (get "tricky_c")))
    (list (source-name->object-name (get "threads_c"))
	  " -lpthread -lrt")
    "") +]

[+ (if-use-threads
    (list "mlucas_DEPENDENCIES="
	  (source-name->object-name (get "normal_c") (get "tricky_c")))
    (source-name->object-name (get "threads_c"))
    "") +]

[+ #| Recipe of object files  |# +]
[+ (if-use-threads
    ""
    (generate-compilation-rule
     (source-name->path-name (get "normal_c") (get "threads_c"))
     "$(NORMALCFLAGS)"
     "NORMAL_O-THREADS_O"
     (source-name->object-name (get "normal_c") (get "threads_c"))
     (source-name->path-name (get "normal_c") (get "threads_c"))
     "echo '  CC       $$NORMAL_O $$THREADS_O'"
     "./")
    (generate-compilation-rule
     (source-name->path-name (get "normal_c"))
     "$(NORMALCFLAGS)"
     "NORMAL_O"
     (source-name->object-name (get "normal_c"))
     (source-name->path-name (get "normal_c"))
     "echo '  CC       $$NORMAL_O'"
     "./")) +]

[+ (generate-compilation-rule
    (source-name->path-name (get "tricky_c"))
    "$(TRICKYCFLAGS)"
    "TRICKY_O"
    (source-name->object-name (get "tricky_c"))
    (source-name->path-name (get "tricky_c"))
    "echo '  CC       $$TRICKY_O'"
    "./") +]
