[+ autogen5 template -*- Scheme -*- am +]
[+ #|
 scheme.tpl - scheme function definitions to be included by main.tpl
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

[+
 ;; generate file name converter
 ;; evaluate to a function
 (define (generate-name-converter text-to-be-appended)
   (define (append-path name)
     (string-substitute name
			" "
			text-to-be-appended))
   (lambda (. name-list)
     (apply string-append (map append-path name-list))))

 ;; append `$(srcdir)/src/' to source file name
 (define source-name->path-name
   (generate-name-converter " $(srcdir)/src/"))

 ;; append `$(srcdir)/html/' to html file name
 (define html-name->path-name
   (generate-name-converter " $(srcdir)/html/"))

 ;; append `$(srcdir)/doc/' to documentation file name
 (define documentation-name->path-name
   (generate-name-converter " $(srcdir)/doc/"))

 ;; append `$(srcdir)/patch/' to patch file name
 (define patch-name->path-name
   (generate-name-converter " $(srcdir)/patch/"))

 ;; append `$(srcdir)/scripts/' to script file name
 (define script-name->path-name
   (generate-name-converter " $(srcdir)/scripts/"))

 ;; append `$(srcdir)/am/' to autogen template and definition file name
 (define template-name->path-name
   (generate-name-converter " $(srcdir)/am/"))

  ;; append `$(srcdir)/COPYING.d/' to license file name
 (define license-name->path-name
   (generate-name-converter " $(srcdir)/COPYING.d/"))

 ;; append directory name to object file name
 (define (object-name->path-name dir-name . name-list)
   (define path (string-append " " dir-name))
   (define (append-path name)
     (string-substitute name
			" "
			path))
   (apply string-append (map append-path name-list)))

 ;; append `../' to the path name
 (define (append-prev-dir . name-list)
   (define (append-path name)
     (string-substitute name
			" "
			" ../"))
   (apply string-append (map append-path name-list)))

 ;; change `.c' extension to `.o' extension
 (define (source-name->object-name . name-list)
   (define (change-extension name)
     (string-substitute name
			".c"
			".o"))
   (apply string-append (map change-extension name-list)))

 ;; argument folding
 ;; "abc"-> "abc"
 ;; '("a" "b" "c") -> "a b c"
 (define (fold-arg arg)
   (if (not (list? arg))
       arg
       (string-join arg)))

 ;; list of arguments used by all rules
 (define arg-list
   '("?target-name?" "?target?" "?ingredient?" "?echo?" "?dir-name?"))

 ;; generic rule generator
 (define (generate-rule . args)
   (string-substitute
    (get "rule")
    (cons "?command?" arg-list)
    (map fold-arg args)))

 ;; compilation rule generator
 (define (generate-compilation-rule . args)
   (define full-arg-list (cons (get "compilation_rule") arg-list))
   (string-substitute
    (apply generate-rule full-arg-list)
    (cons "?ingredient-path-name?" (cons "?compilation-flag?" arg-list))
    (map fold-arg args)))

 ;; srcipt making rule generator
 (define (generate-script-making-rule . args)
   (define sed-scripts
     (string-join (map (lambda (variable)
			 (string-substitute (get "sed_script")
					    "?variable?"
					    variable))
		       (car args))))
   (define sed-command
     (string-substitute (get "make_script_rule")
			"?sed-scripts?"
			sed-scripts))
   (define full-arg-list (cons sed-command arg-list))
   (string-substitute
    (apply generate-rule full-arg-list)
    (cons "?ingredient-path-name?" arg-list)
    (map fold-arg (cdr args))))

 ;; generate automake conditional
 ;; (if-use-threads "I am " "using threads" "not using threads") ->
 ;; "if USE_THREADS
 ;; I am using threads
 ;; else
 ;; I am not using threads
 ;; endif"
 (define (if-use-threads . args)
   (define (am_conditional share use-threads do-not-use-threads)
     (string-append "if USE_THREADS\n"
		    share use-threads
		    "\nelse\n"
		    share do-not-use-threads
		    "\nendif"))
   (apply am_conditional (map fold-arg args)))

 ;; generate sse2, avx and avx2 rules from a single amd64 rule
 ;; (amd64-rules->sse2-avx-avx2-rules "capitalize amd64 will give AMD64") ->
 ;; "capitalize sse2 will give SSE2
 ;; capitalize avx will give AVX
 ;; capitalize avx2 will give AVX2"
 (define (amd64-rules->sse2-avx-avx2-rules rule)
   (define (make-rule instruction)
     (string-substitute rule
			(list "amd64" "AMD64")
			(list instruction (string-upcase instruction))))
   (string-join (map make-rule '("sse2" "avx" "avx2")) "\n")) +]
