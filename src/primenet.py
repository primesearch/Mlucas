#!/usr/bin/env python

# Automatic assignment handler for Mlucas using manual testing forms at mersenne.org

# EWM: adapted from https://github.com/MarkRose/primetools/blob/master/mfloop.py by teknohog and Mark Rose, with help rom Gord Palameta.

# This only handles LL testing (first-time and double-check) for now.
# To-do: Add support for trial factoring work.

# This version can run in parallel with Mlucas, as it uses lockfiles to avoid conflicts when updating files.

################################################################################
#                                                                              #
#   (C) 2017 by Ernst W. Mayer.                                                #
#                                                                              #
#  This program is free software; you can redistribute it and/or modify it     #
#  under the terms of the GNU General Public License as published by the       #
#  Free Software Foundation; either version 2 of the License, or (at your      #
#  option) any later version.                                                  #
#                                                                              #
#  This program is distributed in the hope that it will be useful, but WITHOUT #
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       #
#  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for   #
#  more details.                                                               #
#                                                                              #
#  You should have received a copy of the GNU General Public License along     #
#  with this program; see the file GPL.txt.  If not, you may view one at       #
#  http://www.fsf.org/licenses/licenses.html, or obtain one by writing to the  #
#  Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA     #
#  02111-1307, USA.                                                            #
#                                                                              #
################################################################################

import sys
import os.path
import re
from time import sleep
import os
import math
from optparse import OptionParser

# More python3-backward-incompatibility-breakage-related foo - thanks to Gord Palameta for the workaround:
#import cookielib
#import urllib
#import urllib2
try:
    import http.cookiejar as cookiejar
    from urllib.error import URLError
    from urllib.parse import urlencode
    from urllib.request import build_opener
    from urllib.request import HTTPCookieProcessor
except ImportError:
    import cookielib as cookiejar
    from urllib2 import URLError
    from urllib import urlencode
    from urllib2 import build_opener
    from urllib2 import HTTPCookieProcessor

primenet_baseurl = b"https://www.mersenne.org/"
primenet_login = False

def ass_generate(assignment):
	output = ""
	for key in assignment:
		output += key + "=" + assignment[key] + "&"
	#return output.rstrip("&")
	return output

def debug_print(text):
	if options.debug:
		print(progname + ": " + text)
		sys.stdout.flush()

def greplike(pattern, l):
	output = []
	for line in l:
		s = re.search(r"(" + pattern + ")$", line)
		if s:
			output.append(s.groups()[0])
	return output

def num_to_fetch(l, targetsize):
	num_existing = len(l)
	num_needed = targetsize - num_existing
	return max(num_needed, 0)

def readonly_file(filename):
	# Used when there is no intention to write the file back, so don't
	# check or write lockfiles. Also returns a single string, no list.
	if os.path.exists(filename):
		File = open(filename, "r")
		contents = File.read()
		File.close()
	else:
		contents = ""
	return contents

def read_list_file(filename):
	# Used when we plan to write the new version, so use locking
	lockfile = filename + ".lck"
	try:
		fd = os.open(lockfile, os.O_CREAT | os.O_EXCL)
		os.close(fd)
		if os.path.exists(filename):
			File = open(filename, "r")
			contents = File.readlines()
			File.close()
			return map(lambda x: x.rstrip(), contents)
		else:
			return []
	# This python2-style exception decl gives a syntax error in python3:
	# except OSError, e:
	# https://stackoverflow.com/questions/11285313/try-except-as-error-in-python-2-5-python-3-x
	# gives the fugly but portable-between-both-python2-and-python3 syntactical workaround:
	except OSError:
		_, e, _ = sys.exc_info()
		if e.errno == 17:
			return "locked"
		else:
			raise

def write_list_file(filename, l, mode="w"):
	# Assume we put the lock in upon reading the file, so we can
	# safely write the file and remove the lock
	lockfile = filename + ".lck"
	# A "null append" is meaningful, as we can call this to clear the
	# lockfile. In this case the main file need not be touched.
	if mode != "a" or len(l) > 0:
		content = "\n".join(l) + "\n"
		File = open(filename, mode)
		File.write(content)
		File.close()
	os.remove(lockfile)

def unlock_file(filename):
	lockfile = filename + ".lck"
	os.remove(lockfile)

def primenet_fetch(num_to_get):
	if not primenet_login:
		return []
	# <option value="100">Smallest available first-time tests
	# <option value="101">Double-check tests
	# <option value="102">World record tests
	assignment = {"cores": "1",
		"num_to_get": str(num_to_get),
		"pref": options.worktype,
		"exp_lo": "",
		"exp_hi": "",
	}
	try:
		r = primenet.open(primenet_baseurl + "manual_assignment/?" + ass_generate(assignment) + "B1=Get+Assignments")
		return greplike(workpattern, r.readlines())
	except URLError:
		debug_print("URL open error at primenet_fetch")
		return []

def get_assignment():
	w = read_list_file(workfile)
	if w == "locked":
		return "locked"

	tasks = greplike(workpattern, w)
	num_to_get = num_to_fetch(tasks, int(options.num_cache))

	if num_to_get < 1:
		debug_print(workfile + "already has >= " + str(len(tasks)) + " entries, not getting new work")
		# Must write something anyway to clear the lockfile
		new_tasks = []
	else:
		debug_print("Fetching " + str(num_to_get) + " assignments")
		new_tasks = primenet_fetch(num_to_get)

	write_list_file(workfile, new_tasks, "a")

def mersenne_find(line, complete=True):
	return re.search(r"Program:", line)

def submit_work():
	# Only submit completed work, i.e. the exponent must not exist in worktodo file any more
	files = [resultsfile, sentfile]
	rs = map(read_list_file, files)
	#
	# EWM: Mark Rose comments:
	# This code is calling the read_list_file function for every item in the files list. It's putting the
	# results of the function for the first file, resultsfile, in the first position in the array, rs[0].
	# Inside read_list_file, it's opening the file, calling readlines to get the contents of it into an array,
	# then calling the rstrip function on every line to remove trailing whitespace. It then returns the array.
	#
	# EWM: Note that read_list_file does not need the file(s) to exist - nonexistent files simply yield 0-length rs-array entries.

	if "locked" in rs:
		# Remove the lock in case one of these was unlocked at start
		for i in range(len(files)):
			if rs[i] != "locked":
				debug_print("Calling write_list_file() for" + files[i])
				write_list_file(files[i], [], "a")
		return "locked"

	results = rs[0]
	results = filter(mersenne_find, results)	# remove nonsubmittable lines from list of possibles
	results_send = [line for line in results if line not in rs[1]]	# if a line was previously submitted, discard
	results_send = list(set(results_send))	# In case resultsfile contained duplicate lines for some reason
	debug_print("New-results has " + str(len(results_send)) + " entries.")

	# Only for new results, to be appended to results_sent
	sent = []

	if len(results_send) == 0:
		debug_print("No complete results found to send.")
		# Don't just return here, files are still locked...
	else:
		# EWM: Switch to one-result-line-at-a-time submission to support error-message-on-submit handling:
		for sendline in results_send:
			debug_print("Submitting\n" + sendline)
			try:
				post_data = urlencode({"data": sendline})
				r = primenet.open(primenet_baseurl + "manual_result/default.php", post_data)
				res = r.read()
				if "Error code" in res:
					ibeg = res.find("Error code")
					iend = res.find("</div>", ibeg)
					print("Submission failed: '" + res[ibeg:iend] + "'")
				elif "Accepted" in res:
					sent += sendline
				else:
					print("Submission of results line '" + sendline + "' failed for reasons unknown - please try manual resubmission.")
			except URLError:
				debug_print("URL open error")

	write_list_file(sentfile, results_send, "a")	# EWM: Append entire results_send rather than just sent to avoid resubmitting
													# bad results (e.g. previously-submitted duplicates) every time the script executes.
	unlock_file(resultsfile)	# EWM: don't write anything to resultsfile, but still need to remove lock placed on it by read_list_file

parser = OptionParser()

parser.add_option("-d", "--debug", action="store_true", dest="debug", default=False, help="Display debugging info")

parser.add_option("-u", "--username", dest="username", help="Primenet user name")
parser.add_option("-p", "--password", dest="password", help="Primenet password")
parser.add_option("-w", "--workdir", dest="workdir", default=".", help="Working directory with worktodo.ini and results.txt, default current")

# -t is reserved for timeout, instead use -T for assignment-type preference:
parser.add_option("-T", "--worktype", dest="worktype", default="101", help="Worktype code, default %(default)s for DC, alternatively 100 (smallest available first-time LL) or 102 (world-record-sized first-time LL)")

parser.add_option("-n", "--num_cache", dest="num_cache", default="2", help="Number of assignments to cache, default 2")

parser.add_option("-t", "--timeout", dest="timeout", default="21600", help="Seconds to wait between network updates, default 21600 [6 hours]. Use 0 for a single update without looping.")

(options, args) = parser.parse_args()

progname = os.path.basename(sys.argv[0])
workdir = os.path.expanduser(options.workdir)
timeout = int(options.timeout)

workfile = os.path.join(workdir, "worktodo.ini")
resultsfile = os.path.join(workdir, "results.txt")

# A cumulative backup
sentfile = os.path.join(workdir, "results_sent.txt")

workpattern = r"(DoubleCheck|Test)=.*(,[0-9]+){3}"

# mersenne.org limit is about 4 KB; stay on the safe side
sendlimit = 3000

# adapted from http://stackoverflow.com/questions/923296/keeping-a-session-in-python-while-making-http-requests
primenet_cj = cookiejar.CookieJar()
primenet = build_opener(HTTPCookieProcessor(primenet_cj))

while True:
	# Log in to primenet
	try:
		login_data = {"user_login": options.username,
			"user_password": options.password,
		}

		# This makes a POST instead of GET
		try:
			data = login_data.encode('ascii')	# More python3 hackage to properly encode text data as bytes
		except AttributeError:
			data = urlencode(login_data)

		r = primenet.open(primenet_baseurl + b"default.php", data)
		if not options.username + "<br>logged in" in r.read():
			primenet_login = False
			debug_print("Login failed.")
		else:
			primenet_login = True
			while submit_work() == "locked":
				debug_print("Waiting for results file access...")
				sleep(2)
	except URLError:
		debug_print("Primenet URL open error")

	if primenet_login:
		while get_assignment() == "locked":
			debug_print("Waiting for worktodo.ini access...")
			sleep(2)
	if timeout <= 0:
		break
	sleep(timeout)
