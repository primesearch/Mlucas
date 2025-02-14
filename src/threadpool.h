// EWM: This threadpool file has "more history" than the other Mlucas sources,
// thus I include public-license boilerplate from all 3 authors in the chain below.
// [1] First my standard GPL header covering the code including my customizations:

/*******************************************************************************
*                                                                              *
*   (C) 1997-2021 by Ernst W. Mayer.                                           *
*                                                                              *
*  This program is free software; you can redistribute it and/or modify it     *
*  under the terms of the GNU General Public License as published by the       *
*  Free Software Foundation; either version 2 of the License, or (at your      *
*  option) any later version.                                                  *
*                                                                              *
*  This program is distributed in the hope that it will be useful, but WITHOUT *
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
*  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for   *
*  more details.                                                               *
*                                                                              *
*  You should have received a copy of the GNU General Public License along     *
*  with this program; see the file GPL.txt.  If not, you may view one at       *
*  http://www.fsf.org/licenses/licenses.html, or obtain one by writing to the  *
*  Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA     *
*  02111-1307, USA.                                                            *
*                                                                              *
*******************************************************************************/

// [2] The version is started with was sent by , whose latest analogs of
// same are available at 
// http://sourceforge.net/p/msieve/code/HEAD/tree/trunk/include/thread.h [header]
// http://sourceforge.net/p/msieve/code/HEAD/tree/trunk/common/thread.c [C source].

/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id$
--------------------------------------------------------------------*/

// [3] Jason informs me he first started with Tomer Heber's code at
// http://sourceforge.net/projects/cthreadpool/ , which has the following
// license info (Note that BSD and GPL licenses are more-or-less compatible):

/*--------------------------------------------------------------------
The cthreadpool project is free and open source (BSD License).
If you are not familiar with the thread pool pattern please refer to:
http://en.wikipedia.org/wiki/Thread_pool_pattern

Instructions:
1. In order to use the threadpool add threadpool.c and threadpool.h to
your project and compile it (compiling will create an object(o) file).
2. Check the file threadpool.h for the API.
3. Examples are available in the file example.c.

For questions, suggestions, bug reports or just comments please contact
me at: heber.tomer@gmail.com
--------------------------------------------------------------------*/

#ifndef _THREAD_H_
#define _THREAD_H_

#include "masterdefs.h"
#include "types.h"
#include "mi64.h"	// Sep 2016: Needed for enhanced affinity-setting functionality

#include <pthread.h>

#ifdef __cplusplus
extern "C"
{
#endif

/* mutexes ---------------------------------------------------------*/

#ifdef OS_TYPE_WINDOWS
typedef HANDLE mutex_t;
#else
typedef pthread_mutex_t mutex_t;
#endif
/*
static void mutex_init(mutex_t *m)
{
#ifdef OS_TYPE_WINDOWS
	*m = CreateMutex(NULL, FALSE, NULL);
#else
	pthread_mutex_init(m, NULL);
#endif
}

static void mutex_free(mutex_t *m)
{
#ifdef OS_TYPE_WINDOWS
	CloseHandle(*m);
#else
	pthread_mutex_destroy(m);
#endif
}

static void mutex_lock(mutex_t *m)
{
#ifdef OS_TYPE_WINDOWS
	WaitForSingleObject(*m, INFINITE);
#else
	pthread_mutex_lock(m);
#endif
}

static void mutex_unlock(mutex_t *m)
{
#ifdef OS_TYPE_WINDOWS
	ReleaseMutex(*m);
#else
	pthread_mutex_unlock(m);
#endif
}
*/
/* a thread pool --------------------------------------------------*/

typedef void (*init_func)(void *data, int thread_num);
typedef void (*run_func)(void *data, int thread_num);
typedef void (*shutdown_func)(void *data, int thread_num);

typedef struct {
	init_func init;
	shutdown_func shutdown;
	void *data;
} thread_control_t;

typedef struct {
	init_func init;
	run_func run;
	shutdown_func shutdown;
	void *data;
} task_control_t;

struct threadpool_queue
{
	unsigned int head;
	unsigned int tail;
	unsigned int num_tasks;
	unsigned int max_tasks;
	void **tasks;
};

struct thread_init 
{
	int thread_num;
	struct threadpool *pool;
	thread_control_t control;
};

struct threadpool
{
	struct threadpool_queue tasks_queue;
	struct threadpool_queue free_tasks_queue;

	task_control_t *tasks;

	struct thread_init *thr_init;
	pthread_t *thr_arr;

	unsigned short num_of_threads;
	unsigned short num_of_cores;
	volatile unsigned short stop_flag;

	pthread_mutex_t free_tasks_mutex;
	pthread_cond_t free_tasks_cond;
	pthread_cond_t tasks_done_cond;

	pthread_mutex_t mutex;
	pthread_cond_t new_tasks_cond;
};

struct threadpool* threadpool_init(
			int num_threads, 
			int num_cores, 
			int queue_size, 
			thread_control_t *t);

int threadpool_add_task(struct threadpool *pool, 
			task_control_t *t, 
			int blocking);

void threadpool_free(struct threadpool *pool);

/* returns zero if no pending tasks */
int threadpool_drain(struct threadpool *pool,
			int blocking);

#ifdef __cplusplus
}
#endif

#endif /* !_THREAD_H_ */

