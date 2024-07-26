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

#include "threadpool.h"
#include "util.h"	// This is to get (or not) <hwloc.h>

#ifdef MULTITHREAD	// Wrap contents of this file in flag (set via platform.h at compile time) ensuring no code built in unthreaded mode

// MacOS has its own versions of these, in /usr/include/X11/Xthreads.h:
static void * xmalloc(size_t len) {
	void *ptr = malloc(len);
	if (ptr == NULL) {
		printf("failed to allocate %u bytes\n", (uint32)len);
		exit(-1);
	}
	return ptr;
}

static void * xcalloc(size_t num, size_t len) {
	void *ptr = calloc(num, len);
	if (ptr == NULL) {
		printf("failed to calloc %u bytes\n", (uint32)(num * len));
		exit(-1);
	}
	return ptr;
}

	#define THREAD_POOL_DEBUG	0

	#if THREAD_POOL_DEBUG
		#define REPORT_ERROR(...) fprintf (stderr,"line %d - ",__LINE__); fprintf (stderr, __VA_ARGS__); fprintf (stderr,"\n")
	#else
		#define REPORT_ERROR(...)
	#endif /* THREAD_POOL_DEBUG ? */

	static void threadpool_queue_init(struct threadpool_queue *queue,
					int max_tasks)
	{
		queue->head = 0;
		queue->tail = 0;
		queue->num_tasks = 0;
		queue->max_tasks = max_tasks;
		queue->tasks = (void **)xcalloc(max_tasks, sizeof(void *));
	}

	static void threadpool_queue_free(struct threadpool_queue *queue)
	{
		free(queue->tasks);
	}

	static int threadpool_queue_enqueue(struct threadpool_queue *queue, void *data)
	{
		if (queue->num_tasks == queue->max_tasks) {
			REPORT_ERROR("The queue is full, unable to add data to it.");
			return -1;
		}

		if (queue->tasks[queue->tail] != NULL) {
			REPORT_ERROR("A problem was detected in the queue (expected NULL, but found a different value).");
			return -1;
		}

		queue->tasks[queue->tail] = data;

		queue->num_tasks++;
		queue->tail++;

		if (queue->tail == queue->max_tasks) {
			queue->tail = 0;
		}

		return 0;
	}

	/**
	 * This function removes and returns the head data element in the queue.
	 *
	 * @param queue The queue structure.
	 * @return On success a data element is returned, on failure NULL is returned.
	 */
	static void *threadpool_queue_dequeue(struct threadpool_queue *queue)
	{
		void *data;

		if (queue->num_tasks == 0) {
				REPORT_ERROR("Tried to dequeue from an empty queue.");
				return NULL;
		}

		data = queue->tasks[queue->head];

		queue->tasks[queue->head] = NULL;
		queue->num_tasks--;

		if (queue->num_tasks == 0) {
			queue->head = 0;
			queue->tail = 0;
		}
		else {
			queue->head++;
			if (queue->head == queue->max_tasks) {
				queue->head = 0;
			}
		}

		return data;
	}

	/**
	 * This function checks if a given queue is empty.
	 *
	 * @param queue The queue structure.
	 * @return 1 if the queue is empty, else 0.
	 */
	static int threadpool_queue_is_empty(struct threadpool_queue *queue)
	{
		if (queue->num_tasks == 0) {
			return 1;
		}

		return 0;
	}

	/**
	 * This function checks if a given queue is full.
	 *
	 * @param queue The queue structure.
	 * @return 1 if the queue is full, else 0.
	 */
	static int threadpool_queue_is_full(struct threadpool_queue *queue)
	{
		if (queue->num_tasks == queue->max_tasks) {
			return 1;
		}

		return 0;
	}

	/**
	 * This function queries for the size of the given queue argument.
	 *
	 * @param queue The queue structure.
	 * @return The size of the queue.
	 */
	static int threadpool_queue_getsize(struct threadpool_queue *queue)
	{
		return queue->num_tasks;
	}

	static void threadpool_task_clear(task_control_t *task)
	{
		memset(task, 0, sizeof(task_control_t));
	}

	/**
	 * This function obtains a queued task from the pool and returns it.
	 * If no such task is available the operation blocks.
	 *
	 * @param pool The thread pool structure.
	 * @return A task or NULL on error (or if thread pool should shut down).
	 */
	static task_control_t * threadpool_task_get_task(struct threadpool *pool)
	{
		task_control_t * task;

		if (pool->stop_flag) {
			/* The pool should shut down return NULL. */
			return NULL;
		}

		/* Obtain a task */
		if (pthread_mutex_lock(&(pool->mutex))) {
			perror("pthread_mutex_lock: ");
			return NULL;
		}

		while (threadpool_queue_is_empty(&(pool->tasks_queue)) && !pool->stop_flag) {
			/* Block until a new task arrives. */
			if (pthread_cond_wait(&(pool->new_tasks_cond),&(pool->mutex))) {
				perror("pthread_cond_wait: ");
				if (pthread_mutex_unlock(&(pool->mutex))) {
					perror("pthread_mutex_unlock: ");
				}

				return NULL;
			}
		}

		if (pool->stop_flag) {
			/* The pool should shut down return NULL. */
			if (pthread_mutex_unlock(&(pool->mutex))) {
				perror("pthread_mutex_unlock: ");
			}
			return NULL;
		}

		if ((task = (task_control_t *)threadpool_queue_dequeue(&(pool->tasks_queue))) == NULL) {
			/* Since task is NULL returning task will return NULL as required. */
			REPORT_ERROR("Failed to obtain a task from the jobs queue.");
		}

		if (pthread_mutex_unlock(&(pool->mutex))) {
			perror("pthread_mutex_unlock: ");
			return NULL;
		}

		return task;
	}

	/**
	 * This is the routine the worker threads do during their life.
	 *
	 * @param data Contains a pointer to the startup data
	 * @return NULL.
	 */
	#if defined(__GNUC__) && (__GNUC__ > 4 || __GNUC__ == 4 && __GNUC_MINOR__>1) && defined(CPU_IS_X86)

	/* gcc on win32 needs to force 16-byte stack alignment on
	   thread entry, as this exceeds what windows may provide; see

	   http://sourceware.org/ml/pthreads-win32/2008/msg00053.html
	*/
	__attribute__((force_align_arg_pointer))
	#endif
	static void *worker_thr_routine(void *data)
	{
	#if INCLUDE_HWLOC
		char cbuf[STR_MAX_LEN*2];
		char str[80];
	#endif
		struct thread_init *init = (struct thread_init *)data;
		int my_id = init->thread_num;
		struct threadpool *pool = init->pool;
		thread_control_t *t = &init->control;
		task_control_t *task;

		// Set CPU affinity masks of the thread:
	#ifdef __OpenBSD__

	  #warning OpenBSD has no user-affinity-setting support ... affinity-setting is up to the OS.

	#elif defined(__FreeBSD__)

	  #warning FreeBSD affinity code needs debug & test!

		int i,errcode;
		cpuset_t cpu_set;

		CPU_ZERO (&cpu_set);
		i = my_id % pool->num_of_cores;	// get cpu mask using sequential thread ID modulo #available cores
		i = mi64_ith_set_bit(CORE_SET, i+1, MAX_CORES>>6);	// Remember, [i]th-bit index in arglist is *unit* offset, i.e. must be in [1,MAX_CORES]
		if(i < 0) {
			fprintf(stderr,"Affinity CORE_SET does not have a [%u]th set bit!",my_id % pool->num_of_cores);
			ASSERT(0, "Aborting.");
		}
		CPU_SET(i, &cpu_set);
		errcode = cpuset_setaffinity(CPU_LEVEL_WHICH, CPU_WHICH_TID, -1, sizeof (cpu_set), &cpu_set);
		if (errcode) {
			perror("cpuset_setaffinity");
		}

	#elif defined(OS_TYPE_LINUX) && !defined(__MINGW32__)

	  #if 0
	  /*
		// This is the affinity API tied to pthread library ... interestingly, it's less portable than the
		// Linux system-centric one below; e.g. GCC gives "error: unknown type name ‘cpuset_t’; did you mean ‘cpu_set_t’?" here:
		int i,errcode;
		cpuset_t *cset;
		pthread_t pth;

		cset = cpuset_create();
		if (cset == NULL) {
			err(EXIT_FAILURE, "cpuset_create");
		}
		i = my_id % pool->num_of_cores;	// get cpu mask using sequential thread ID modulo #available cores
		cpuset_set((cpuid_t)i, cset);

		pth = pthread_self();
		errcode = pthread_setaffinity_np(pth, cpuset_size(cset), cset);
		if (errcode) {
			perror("pthread_setaffinity_np");
		}
		cpuset_destroy(cset);
	  */
	  #else

		cpu_set_t cpu_set;
		int i,errcode;
		pid_t thread_id = syscall (__NR_gettid);
	  #if THREAD_POOL_DEBUG
		printf("executing worker thread id %u, syscall_id = %u\n", my_id, thread_id);
	  #endif

		i = my_id % pool->num_of_cores;
		i = mi64_ith_set_bit(CORE_SET, i+1, MAX_CORES>>6);	// Remember, [i]th-bit index in arglist is *unit* offset, i.e. must be in [1,MAX_CORES]
		if(i < 0) {
			fprintf(stderr,"Affinity CORE_SET does not have a [%u]th set bit!",my_id % pool->num_of_cores);
			ASSERT(0, "Aborting.");
		}

	 #if INCLUDE_HWLOC

	  if(HWLOC_AFFINITY) {	// Global, declared in Mdata.h, defined in Mlucas.c, set in util.c::host_init()
		hwloc_bitmap_t cpuset = hwloc_bitmap_alloc();
		hwloc_obj_t obj = hwloc_get_obj_by_type(hw_topology, HWLOC_OBJ_PU, i);
		if (obj) {
			hwloc_bitmap_or(cpuset, cpuset, obj->cpuset);
		} else {
			snprintf(cbuf,STR_MAX_LEN*2,"[hwloc] Error: HWLOC_OBJ_PU[%u] not found.\n",i);
			fprintf(stderr,"%s",cbuf);
		}
		// Set affinity to specified logical CPUs:
		if (hwloc_set_cpubind(hw_topology, cpuset, HWLOC_CPUBIND_THREAD)) {
			int error = errno;
			hwloc_bitmap_snprintf (str, sizeof (str), cpuset);
			snprintf(cbuf,STR_MAX_LEN*2,"[hwloc] Warning: Unable to set affinity to cpuset %s: %s; leaving up to OS to manage thread/core binding.\n",str,strerror(error));
			fprintf(stderr,"%s",cbuf);
	  #if THREAD_POOL_DEBUG
		} else {
			printf("[hwloc] tid = %d: HWLOC_OBJ_PU[%u], lidx %u, pidx %u: setaffinity[%d] to cpuset %s\n",my_id,i,obj->logical_index,obj->os_index,str);
	  #endif
		}
		hwloc_bitmap_free(cpuset);
	  }	// HWLOC_AFFINITY = True?

	 #else	// INCLUDE_HWLOC = False:

		// get cpu mask using sequential thread ID modulo #available cores in runtime-specified affinity set:
		CPU_ZERO (&cpu_set);
		CPU_SET(i, &cpu_set);
		errcode = sched_setaffinity(thread_id, sizeof(cpu_set), &cpu_set);
	  #if THREAD_POOL_DEBUG
		printf("syscall_id = %u, tid = %d, setaffinity[%d] = %d, ISSET[%d] = %d\n", thread_id,my_id,i,errcode,i,CPU_ISSET(i, &cpu_set));
	  #endif
		if (errcode) {
			perror("sched_setaffinity");
			fprintf(stderr,"INFO: Your run should be OK, but leaving up to OS to manage thread/core binding.\n");
		}
	  #endif

	 #endif	// INCLUDE_HWLOC?

	#elif defined(OS_TYPE_MACOSX)

		thread_t thr = mach_thread_self();
		thread_extended_policy_data_t epolicy;
		epolicy.timeshare = FALSE;
		kern_return_t ret = thread_policy_set(
			thr, THREAD_EXTENDED_POLICY,
			(thread_policy_t) &epolicy, THREAD_EXTENDED_POLICY_COUNT);
		if (ret != KERN_SUCCESS) {
			printf("thread_policy_set returned %d", ret);
			exit(-1);
		}

		thread_affinity_policy_data_t apolicy;
		int i = my_id % pool->num_of_cores;	// get cpu mask using sequential thread ID modulo #available cores
		i = mi64_ith_set_bit(CORE_SET, i+1, MAX_CORES>>6);	// Remember, [i]th-bit index in arglist is *unit* offset, i.e. must be in [1,MAX_CORES]
		if(i < 0) {
			fprintf(stderr,"Affinity CORE_SET does not have a [%u]th set bit!",my_id % pool->num_of_cores);
			ASSERT(0, "Aborting.");
		}
		apolicy.affinity_tag = i; // set affinity tag
	  #if THREAD_POOL_DEBUG
		printf("Setting CPU = %d affinity of worker thread id %u, mach_id = %u\n", i, my_id, thr);
	  #endif

		ret = thread_policy_set(
			thr, THREAD_EXTENDED_POLICY,
			(thread_policy_t) &apolicy, THREAD_EXTENDED_POLICY_COUNT);
		if (ret != KERN_SUCCESS) {
			printf("thread_policy_set returned %d", ret);
			exit(-1);
		}

	#else

		printf("executing worker thread id %u, #cores = %u\n", my_id, pool->num_of_cores);

	#endif

		/* initialize thread-local state */

		if (t->init != NULL) {
			t->init(t->data, my_id);
		}

		while (1) {
			task = threadpool_task_get_task(pool);
			if (task == NULL) {
				if (pool->stop_flag) {
					/* Worker thr needs to exit (thread pool was shutdown). */
					break;
				}
				else {
					/* An error has occurred. */
					REPORT_ERROR("Warning an error has occurred when trying to obtain a worker task.");
					REPORT_ERROR("The worker thread has exited.");
					break;
				}
			}

			/* Execute task */

			if (task->init != NULL)
				task->init(task->data, my_id);

			if (task->run != NULL)
				task->run(task->data, my_id);

			if (task->shutdown != NULL)
				task->shutdown(task->data, my_id);

			/* Release the task by returning it to the free_task_queue. */
			threadpool_task_clear(task);
			if (pthread_mutex_lock(&(pool->free_tasks_mutex))) {
				perror("pthread_mutex_lock: ");
				REPORT_ERROR("The worker thread has exited.");
				break;
			}

			if (threadpool_queue_enqueue(&(pool->free_tasks_queue),task)) {
				REPORT_ERROR("Failed to enqueue a task to free tasks queue.");
				if (pthread_mutex_unlock(&(pool->free_tasks_mutex))) {
					perror("pthread_mutex_unlock: ");
				}

				REPORT_ERROR("The worker thread has exited.");
				break;
			}

			if (threadpool_queue_getsize(&(pool->free_tasks_queue)) == 1) {
				/* Notify all waiting threads that new tasks can added. */
				if (pthread_cond_broadcast(&(pool->free_tasks_cond))) {
					perror("pthread_cond_broadcast: ");
					if (pthread_mutex_unlock(&(pool->free_tasks_mutex))) {
						perror("pthread_mutex_unlock: ");
					}

					break;
				}
			}

			if (threadpool_queue_is_full(&(pool->free_tasks_queue)) == 1) {
				/* Notify any waiting threads that threadpool is not busy */

				if (pthread_cond_broadcast(&(pool->tasks_done_cond))) {
					perror("pthread_cond_broadcast: ");
					if (pthread_mutex_unlock(&(pool->free_tasks_mutex))) {
						perror("pthread_mutex_unlock: ");
					}

					break;
				}
			}

			if (pthread_mutex_unlock(&(pool->free_tasks_mutex))) {
				perror("pthread_mutex_unlock: ");
				REPORT_ERROR("The worker thread has exited.");
				break;
			}
		}

		/* tear down thread-local state */

		if (t->shutdown != NULL) {
			t->shutdown(t->data, my_id);
		}

		return NULL;
	}

	/**
	 * This function does the following steps:
	 * 1. It raises a flag that notifies the worker threads to stop working.
	 * 2. It waits until all worker threads are done with their execution.
	 * 3. It frees all the allocated memory of the threadpool struct.
	 *
	 * @param ptr The pool to stop its worker threads.

	 * @return 0.
	 */
	void threadpool_free(struct threadpool *pool)
	{
		int i;

		pool->stop_flag = 1;

		/* Wakeup all worker threads (broadcast operation). */
		if (pthread_mutex_lock(&(pool->mutex))) {
			perror("pthread_mutex_lock: ");
			REPORT_ERROR("Warning: Memory was not released.");
			REPORT_ERROR("Warning: Some of the worker threads may have failed to exit.");
			return;
		}

		if (pthread_cond_broadcast(&(pool->new_tasks_cond))) {
			perror("pthread_cond_broadcast: ");
			REPORT_ERROR("Warning: Memory was not released.");
			REPORT_ERROR("Warning: Some of the worker threads may have failed to exit.");
			return;
		}

		if (pthread_mutex_unlock(&(pool->mutex))) {
			perror("pthread_mutex_unlock: ");
			REPORT_ERROR("Warning: Memory was not released.");
			REPORT_ERROR("Warning: Some of the worker threads may have failed to exit.");
			return;
		}

		/* Wait until all worker threads are done. */
		for (i = 0; i < pool->num_of_threads; i++) {
			if (pthread_join(pool->thr_arr[i],NULL)) {
				perror("pthread_join: ");
			}
		}

		/* shut down any tasks that are still waiting */
		while (threadpool_queue_getsize(&(pool->tasks_queue))) {

			task_control_t *task = (task_control_t *)
					threadpool_queue_dequeue(&(pool->tasks_queue));

			if (task != NULL && task->shutdown != NULL) {
				task->shutdown(task->data, 0);
			}
		}

		/* Free all allocated memory. */
		threadpool_queue_free(&(pool->tasks_queue));
		threadpool_queue_free(&(pool->free_tasks_queue));
		free(pool->tasks);
		free(pool->thr_arr);
		free(pool->thr_init);
		free(pool);
	}

	struct threadpool* threadpool_init(
				int num_threads,
				int num_cores,
				int queue_size,
				thread_control_t *t)
	{
		int i;
		struct threadpool *pool = (struct threadpool *)xcalloc(1,
						sizeof(struct threadpool));

		/* Init the mutex and cond vars. */
		if (pthread_mutex_init(&(pool->free_tasks_mutex),NULL)) {
			perror("pool->free_tasks_mutex init: ");
			free(pool);
			return NULL;
		}
		if (pthread_mutex_init(&(pool->mutex),NULL)) {
			perror("pool->mutex init: ");
			free(pool);
			return NULL;
		}
		if (pthread_cond_init(&(pool->free_tasks_cond),NULL)) {
			perror("pool->free_tasks_cond init: ");
			free(pool);
			return NULL;
		}
		if (pthread_cond_init(&(pool->tasks_done_cond),NULL)) {
			perror("pool->tasks_done_cond init: ");
			free(pool);
			return NULL;
		}
		if (pthread_cond_init(&(pool->new_tasks_cond),NULL)) {
			perror("pool->new_tasks_cond init: ");
			free(pool);
			return NULL;
		}

		/* Init the queues. */
		threadpool_queue_init(&(pool->tasks_queue), queue_size);
		threadpool_queue_init(&(pool->free_tasks_queue), queue_size);
		pool->tasks = (task_control_t *)xmalloc(queue_size *
						sizeof(task_control_t));

		/* Add all the free tasks to the free tasks queue. */
		for (i = 0; i < queue_size; i++) {
			threadpool_task_clear((pool->tasks) + i);
			if (threadpool_queue_enqueue(&(pool->free_tasks_queue),(pool->tasks) + i)) {
				REPORT_ERROR("Failed to a task to the free tasks queue during initialization.");
				return NULL;
			}
		}

		/* Create the thr_arr. */
		if ((pool->thr_arr = malloc(sizeof(pthread_t) * num_threads)) == NULL) {
			perror("malloc: ");
			free(pool);
			return NULL;
		}

		if ((pool->thr_init = malloc(sizeof(struct thread_init) * num_threads)) == NULL) {
			perror("malloc: ");
			free(pool);
			return NULL;
		}

		/* Start the worker threads. */
		for (i = 0; i < num_threads; i++)
		{
			pool->num_of_threads = num_threads;
			pool->num_of_cores = num_cores;

			pool->thr_init[i].thread_num = i;
			pool->thr_init[i].pool = pool;
			pool->thr_init[i].control = *t;

			if (pthread_create(&(pool->thr_arr[i]),NULL,
					worker_thr_routine,
					&(pool->thr_init[i]))) {
				perror("pthread_create:");

				threadpool_free(pool);
				return NULL;
			}
		}

		// Set CPU affinity masks of the threads:
	#if 0//def OS_TYPE_LINUX

		cpu_set_t cpu_set;
		pid_t thread_id;
		int errcode;

		CPU_ZERO (&cpu_set);
		for (i = 0; i < num_threads; i++)
		{
			CPU_SET(i, &cpu_set);
			thread_id = (pid_t) syscall (__NR_gettid);
			errcode = sched_setaffinity (thread_id, sizeof(cpu_set), &cpu_set);
		  #if THREAD_POOL_DEBUG
			printf("syscall_id = %u, setaffinity[%d] = %d, ISSET[%d] = %d\n", thread_id,i,errcode,i,CPU_ISSET(i, &cpu_set));
		  #endif
		}

	#elif 0//OS_TYPE_MACOSX

	// above code gives these diagnostcs for 2-threaded:
	syscall_id = 32252, setaffinity[0] = 0, ISSET[0] = 1
	syscall_id = 32252, setaffinity[1] = 0, ISSET[1] = 1
	Fermat_mod_square: Init threadpool of 2 threads
	Thread 0, self_id = 2683893504, sys_id = 32267
	Thread 1, self_id = 2692286208, sys_id = 32268
	Thread 2, self_id = 2683893504, sys_id = 32267
	Thread 3, self_id = 2692286208, sys_id = 32268

		#error To-Do: Add thread-affinity support for MacOS!

	#endif

		return pool;
	}

	int threadpool_add_task(struct threadpool *pool, task_control_t *new_task, int blocking)
	{
		task_control_t *task;

		if (pool == NULL) {
			REPORT_ERROR("The threadpool received as argument is NULL.");
			return -1;
		}

		if (pthread_mutex_lock(&(pool->free_tasks_mutex))) {
			perror("pthread_mutex_lock: ");
			return -1;
		}

		/* Check if the free task queue is empty. */
		while (threadpool_queue_is_empty(&(pool->free_tasks_queue))) {
			if (!blocking) {
				/* Return immediately if the command is non blocking. */
				if (pthread_mutex_unlock(&(pool->free_tasks_mutex))) {
					perror("pthread_mutex_unlock: ");
					return -1;
				}

				return -2;
			}

			/* blocking is set to 1, wait until free_tasks queue has a task to obtain. */
			if (pthread_cond_wait(&(pool->free_tasks_cond),&(pool->free_tasks_mutex))) {
				perror("pthread_cond_wait: ");
				if (pthread_mutex_unlock(&(pool->free_tasks_mutex))) {
					perror("pthread_mutex_unlock: ");
				}

				return -1;
			}
		}

		/* Obtain an empty task. */
		if ((task = (task_control_t *)threadpool_queue_dequeue(&(pool->free_tasks_queue))) == NULL) {
			REPORT_ERROR("Failed to obtain an empty task from the free tasks queue.");
			if (pthread_mutex_unlock(&(pool->free_tasks_mutex))) {
				perror("pthread_mutex_unlock: ");
			}

			return -1;
		}

		if (pthread_mutex_unlock(&(pool->free_tasks_mutex))) {
			perror("pthread_mutex_unlock: ");
			return -1;
		}

		*task = *new_task;

		/* Add the task, to the tasks queue. */
		if (pthread_mutex_lock(&(pool->mutex))) {
			perror("pthread_mutex_lock: ");
			return -1;
		}

		if (threadpool_queue_enqueue(&(pool->tasks_queue),task)) {
			REPORT_ERROR("Failed to add a new task to the tasks queue.");
			if (pthread_mutex_unlock(&(pool->mutex))) {
				perror("pthread_mutex_unlock: ");
			}
			return -1;
		}

		if (threadpool_queue_getsize(&(pool->tasks_queue)) == 1) {
			/* Notify all worker threads that there are new jobs. */
			if (pthread_cond_broadcast(&(pool->new_tasks_cond))) {
				perror("pthread_cond_broadcast: ");
				if (pthread_mutex_unlock(&(pool->mutex))) {
					perror("pthread_mutex_unlock: ");
				}

				return -1;
			}
		}

		if (pthread_mutex_unlock(&(pool->mutex))) {
			perror("pthread_mutex_unlock: ");
			return -1;
		}

		return 0;
	}

	int threadpool_drain(struct threadpool *pool, int blocking)
	{
		if (pthread_mutex_lock(&(pool->free_tasks_mutex))) {
			perror("pthread_mutex_lock: ");
			return -1;
		}

		/* Check if the free task queue is full. */
		while (!threadpool_queue_is_full(&(pool->free_tasks_queue))) {
			if (!blocking) {
				/* Return immediately if the command is non blocking. */
				if (pthread_mutex_unlock(&(pool->free_tasks_mutex))) {
					perror("pthread_mutex_unlock: ");
					return -1;
				}

				return 1;
			}

			/* blocking is set to 1, wait until free_tasks queue is full */
			if (pthread_cond_wait(&(pool->tasks_done_cond),&(pool->free_tasks_mutex))) {
				perror("pthread_cond_wait: ");
				if (pthread_mutex_unlock(&(pool->free_tasks_mutex))) {
					perror("pthread_mutex_unlock: ");
				}

				return -1;
			}
		}

		if (pthread_mutex_unlock(&(pool->free_tasks_mutex))) {
			perror("pthread_mutex_unlock: ");
			return -1;
		}

		return 0;
	}

#endif	// ifdef MULTITHREAD ?
