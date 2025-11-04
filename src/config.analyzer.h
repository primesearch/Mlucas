/* Dummy, minimal config.h for analyzer purposes. */
#include <string.h>
#ifndef _GNU_SOURCE
void *memrchr(const void *s, int c, size_t n);
#endif
