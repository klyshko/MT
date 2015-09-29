#include "wrapper.h"
#include <malloc.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <execinfo.h>
#include <stdarg.h>

/* Obtain a backtrace and print it to stdout. */
/* Taken from http://www.gnu.org/software/libc/manual/html_node/Backtraces.html */
void print_trace(FILE *f) {
    void *array[70];
    size_t size;
    char **strings;
    size_t i;

    size = backtrace (array, 30);
    strings = backtrace_symbols (array, size);

    fprintf(f, "Stack trace (%zd frames):\n", size);

    for (i = 0; i < size; i++)
        fprintf(f, "%s\n", strings[i]);

    free(strings);
}

__attribute__((noinline,noreturn)) void die(const char *msg, ... ) {
    fflush(stdout); fflush(stderr);
    fprintf(stderr, "=====================================\n");
    fprintf(stderr, "Fatal error!\n");
    if (msg) {
        va_list args;
        va_start(args, msg);
        vfprintf(stderr, msg, args);
        va_end(args);
    }
    fprintf(stderr, "\n");
    print_trace(stderr);
    fflush(stdout); fflush(stderr);
    exit (-1);
}

inline void check_stream(FILE *stream) {
    check(!ferror(stream));
}

size_t safe_fread(void *ptr, size_t size, size_t nmemb, FILE *stream) {
    size_t ret;
    ret = fread(ptr, size, nmemb, stream);
    check_stream(stream);
    return ret;
}

char *safe_fgets(char *s, int size, FILE *stream) {
    char *ret;
    ret = fgets(s, size, stream);
    check_stream(stream);
    return ret;
}

FILE *safe_fopen(const char *path, const char *mode) {
    FILE *ret;
    ret = fopen(path, mode);
    if (!ret)
        DIE("Opening file '%s'", path);
    return ret;
}

