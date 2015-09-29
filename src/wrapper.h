#pragma once

#ifdef __cplusplus
# include <cstdio>
# define _DEF_NULL = NULL
#else
# include <stdio.h>
# define _DEF_NULL 
#endif

// Die and print some message (and may be stacktrace!)
void die(const char *msg _DEF_NULL, ...) __attribute__ ((__noreturn__));

#if __STDC_VERSION__ < 199901L
# if __GNUC__ >= 2
#  define __func__ __FUNCTION__
# else
#  define __func__ "<unknown>"
# endif
#endif

// Macro for priniting pretty post-mortem
#define DIE(format, ...) die("%s:%d [%s]: " #format, __FILE__, __LINE__, __func__, ##__VA_ARGS__)

#ifdef __cplusplus
// x != NULL or die
template <typename T> 
inline void check(T x, const char *msg = NULL) {
    if (!x) die(msg);
};
#endif

size_t safe_fread(void *ptr, size_t size, size_t nmemb, FILE *stream);
char *safe_fgets(char *s, int size, FILE *stream);
FILE *safe_fopen(const char *path, const char *mode);

