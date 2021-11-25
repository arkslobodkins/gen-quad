#ifndef GET_TIME_H
#define GET_TIME_H

#ifdef __cplusplus
extern "C" {
#endif

// from https://github.com/pjreddie/darknet/blob/master/src/utils.h
#define TIME(a) \
do {                                                                    \
   double start = get_cur_time();                                       \
   a;                                                                   \
   printf("%s took: %f seconds\n", #a, get_cur_time - start);           \
   } while (0)


double get_cur_time();

#ifdef __cplusplus
}
#endif


#endif
