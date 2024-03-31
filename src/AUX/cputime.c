/*
   report statistics of CPU, memory usage
*/

#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

double cputime_()
{
  struct rusage *rusage;
  double t;

  rusage = (struct rusage *) malloc(sizeof(struct rusage));
  (void) getrusage(RUSAGE_SELF, rusage);

  if(rusage != (struct rusage *) 0) {
    t = rusage->ru_utime.tv_sec + rusage->ru_utime.tv_usec/1.e6;
  }
  else {
    t = 0.0;
  }
  free(rusage);

  return (t);
}
