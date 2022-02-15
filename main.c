#include <stdio.h>

#define _In_
#define _Out_

#define GAMMA 0.91
#define TIME1 2245
#define TIME2 342
#define SCOPE_SIZE 2

#define DIV_SIZE 10

static double sample_get_t_avg(_In_ const int *sample, _In_ int sample_len) {
  double t_avg = 0;
  for (int i = 0; i < sample_len; ++i) t_avg += sample[i];
  return t_avg / sample_len;
}

static void sample_sort(_In_ int *sample, _In_ int sample_size) {
  for (int i = 0; i < sample_size; ++i) {
    for (int j = i + 1; j < sample_size; ++j) {
      if (sample[i] > sample[j]) {
        sample[i] ^= sample[j];
        sample[j] ^= sample[i];
        sample[i] ^= sample[j];
      }
    }
  }
}

static void sample_print(_In_ const int *sample, _In_ int sample_size) {
  printf("sample:\n");
  for (int i = 0; i < sample_size; ++i) printf("%d ", sample[i]);
  printf("\n");
}

static double *sample_create_intervals(_In_ int max_val,
                                       _In_ double interval_size,
                                       _Out_ double *intervals) {
  int j = 1;
  for (double i = 0; i < max_val; i += interval_size, j++) intervals[j - 1] = i;
  return intervals;
}

static int sample_f_interval_val(_In_ const int *sample,
                                 _In_ int sample_len,
                                 _In_ double from,
                                 _In_ double to) {
  int count = 0;
  for (int i = 0; i < sample_len; ++i) {
    if (sample[i] > from && sample[i] <= to) {
      count++;
    }
  }
  return count;
}

static void sample_f_t(_In_ int *sample,
                       _In_ int sample_len,
                       _In_ double *intervals,
                       _In_ double interval_size,
                       _Out_ double *f) {
  for (int i = 0; i < DIV_SIZE; i++) {
    f[i] = sample_f_interval_val(sample, sample_len, intervals[i], intervals[i] + interval_size);
    f[i] = f[i] / (sample_len * interval_size);
  }
}

static void sample_p_t(_In_ const double *f, _In_ double interval_size, _Out_ double *p) {
  p[0] = 1.0;
  double pi = 0;
  for (int i = 0; i < DIV_SIZE; i++) {
    pi += f[i] * interval_size;
    p[i + 1] = 1.0 - pi;
  }
}

inline static double sample_t_gamma(_In_ const double *p,
                                    _In_ double interval_size,
                                    _In_ double gamma,
                                    _In_ int from,
                                    _In_ int to) {
  double d = (p[to] - gamma) / (p[to] - p[from]);
  return interval_size - interval_size * d;
}

static int sample_index(_In_ const double *intervals, _In_ double interval_size, _In_ double time) {
  int index = 0;
  for (int i = 0; i < DIV_SIZE; i++) {
    if (time >= intervals[i] && time <= intervals[i] + interval_size) {
      index = i;
      break;
    }
  }
  return index;
}

static double sample_p_time(_In_ const double *f,
                            _In_ const double *intervals,
                            _In_ double interval_size,
                            _In_ double time) {
  int index = sample_index(intervals, interval_size, time);
  double pTime = 0;
  for (int i = 0; i < index; i++) pTime += f[i] * interval_size;
  pTime += f[index] * (time - intervals[index]);
  return 1.0 - pTime;
}

static void sample_d(_In_ const double *arr, int count, _Out_ int *d) {
  for (int i = 0; i < count - 1; i++) {
    double from = arr[i];
    double to = arr[i + 1];
    if (GAMMA <= from && GAMMA > to) {
      d[0] = i;
      d[1] = i + 1;
    }
  }
}

int main(void) {
  int sample[] = {
          577, 566, 251, 868, 2525, 26, 1040, 580, 57,
          376, 906, 304, 1407, 134, 62, 230, 531, 653,
          2143, 515, 701, 88, 762, 337, 278, 87, 85,
          10, 484, 23, 454, 61, 508, 15, 469, 309, 434,
          2323, 219, 120, 66, 212, 220, 498, 522, 267,
          820, 648, 18, 417, 519, 34, 662, 756, 461,
          221, 291, 224, 234, 365, 416, 114, 448, 395,
          71, 4, 44, 462, 537, 226, 125, 157, 18, 216,
          844, 438, 730, 274, 395, 40, 481, 118, 14,
          1297, 102, 249, 34, 758, 611, 159, 1115,
          353, 117, 15, 158, 259, 489, 38, 1416, 2344};

  int sample_len = sizeof(sample) / sizeof(sample[0]);
  double t_avg = sample_get_t_avg(sample, sample_len);
  printf("t_avg: %f\n", t_avg);

  sample_sort(sample, sample_len);
  sample_print(sample, sample_len);

  int max_val = sample[sample_len - 1];
  double interval_size = (double) max_val / DIV_SIZE;
  printf("%d", max_val);

  double intervals[DIV_SIZE];
  sample_create_intervals(max_val, interval_size, intervals);

  double f[DIV_SIZE];
  sample_f_t(sample, sample_len, intervals, interval_size, f);

  double p[DIV_SIZE + 1];
  sample_p_t(f, interval_size, p);

  int scope[SCOPE_SIZE];
  sample_d(p, DIV_SIZE + 1, scope);

  double t = sample_t_gamma(p, interval_size, GAMMA, scope[0], scope[1]);
  double p_time = sample_p_time(f, intervals, interval_size, TIME1);
  double p_time_int = sample_p_time(f, intervals, interval_size, TIME2);
  int index2 = sample_index(intervals, interval_size, TIME2);
  printf("\nT: %f", t);
  printf("\nP(%d): %f", TIME1, p_time);
  printf("\nP(%d): %f", TIME2, p_time_int);
  printf("\nGamma(%d): %f", TIME2, f[index2] / p_time_int);
  return 0;
}
