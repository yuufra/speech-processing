#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv)
{
  /* check the format of input */
  if (argc != 4)
  {
    fprintf(stderr, "Usage: %s DATfile skip[sample] frame_length[sample]\n", argv[0]);
    exit(1);
  }
  FILE *fpDAT;
  int nskip;
  int framelen;
  int i,j;

  /* check the validity of input */
  if ((fpDAT = fopen(argv[1], "r")) == NULL)
    exit(1);
  if ((nskip = atoi(argv[2])) < 0)
    exit(1);
  if ((framelen = atoi(argv[3])) < 0)
    exit(1);

  fprintf(stderr, "# DATfile = %s\n", argv[1]);
  fprintf(stderr, "# %d samples are skipped.\n", nskip);
  fprintf(stderr, "# 1 frame contains %d sampels.\n", framelen);

  /* memory allocation & initilization */
  /* calloc() puts zero-values for assigned memories. */
  short *sdata = (short *)calloc(framelen, sizeof(short));
  double *cor = (double *)calloc(framelen, sizeof(double));
  if (sdata == NULL)
    exit(1);

  fseek(fpDAT, nskip * sizeof(short), SEEK_SET);
  fread(sdata, sizeof(short), framelen, fpDAT);
  fclose(fpDAT);

  for (i = 0; i < framelen; i++){
    for (j = 0; j < framelen-1-i; j++){
      cor[i] += (double)sdata[j] * sdata[j+i];
    }
    if (i != 0) 
      cor[i] /= cor[0];
  }
  cor[0] = 1.0;

  for (i = 0; i < framelen; i++){
      printf("%lf %lf\n", (double)i/16000, cor[i]);
  }

  int argmax = 1;
  for (i = 2; i < framelen; i++){
    if (cor[i] > cor[argmax])
      argmax = i;
  }

  return 0;
}
