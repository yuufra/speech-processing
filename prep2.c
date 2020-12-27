/*
usage
gcc prep2.c (ハミング窓をかけないとき)
gcc prep2.c -DWINDOW -lm (ハミング窓をかけるとき)
./a.out
gnuplot
plot "dat_prep2.txt" with lines linetype 1
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){
    FILE* fp = fopen("speech_sample/A_a.wav","rb");
    int n = 1000000;
    short* data = (short*) malloc(sizeof(short)*n);
    int m = fread(data,sizeof(short),n,fp);
    fclose(fp);

    // FILE* gp = popen("gnuplot --persist","w");
    // fprintf(gp,"plot '-' with lines linetype 1\n");
    FILE* gp = fopen("dat_prep2.txt","w");
#ifdef WINDOW
    for (int i=0;i<n;i++)
        data[i] *= 0.54-0.46*cos(2*i*M_PI/(n-1));
#endif
    for (int i=22;i<m;i++) //wavでは先頭22バイトを読み飛ばす
        fprintf(gp,"%d\t%d\n",i-22,data[i]);
    fprintf(gp,"e\n");
    // pclose(gp);
    fclose(gp);

    return 0;
}