#include <stdio.h>
#include <stdlib.h>
FILE *fenkf;

int openenkffile_(int *orbNumber)
{
  char fname[100],s[6];
  sprintf(&s[0],"%6.6i",*orbNumber);
  strcpy(&fname[0],"US/enk2File");
  strcat(&fname[0],s);
  printf("%s \n",fname);
  fenkf=fopen(fname,"rb");
}

int main(void)
{
  int obNumber=2665;
  openenkffile_(&obNumber);
}
int closeenkffile_()
{
  fclose(fenkf);
}

void enkfwi5_(int *n5)
{
  fwrite(n5,sizeof(int),5,fenkf);
}

void enkfwi1_(int *n1)
{
  fwrite(n1,sizeof(int),1,fenkf);
}

void enkfwf_(float *f)
{
  fwrite(f,sizeof(float),1,fenkf);
}


void enkfri5_(int *n5)
{
  fread(n5,sizeof(int),5,fenkf);
}

void enkfri1_(int *n1)
{
  fread(n1,sizeof(int),1,fenkf);
}

void enkfrf_(float *f)
{
  fread(f,sizeof(float),1,fenkf);
  // printf("%g \n",*f);
}
