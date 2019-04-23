#include <iostream>
#include <sstream>
#include <stdio.h>
using namespace std;
 
void main()
{
	const int Nquad = 19;
    FILE *fp;
    stringstream filename;
    filename << Nquad << ".txt";
    fp = fopen(filename.str().c_str(), "r");     /*  読み込みモードでファイルをオープン  */

    if(fp == NULL) {
        printf("ファイルを開くことが出来ませんでした．\n");
        return;
    }
 
 	double w, x, y, z;
    for(int i=0; i<Nquad; i++){
        fscanf(fp, "%lf %lf %lf %lf", &w, &x, &y, &z );
        printf("%d %e %e %e %e \n", i, w, x, y, z);
    }
 
    fclose(fp);
}