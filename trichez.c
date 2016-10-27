/*
Código feito por Matheus Henrique Trichez,
Matrícula - 1211100047,
UFFS - Universidade Federal da Fronteira Sul,
Curso - Ciência da Computação,
CCR - Álgebra Linear.
----------------------------------
Método de ortogonalização de Gram-Schimidt para vetores de dimensão 5
----------------------------------
README:
  - O programa recebe um arquivo como entrada e escreve no mesmo como saída;
  - Se o conjunto de vetores encontrados no arquivo de entrada for LD,
      ele informará e não apresentará vetores como saída. Caso contrário,
      os vetores ortogonalizados serão escritos ao final do arquivo.

  - SINOPSE DE COMPILAÇÃO:
      gcc trichez.c -lm

  - SINOPSE DE EXECUÇÃO:
      ./a.out entrada.txt

  - entrada.txt é um arquivo de texto contendo vetores de dimensão 5, segue exemplo:
      //Estes são os vetores que devem ser usados

      (1,2,3,4,5)

      (1,2,3,4,4)

      (1,1,3,4,3)

      (1,1,2,4,2)

      (1,1,1,1,1)

      Fim do arquivo

  - p.s: O programa procura por numeros no arquivo ignorando outros caracteres.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#define DIM 5                          //DIMENSÃO dos vetores. também usado paradeterminar o número de vetores na base

typedef struct vet{

  double * v;
  int dimensao;
  double norma;

}Vetor;

//------------------------------------------------------------

Vetor * iniciliazaVetor(Vetor *v);
double * projecao(Vetor *u, Vetor *v);  // retorna o vetor projeção
double calculaNorma(Vetor *v);          // calcula e seta a norma do vetor passado por parâmetro

/*---- /\ -------- manipulação de vetores ---------- /\ ------*/


int main(int argc, char const *argv[]) {

  int pch = NULL;
  int n = 0; // não começar em 0 (zero), pois ele é usado em divisão
  int i = 0;
  char c;
  char str[150];
  double *vetproj1 = NULL;
  double *vetproj2 = NULL;
  double *vetproj3 = NULL;
  double *vetproj4 = NULL;

  FILE *f = NULL;

  //instanciar vetores
  Vetor * v1 = NULL;
  Vetor * v2 = NULL;
  Vetor * v3 = NULL;
  Vetor * v4 = NULL;
  Vetor * v5 = NULL;
  v1 = iniciliazaVetor(v1);
  v2 = iniciliazaVetor(v2);
  v3 = iniciliazaVetor(v3);
  v4 = iniciliazaVetor(v4);
  v5 = iniciliazaVetor(v5);

  Vetor *proj1 = NULL;
  Vetor *proj2 = NULL;
  Vetor *proj3 = NULL;
  Vetor *proj4 = NULL;
  Vetor *proj5 = NULL;
  proj1 = iniciliazaVetor(proj1);
  proj2 = iniciliazaVetor(proj2);
  proj3 = iniciliazaVetor(proj3);
  proj4 = iniciliazaVetor(proj4);
  proj5 = iniciliazaVetor(proj5);

/*------ \/\/ ------ copia do arquivo p os vetores ------- \/\/ ----------*/

  f = fopen(argv[1] ,"r+");
  do{
    pch = getc(f);
    if(isdigit(pch)){
      c = pch;
      if(i < 5)
      v1->v[i] = atoi(&c);
      else if(i < 10)
      v2->v[i%5] = atoi(&c);
      else if(i < 15)
      v3->v[i%5] = atoi(&c);
      else if(i < 20)
      v4->v[i%5] = atoi(&c);
      else if(i < 25)
      v5->v[i%5] = atoi(&c);
      i++;
    }
  }while(pch != EOF);

  v1->norma = calculaNorma(v1);
  v2->norma = calculaNorma(v2);
  v3->norma = calculaNorma(v3);
  v4->norma = calculaNorma(v4);
  v5->norma = calculaNorma(v5);

  /*------- \/\/\/ ------ Gram-Schimidt ------ \/\/ ------------*/

  for (n = 0; n < DIM; n++) {
    proj1->v[n] = v1->v[n];                                                       //primeiro, só copia
  }
  proj1->norma = calculaNorma(proj1);

  vetproj1 =  projecao(v2,proj1);
  for (n = 0; n < DIM; n++) {
    proj2->v[n] = v2->v[n] - vetproj1[n];                                         // segundo
  }
  proj2->norma = calculaNorma(proj2);


  vetproj1 = projecao(v3,proj1);
  vetproj2 = projecao(v3,proj2);
  for (n = 0; n < DIM; n++) {
    proj3->v[n] = v3->v[n] - vetproj1[n] - vetproj2[n];                           //terceiro
  }
  proj3->norma = calculaNorma(proj3);

  vetproj1 = projecao(v4,proj1);
  vetproj2 = projecao(v4,proj2);
  vetproj3 = projecao(v4,proj3);
  for (n = 0; n < DIM; n++) {
    proj4->v[n] = v4->v[n] - vetproj1[n] - vetproj2[n] - vetproj3[n];             //quarto,
  }
  proj4->norma = calculaNorma(proj4);

  vetproj1 = projecao(v5,proj1);
  vetproj2 = projecao(v5,proj2);
  vetproj3 = projecao(v5,proj3);
  vetproj4 = projecao(v5,proj4);
  for (n = 0; n < DIM; n++) {
    proj5->v[n] = v5->v[n] - vetproj1[n] - vetproj2[n] - vetproj3[n] - vetproj4[n];//quinto vetor ortogonalizado
  }
  proj5->norma = calculaNorma(proj5);

  /* ---------- \/\/\/\/ -------- parte de escrita no arquivo ----------- \/\/\/\/ -------------*/
  
  if( proj1->norma == 0 || proj2->norma == 0 || proj3->norma == 0 || proj4->norma == 0 || proj5->norma == 0)
    fputs("\n============== O CONJUNTO ACIMA É LD ==============",f);
  else{
    fputs("\n--------------------------------------\n(",f);
    for (n = 0; n < DIM; n++) {
      fprintf(f,"%f, ", proj1->v[n]);
    }
    fprintf(f,") norma: [%f]\n", proj1->norma);

    fputs("\n--------------------------------------\n(",f);
    for (n = 0; n < DIM; n++) {
      fprintf(f,"%f, ", proj2->v[n]);
    }
    fprintf(f,") norma: [%f]\n", proj2->norma);

    fputs("\n--------------------------------------\n(",f);
    for (n = 0; n < DIM; n++) {
      fprintf(f,"%f, ", proj3->v[n]);
    }
    fprintf(f,") norma: [%f]\n", proj3->norma);

    fputs("\n--------------------------------------\n(",f);
    for (n = 0; n < DIM; n++) {
      fprintf(f,"%f, ", proj4->v[n]);
    }
    fprintf(f,") norma: [%f]\n", proj4->norma);

    fputs("\n--------------------------------------\n(",f);
    for (n = 0; n < DIM; n++) {
      fprintf(f,"%f, ", proj5->v[n]);
    }
    fprintf(f,") norma: [%f]\n", proj5->norma);
  }

  fclose(f);
  return 0;
}


Vetor * iniciliazaVetor(Vetor *v){

  v = malloc(sizeof(Vetor));
  v->v = malloc(DIM * sizeof(double));

  v->dimensao = DIM;
  v->norma = 0;

  return v;
}//DEVOLVE o vetor inicializado

double * projecao(Vetor *u, Vetor *v){                 // projeção de u sobre v
  int i;
  double delta = 0;
  double * proj = malloc(DIM * sizeof(double));

  for(i = 0; i < DIM; i++)
    delta += (v->v[i] * u->v[i]);                       // produto interno (escalar)

  for(i = 0; i < DIM; i++)
    proj[i] = (delta / (pow(v->norma, 2))) * v->v[i];

  return proj;
} //DEVOLVE o vetor projeção

double calculaNorma(Vetor *v){
  int i;

  for(i = 0; i < v->dimensao; i++)
  v->norma += pow(v->v[i],2);
  v->norma = sqrt(v->norma);

  return v->norma;
}//DEVOLVE a norma calculada
