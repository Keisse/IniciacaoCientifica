#include <cstdlib>
#include <iostream>
#include "glpk.h"

#define INF (2000111222)

int nCidades, nEstradas;
int G[123][123];
struct Estrada{
  int ori, dest, dist;     
};
Estrada estradas[12345];

void leEstradas() {
    scanf("%d%d", &nCidades, &nEstradas);
    for (int ori = 1; ori <= nCidades; ++ori) {
        for (int dest = 1; dest <= nCidades; ++dest) {
            G[ori][dest] = INF;
        }
    }
    for(int estrada = 1; estrada <= nEstradas; ++estrada) {
            int ori, dest, dist;
            scanf("%d%d%d", &ori, &dest, &dist);
            estradas[estrada].ori = ori;
            estradas[estrada].dest = dest;
            estradas[estrada].dist = dist;
            G[ori][dest] = G[dest][ori] = dist;
    }    
}

void NADA() {
     glp_prob *lp;
     lp = glp_create_prob();
     glp_add_cols(lp, 3);
     glp_set_col_bnds(lp, 1, GLP_LO, 0.0, 3.0); // 0 <= x1
     glp_set_col_bnds(lp, 2, GLP_LO, 0.0, 2.0); // 0 <= x2
     glp_set_col_bnds(lp, 3, GLP_LO, 0.0, 0.0); // 0 <= x3

     glp_set_obj_dir(lp, GLP_MAX); //max
     glp_set_obj_coef(lp, 1, -3.0); // -3x1
     glp_set_obj_coef(lp, 2, 4.0); // +4x2
     glp_set_obj_coef(lp, 3, 11.0); // +11x3
     // max -3x1 + 4x2 + 11 x3.
     
     int indCol[123];
     double val[123];

     glp_add_rows(lp, 1);
     indCol[1] = 1; val[1] = 10; // 10*x1
     indCol[2] = 2; val[2] = 3; // 3*x2
     glp_set_row_bnds(lp, 1, GLP_UP, 0.0, 15.0);// <=15
     glp_set_mat_row(lp, 1, 2, indCol, val);//  10 x1 + 3 x2 <= 15

     glp_add_rows(lp, 1);
     indCol[1] = 3; val[1] = 9; // 9*x3
     indCol[2] = 1; val[2] = 7; // 7*x1
     glp_set_row_bnds(lp, 2, GLP_UP, 0.0, 38.0);// <=38
     glp_set_mat_row(lp, 2, 2, indCol, val);//  7x1+9x2<=38

     glp_add_rows(lp, 1);
     indCol[1] = 3; val[1] = 15; // 15*x3
     glp_set_row_bnds(lp, 3, GLP_LO, 0.0, 25.0);// >=25
     glp_set_mat_row(lp, 3, 1, indCol, val);//  15x3 >=25
     
     glp_set_col_kind(lp, 1, GLP_IV);// X1 EH INTEIRO
     glp_set_col_kind(lp, 2, GLP_IV);// X2 EH INTEIRO
     glp_set_col_kind(lp, 3, GLP_IV);// X3 EH INTEIRO
     
     glp_intopt(lp, NULL); // acha solucao com restricao de integralidade
//     glp_simplex(lp, NULL);

//     printf("Solucao Otima: %.3f\n", glp_get_obj_val(lp));     
//     printf("X1: %.3f\n", glp_get_col_prim(lp, 1));     
//     printf("X2: %.3f\n", glp_get_col_prim(lp, 2));     
//     printf("X3: %.3f\n", glp_get_col_prim(lp, 3));

     printf("Solucao Otima: %.3f\n", glp_mip_obj_val(lp));     
     printf("X1: %.3f\n", glp_mip_col_val(lp, 1));     
     printf("X2: %.3f\n", glp_mip_col_val(lp, 2));     
     printf("X3: %.3f\n", glp_mip_col_val(lp, 3));
          
//     for (int est = 1; est <= nEstradas; ++est) {
//          glp_set_col_bnds(lp, est, GLP_LO, 0.0, 0.0);
//     }
}

glp_prob * montarModeloInicial() {
     glp_prob *lp;
     lp = glp_create_prob();//cria problema
     
     glp_add_cols(lp, nEstradas); // cria uma variavel por estrada
     for(int est = 1; est <= nEstradas; ++est) { // para cada estrada
          glp_set_col_bnds(lp, est, GLP_DB, 0.0, 1.0); // estrada entre 0 e 1
          glp_set_col_kind(lp, est, GLP_BV); // estrada binaria
     }
     
     glp_set_obj_dir(lp, GLP_MIN); //MIN
     for(int est=1; est <= nEstradas; ++est) { // para cada estrada
          glp_set_obj_coef(lp, est, estradas[est].dist);
           // custo da estrada eh coeficiente da objetivo
     }
     
     for (int cid = 1; cid <= nCidades; ++cid) {  //para cada cidade
          int indCol[123];
          double val[123];
          int nCoef = 0;
          for (int est = 1; est < nEstradas; ++est) { //para cada estrada
              Estrada estrada = estradas[est]; 
              if (estrada.ori == cid || estrada.dest == cid) {
                               //se cidade toca estrada est
                 indCol[nCoef + 1] = est; // a est-esima estrada
                 val[nCoef + 1] = 1.0; // com coeficiente 1 na linha
                 nCoef++; //incremente numero de coeficiente
              }
          }
          glp_add_rows(lp, 1);
          glp_set_mat_row(lp, cid, nCoef, indCol, val); // adiciona coeficientes da linha
          glp_set_row_bnds(lp, cid, GLP_DB, 2.0, 2.0); // restringe linha = 2.
     }
}

void busca(int v, int aE[], int vA[]) {
     if (vA[v]) return;
     vA[v] = 1;
     for(int e = 1; e <= nEstradas; ++e) {
             if (aE[e]) {
                   Estrada estrada = estradas[e];
                   if (estrada.ori == v) {
                      busca(estrada.dest);                
                   } else if (estrada.dest == v) {
                     busca(estrada.ori);
                   }
             }
     }
}

void encontraVerticesAlcancaveis(int arestaEscolhidas[], int verticesAlcancados[]) {
     for (int v = 1;) verticesAlcancados[v] = 0;
     busca(1, arestaEscolhidas, verticesAlcancados);
} 

int main(int argc, char *argv[])
{
    leEstradas();
    glp_prob *lp = montarModeloInicial();

while(1){     
     glp_intopt(lp, NULL); // acha solucao com restricao de integralidade

     printf("Solucao Otima: %.3f\n", glp_mip_obj_val(lp));     
     printf("X1: %.3f\n", glp_mip_col_val(lp, 1));     
     printf("X2: %.3f\n", glp_mip_col_val(lp, 2));     
     printf("X3: %.3f\n", glp_mip_col_val(lp, 3));
     int arestaEscolhida[1234];
     for(int est = 1; est) {
             arestaEscolhida[est] = glp_mip_col_val(lp, est);
     }
     int verticesAlcancados[123];
     for( ) ; // para contar se todos os vertices foram alcancados
     // se tiverem sido, de um break e mostre a resposta;
     encontraVerticesAlcancaveis(arestaEscolhida, verticesAlcancados);

     glp_add_row(lp, 1);
          int indCol[123];
          double val[123];
          int nCoef = 0;
     for (int e = 1; e <= nArestas; ++e) {
         Estrada estrada = estradas[e];
         int nextremosAlcancados = verticesAlcancados[estrada.ori] + 
             verticesAlcancados[estrada.dest];
         if (nextremosAlcancados == 1) {
            indCol[nCoef + 1] = e;
            val[nCoef + 1] = 1.0;
         } 
     }
     glp_set_mat_row(lp, glp_get_num_rows(lp), nCoef, indCol, val);
     glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_LO, 2.0, 2.0);
    system("PAUSE");
    return EXIT_SUCCESS;
}
