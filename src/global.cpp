//$Source$
//------------------------------------------------------------------------------
// global
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file global.cpp
 * @brief Programacion de operacion global.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#include "..\include\global.hpp"
using namespace std;

Matrix eopdata;
Matrix PC;
Matrix Cnm;
Matrix Snm;

void eop19620101(int c)
{
    eopdata = zeros(13, c);
    FILE *fid = fopen("../data/eop19620101.txt", "r");
    if (fid == NULL)
    {
        printf("Fail open file\n");
        exit(EXIT_FAILURE);
    }
    for (int j = 1; j <= c; j++)
    {
        fscanf(fid, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &(eopdata(1, j)), &(eopdata(2, j)), &(eopdata(3, j)),
            &(eopdata(4, j)), &(eopdata(5, j)), &(eopdata(6, j)),
            &(eopdata(7, j)), &(eopdata(8, j)), &(eopdata(9, j)),
            &(eopdata(10, j)), &(eopdata(11, j)), &(eopdata(12, j)),
            &(eopdata(13, j)));
    }
    fclose(fid);
}

void DE430Coeff(int f, int c)
{
    PC = zeros(f, c);
    FILE *fid = fopen("../data/DE430Coeff.txt", "r");
    if (fid == NULL)
    {
        printf("Fail open file\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 1; i <= f; i++)
    {
        for (int j = 1; j <= c; j++)
        {
            fscanf(fid, "%lf", &(PC(i, j)));
        }
    }
    fclose(fid);
}

void GGM03S(int c)
{
    Cnm = zeros(c, c);
    Snm = zeros(c, c);
    FILE *fid = fopen("../data/GGM03S.txt", "r");
    if (fid == NULL)
    {
        printf("Fail open file\n");
        exit(EXIT_FAILURE);
    }
    double x;
    for (int n = 0; n <= c - 1; n++)
    {
        for (int m = 0; m <= n; m++)
        {
            fscanf(fid, "%lf %lf %lf %lf %lf %lf",
                   &(x), &(x), &(Cnm(n + 1, m + 1)), &(Snm(n + 1, m + 1)), &(x), &(x));
        }
    }
    fclose(fid);
}