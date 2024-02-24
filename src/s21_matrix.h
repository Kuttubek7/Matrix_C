#ifndef SRC_MATRIX_H_
#define SRC_MATRIX_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define SUCCESS 1
#define FAILURE 0

typedef struct matrix_struct {
  double **matrix;
  int rows;     // ряды
  int columns;  // столбцы
} matrix_t;

// создание матриц
int s21_create_matrix(int rows, int columns, matrix_t *result);
// очищает матрицу
void s21_remove_matrix(matrix_t *A);
// сравнение матриц
int s21_eq_matrix(matrix_t *A, matrix_t *B);
// сложение матрицы
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
// вычитание матрицы
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
// умножение матрицы на число
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
// умножение двух матриц
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
// транспонирование матрицы
int s21_transpose(matrix_t *A, matrix_t *result);
// определитель матрицы
int s21_determinant(matrix_t *A, double *result);
// минор и алгебраическое дополнение
int s21_calc_complements(matrix_t *A, matrix_t *result);
// обратная матрица
int s21_inverse_matrix(matrix_t *A, matrix_t *result);

// ДОПОЛНИТЕЛЬНЫЕ ФУНКЦИИ
// заполняет рандомными числами
int s21_zapolnenie(matrix_t *result);
// принтует
void s21_output(matrix_t *result);
// проверяет существует ли матрица
int s21_check_mtrx(matrix_t *A);
// инициализирует нужные матрицы
void s21_init(matrix_t *myMat, int num);
int s21_summator(int num, int pow);
// Убирает ряд и столбец из матрицы A и записывает в result
void s21_small_mtrx(matrix_t A, int row, int column, matrix_t *result);
// считает определитель
double s21_recursion_det(matrix_t A);

#endif  // SRC_MATRIX_H_