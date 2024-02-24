#include "s21_matrix.h"

// проверяет существует ли матрица
int s21_check_mtrx(matrix_t *A) {
  return A && A->matrix != NULL && (A->columns > 0 && A->rows > 0);
}

// принтует матрицу
// void s21_output(matrix_t *result) {
//   printf("  %d %d\n", result->rows, result->columns);
//   for (int i = 0; i < result->rows; i++) {
//     for (int j = 0; j < result->columns; j++) {
//       printf("%5.1f", result->matrix[i][j]);
//     }
//     printf("\n");
//   }
//   printf("\n");
// }

// считает определитель
double s21_recursion_det(matrix_t A) {
  matrix_t tmp = {0};
  double res = 0;
  if (A.rows == 1) {
    res = A.matrix[0][0];
  } else if (A.rows == 2) {
    res = (A.matrix[0][0] * A.matrix[1][1]) - (A.matrix[0][1] * A.matrix[1][0]);
  } else if (A.rows > 2) {
    int sign = 1;
    for (int i = 0; i < A.rows; i++) {
      s21_small_mtrx(A, i, 0, &tmp);
      double det = s21_recursion_det(tmp);
      res += sign * A.matrix[i][0] * det;
      sign *= -1;
      s21_remove_matrix(&tmp);
    }
  }
  return res;
}

// доп функция для инициализации
// int s21_summator(int num, int pow) {
//   int res = 0;
//   for (size_t i = 0; i < pow; i++) {
//     res += num;
//   }
//   return res;
// }

// инициализация нужных матриц
// void s21_init(matrix_t *myMat, int num) {
//   // int A[] = {1, 2, 3, 0, 4, 2, 5, 2, 1};
//   int A[] = {2, 5, 7, 6, 3, 4, 5, -2, -3};
//   // int A[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
//   // int A[] = {3, 1, 0, 1, 2, 1, 2, 4, 3};
//   int B[] = {-1, 4, -2, 5};

//   switch (num) {
//     case 'B':
//       for (int i = 0; i < myMat->rows; i++) {
//         for (int j = 0; j < myMat->columns; j++) {
//           myMat->matrix[i][j] = B[(summator(i, myMat->columns) + j)];
//         }
//       }
//       break;
//     case 'A':
//       for (int i = 0; i < myMat->rows; i++) {
//         for (int j = 0; j < myMat->columns; j++) {
//           myMat->matrix[i][j] = A[(summator(i, myMat->columns) + j)];
//         }
//       }
//       break;

//     default:
//       for (int i = 0; i < myMat->rows; i++) {
//         for (int j = 0; j < myMat->columns; j++) {
//           myMat->matrix[i][j] = rand() % 5;
//         }
//       }
//       break;
//   }
// }

// Убирает ряд и столбец из матрицы A и записывает в result
void s21_small_mtrx(matrix_t A, int row, int column, matrix_t *result) {
  int true_i = 0;
  int true_j = 0;
  if (!s21_create_matrix(A.rows - 1, A.columns - 1, result)) {
    for (int i = 0; i < A.rows; i++) {
      true_j = 0;
      if (i == row) {
        continue;
      }
      for (int j = 0; j < A.columns; j++) {
        if (j == column) {
          continue;
        }
        result->matrix[true_i][true_j] = A.matrix[i][j];
        true_j++;
      }
      true_i++;
    }
  }
}