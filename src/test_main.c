#include <stdio.h>
#include <stdlib.h>

#include "unit_test/unit.h"

void run_test_case(Suite *testcase);
void run_all_tests();

int main() {
  run_all_tests();
  return 0;
}

void run_test_case(Suite *testcase) {
  static int count = 1;
  printf("\nTEST CASE №%d\n", count++);

  SRunner *sr;
  sr = srunner_create(testcase);

  srunner_set_fork_status(sr, CK_NOFORK);
  srunner_run_all(sr, CK_NORMAL);
  srunner_free(sr);
}

void run_all_tests() {
  Suite *list_cases[] = {
      suite_s21_create_matrix(), suite_s21_eq_matrix(),
      suite_s21_sum_matrix(),    suite_s21_sub_matrix(),
      suite_s21_mult_number(),   suite_s21_mult_matrix(),
      suite_s21_transpose(),     suite_s21_calc_complements(),
      suite_s21_determinant(),   suite_s21_inverse_matrix(),
  };
  for (int i = 0; i < (int)(sizeof(list_cases) / sizeof(Suite *)); i++)
    run_test_case(list_cases[i]);
}