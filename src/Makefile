CC=gcc 
CFLAGS=-g -c -Wall -Wextra -Werror -std=c11
TEST_FLAGS=-lcheck -lm -lpthread
GCOV_FLAGS=--coverage -Wall -Wextra -Werror
SOURCES=$(wildcard s21_*.c)
OBJECTS=$(SOURCES:.c=.o)
TESTS_SOURCES=test_main.c $(wildcard unit_test/*.c)
TESTS_OBJECTS=$(TESTS_SOURCES:.c=.o)
EXECUTABLE=s21_matrix

# флаги для линукса
ifeq ($(shell uname -s),Linux)
	TEST_FLAGS += -lrt -lsubunit
endif

all: test

$(EXECUTABLE): $(OBJECTS)
	$(CC) $^ -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

rebuild: clean all

s21_matrix.a: $(OBJECTS)
	ar -rc s21_matrix.a $(OBJECTS)
	ranlib s21_matrix.a

style:
	mv ../materials/build/scripts/.clang-format ../src
	clang-format -i *.c *.h unit_test/*.c unit_test/*.h
	clang-format -n *.c *.h unit_test/*.c unit_test/*.h
	mv .clang-format ../materials/build/scripts

clean:
	rm -rf $(OBJECTS) $(EXECUTABLE) $(TESTS_OBJECTS) *.gcov *.gcno *.gcda *.info test *.dSYM a.out report s21_matrix.a .DS_Store

test: $(TESTS_OBJECTS) s21_matrix.a
	$(CC) $^ -o test $(TEST_FLAGS)
	./test

gcov_report: $(TESTS_OBJECTS)
	$(CC) $(GCOV_FLAGS) $(TESTS_OBJECTS) $(SOURCES) $(TEST_FLAGS) -o test
	./test
	gcov *.gcno
	lcov -t "s21_matrix" -o s21_matrix.info -c -d . --rc lcov_coverage
	genhtml -o report s21_matrix.info --rc lcov_coverage
	open report/index.html

main:
	$(CC) s21_matrix.c s21_dop_function.c
	./a.out 

# leaks: test
# 	leaks -atExit -- ./test

# vlagrind: test
# 	# valgring --leak-check=full ./test