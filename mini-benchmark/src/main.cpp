#include "test_range.hpp"
#include "test_mdrange.hpp"
#include "test_stencil.hpp"

int main (int argc, char *argv[]) {

  int TIMES       = 1;
  int START_NUM   = 1024;
  int RUNTIMES    = 1;
  int PRINTMATRIX = 0;

  printf ("argc = %d\n", argc);
  if (argc == 5) {
    TIMES       = atoi (argv[1]);
    START_NUM   = atoi (argv[2]);
    RUNTIMES    = atoi (argv[3]);
    PRINTMATRIX = atoi (argv[4]);
  }

  using T = double;
  // using T = float;
  // test_range_axpy<T> ();
  // test_mdrange_for_2d<T> ();
  // test_mdrange_for_3d<T> ();
  test_dot_product_2d<T> ();
  // test_5_point_stencil<T> ();
  printf ("==============================================\n");
  printf ("End main\n");

  return 0;
}
