#include <libint2.h>

void test_c_api() {
  libint2_static_init();

  libint2_static_cleanup();
}