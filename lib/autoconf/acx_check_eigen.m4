AC_DEFUN([ACX_CHECK_EIGEN], [
  AC_LANG_SAVE
  AC_LANG([C++])

  acx_have_eigen=no
  AC_CHECK_HEADER(Eigen/Core,[
    AC_MSG_CHECKING([if Eigen library is usable])
    AC_LINK_IFELSE(
     [AC_LANG_PROGRAM(
      [[#include <Eigen/Dense>
        #include <Eigen/Eigenvalues>
      ]],
      [[
       typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
       Matrix h; h << 2.0, 1.0, 1.0, 2.0;
       Matrix s; s << 1.0, 0.0, 0.0, 1.0;
       Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(h, s);
       auto evals = gen_eig_solver.eigenvalues();
      ]]
      )
     ],
     [acx_have_eigen=yes
      AC_DEFINE(LIBINT_HAS_EIGEN)
     ]
    )
    AC_MSG_RESULT([$acx_have_eigen])
   ]
  )

  AC_LANG_RESTORE  
])
