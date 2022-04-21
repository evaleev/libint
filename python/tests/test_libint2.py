import unittest
import libint2
from numpy.linalg import norm

from libint2 import Shell, BasisSet

libint2.Engine.num_threads = 1

s = Shell(0, [(1,10)])
p = Shell(1, [(1,10)])

h2o = [
  (8, [  0.00000, -0.07579, 0.00000 ]),
  (1, [  0.86681,  0.60144, 0.00000 ]),
  (1, [ -0.86681,  0.60144, 0.00000 ]),
]

class TestLibint(unittest.TestCase):

  def test_core(self):
    self.assertTrue(libint2.MAX_AM > 0)    

    basis = BasisSet('6-31g', h2o)
    self.assertEqual(len(basis), 9)
    basis.pure = False
    pure = [False]*len(basis)
    self.assertEqual([ s.pure for s in basis], pure)
    basis[0].pure = True
    pure[0] = True
    self.assertEqual([ s.pure for s in basis], pure)
  
  def test_integrals(self):
    self.assertAlmostEqual(norm(libint2.kinetic().compute(s,s)), 1.5)
    self.assertAlmostEqual(norm(libint2.overlap().compute(s,s)), 1.0)
    self.assertAlmostEqual(norm(libint2.nuclear(h2o).compute(s,s)), 14.54704336519)

    self.assertAlmostEqual(norm(libint2.coulomb().compute(p,p,s,s)), 1.62867503968)

    self.assertAlmostEqual(
      norm(libint2.Engine(libint2.Operator.coulomb, braket=libint2.BraKet.XXXS).compute(s,s,s)),
      3.6563211198
    )

    basis = [ s, p, s, p ]
    self.assertAlmostEqual(libint2.overlap().compute(basis, basis).sum(), 16.0)
    self.assertAlmostEqual(
      norm(libint2.coulomb().compute(basis, basis, basis, basis)),
      14.7036075402
    )

if __name__ == '__main__':
  unittest.main()
