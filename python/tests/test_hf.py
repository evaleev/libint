import libint2

import numpy as np
import scipy.linalg

def compute_1body_ints(oper, basis, params = None):
  engine = libint2.Engine(oper)
  if params: engine.set_params(params)
  return engine.compute_1body_ints(basis)

  # h = np.zeros([basis.nbf,basis.nbf])
  # bf = basis.functions
  # for i,p in enumerate(bf):
  #   for j,q in enumerate(bf):
  #     h[p,q] = engine.compute(*basis.shells(i,j))
  # return h


def compute_2body_fock(D, basis):
  engine = libint2.Engine(libint2.Operator.coulomb)
  return engine.compute_2body_fock(D, basis)

  # F = np.zeros([basis.nbf,basis.nbf])
  # bf = basis.functions
  # for i in range(len(basis)):
  #   for j in range(0,i+1):
  #     for k in range(len(basis)):
  #       for l in range(0,k+1):

  #         p,q,r,s = [bf[idx] for idx in (i,j,k,l)]
  #         eri = engine.compute(*basis.shells(i,j,k,l))
  #         if eri is None: continue

  #         symm = 1.0
  #         if i == j: symm /= 2
  #         if k == l: symm /= 2

  #         F[p,q] += 2*np.einsum('pqrs,rs->pq', eri, D[r,s])*2*symm

  #         F[p,r] -= np.einsum('pqrs,qs->pr', eri, D[q,s])/2*symm
  #         F[q,r] -= np.einsum('pqrs,ps->qr', eri, D[p,s])/2*symm

  #         F[p,s] -= np.einsum('pqrs,qr->ps', eri, D[q,r])/2*symm
  #         F[q,s] -= np.einsum('pqrs,pr->qs', eri, D[p,r])/2*symm

  #         # J = engine.compute(*basis.shells(i,k,j,l))
  #         # F[p,q] += 2*np.einsum('pqrs,rs->pq', J, D[r,s])
  #         # K = engine.compute(*basis.shells(i,k,j,l))
  #         # F[p,q] -= np.einsum('prqs,rs->pq', K, D[r,s])

  # return F+F.T

class RHF:

  def __init__(self, basis, atoms):
    if not isinstance(basis, libint2.BasisSet):
      basis = libint2.BasisSet(basis, atoms)
    assert basis and basis.nbf
    self.basis = basis

    self.nelec = sum([z for z,r in atoms])
    self.ndocc = self.nelec//2
    self.enuc = 0
    for i,(qi,ri) in enumerate(atoms):
      for j,(qj,rj) in enumerate(atoms[:i]):
        self.enuc += qi*qj/np.linalg.norm(np.array(ri)-np.array(rj), ord=2)

    self.S = compute_1body_ints(libint2.Operator.overlap, basis)
    self.T = compute_1body_ints(libint2.Operator.kinetic, basis)
    self.V = compute_1body_ints(libint2.Operator.nuclear, basis, atoms)
    self.H = self.V + self.T
    self.compute_density(self.H)
    self.converged = False
    self.F = None
    self.ehf = None

  def compute_density(self, F):
    ndocc = self.ndocc
    eig, C = scipy.linalg.eigh(F,self.S)
    D = np.matmul(C[:,:ndocc], C[:,:ndocc].T)
    self.C = C
    self.D = D

  def energy(self):
    if not self.converged:
      self.converge()
    return (self.enuc + self.ehf)

  def converge(self, iterations=30, tol=1e-6):
    self.converged = False
    for i in range(iterations):
      F = compute_2body_fock(self.D, self.basis)
      F += self.H
      ehf = self.ehf or 0.0
      self.compute_density(F)
      self.F = F
      self.ehf = np.sum(np.multiply(self.D, self.H+F))
      ediff = ehf-self.ehf
      #print("%i  %f  %f"%(i, (self.enuc + self.ehf), ediff))
      if abs(ediff) < tol:
        self.converged = True
        break


import unittest

h2o = [
  (8, [  0.00000, -0.07579, 0.00000 ]),
  (1, [  0.86681,  0.60144, 0.00000 ]),
  (1, [ -0.86681,  0.60144, 0.00000 ]),
]

class TestHF(unittest.TestCase):
  def test_hf(self):
    basis = '6-31g'
    rhf = RHF(basis, h2o)
    self.assertAlmostEqual(rhf.energy(), -75.1903033978)

if __name__ == '__main__':
  unittest.main()
      
