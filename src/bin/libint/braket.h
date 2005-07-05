
#include <polyconstr.h>
#include <bfset.h>

#ifndef _libint2_src_bin_libint_braket_h_
#define _libint2_src_bin_libint_braket_h_



namespace libint2 {

  /** VectorBraket is a std::vector-based type that can be used as a BraSetType or a KetSetType parameter
      to construct an instance of GenIntegralSet. VectorBraket is the most general implementation of Braket
      concept, but is fairly heavy. For better efficiency use ArrayBraket.
  */
  template <class BFS> class VectorBraket {

  public:
    typedef BFS bfs_type;
    typedef BFS bfs_stor;
    typedef bfs_stor& bfs_ref;
    typedef const bfs_stor& bfs_cref;
    typedef vector< bfs_stor > BFSVector;
    typedef vector< BFSVector > BFSMatrix;
    typedef VectorBraket<typename BFS::iter_type> iter_type;
    typedef struct{} parent_type;

    /** This one is a very dangerous constructor -- do not to use it if at all possible.
      Provided only for compatibility for generic subiterator algorithms */
    VectorBraket();
    VectorBraket(const BFSMatrix&);
    VectorBraket(const VectorBraket&);
    ~VectorBraket();

    /// Comparison function
    bool operator==(const VectorBraket&) const;
    /// Returns pointer to the i-th function for particle p
    bfs_ref member(unsigned int p, unsigned int i);
    /// Returns pointer to the i-th function for particle p
    bfs_cref member(unsigned int p, unsigned int i) const;
    /// Returns pointer to the SubIterator for i-th BFS of particle p
    SubIterator* member_subiter(unsigned int p, unsigned int i) const;
    /// Sets i-th function for particle p
    void set_member(bfs_cref, unsigned int p, unsigned int i);
    /// Sets i-th function for particle p (does a dynamic cast inside)
    void set_member(const ConstructablePolymorphically&, unsigned int p, unsigned int i);
    /// Returns the number of BFS for particle p
    const unsigned int num_members(unsigned int p) const;
    /// Returns the number of particles
    const unsigned int num_part() const;

  private:
    
    BFSMatrix bfs_;

  };

  template <class BFS>
    VectorBraket<BFS>::VectorBraket() :
    bfs_(0)
    {
    }
  
  template <class BFS>
    VectorBraket<BFS>::VectorBraket(const BFSMatrix& bfs) :
    bfs_(bfs)
    {
    }

  template <class BFS>
    VectorBraket<BFS>::VectorBraket(const VectorBraket& a) :
    bfs_(a.bfs_)
    {
    }

  template <class BFS>
    VectorBraket<BFS>::~VectorBraket()
    {
    }

  template <class BFS>
    typename VectorBraket<BFS>::bfs_ref
    VectorBraket<BFS>::member(unsigned int p, unsigned int i)
    {
      return bfs_.at(p).at(i);
    }
  
  template <class BFS>
    typename VectorBraket<BFS>::bfs_cref
    VectorBraket<BFS>::member(unsigned int p, unsigned int i) const
    {
      return bfs_.at(p).at(i);
    }

  template <class BFS>
    SubIterator*
    VectorBraket<BFS>::member_subiter(unsigned int p, unsigned int i) const
    {
      return static_cast<SubIterator*>(new SubIteratorBase<BFS>( member(p,i) ) );
    }
  
  template <class BFS>
    void
    VectorBraket<BFS>::set_member(bfs_cref bfs, unsigned int p, unsigned int i)
    {
      if (p >= bfs_.size())
        bfs_.resize(p+1);
      if (i >= bfs_[p].size())
        bfs_[p].resize(i+1);
      BFS bfs_tmp(bfs);
      bfs_[p][i] = bfs_tmp;
    }

  template <class BFS>
    void
    VectorBraket<BFS>::set_member(const ConstructablePolymorphically& bfs, unsigned int p, unsigned int i)
    {
      // WARNING : can be VERY dangerous
      // try constructing BFS from bfs.
      BFS bfs_cast(bfs);
      
      if (p >= bfs_.size())
        bfs_.resize(p+1);
      if (i >= bfs_[p].size())
        bfs_[p].resize(i+1);
      bfs_[p][i] = bfs_cast;
    }
  
  template <class BFS>
    const unsigned int
    VectorBraket<BFS>::num_members(unsigned int p) const
    {
      return bfs_.at(p).size();
    }

  template <class BFS>
    const unsigned int
    VectorBraket<BFS>::num_part() const
    {
      return bfs_.size();
    }

  template <class BFS>
    bool
    VectorBraket<BFS>::operator==(const VectorBraket<BFS>& a) const
    {
      const BFSMatrix& bfs0 = bfs_;
      const BFSMatrix& bfs1 = a.bfs_;
      if (bfs0.size() != bfs1.size())
        return false;

      // compare each row
      const int size1 = bfs_.size();
      for(int i=0; i<size1; i++) {
        const BFSVector& row0 = bfs0[i];
        const BFSVector& row1 = bfs1[i];

        if (row0.size() != row1.size())
          return false;

        // compare each element
        const int size2 = row0.size();
        for(int j=0; j<size2; j++)
          if (!PtrEquiv<BFS>::equiv(row0[j],row1[j]))
            return false;
      }
      return true;
    }

  //////////////////////////////////

  /** ArrayBraket is a lightweight implementation of Braket concept.
      BFS specifies the BasisFunctionSet. NP is the number of particles (>= 1).
      ArrayBraket assumes 1 BFS per particle per bra/ket.
  */
  template <class BFS, unsigned int NP> class ArrayBraket {
  public:
    /// This type
    typedef ArrayBraket<BFS,NP> this_type;
    /// There's no parent
    typedef struct{} parent_type; 
    /// The iterator through ArrayBraket
    typedef ArrayBraket<typename BFS::iter_type,NP> iter_type;

    /** This one is a very dangerous constructor -- do not to use it if at all possible.
        Provided only for compatibility for generic subiterator algorithms */
    ArrayBraket();
    ArrayBraket(const BFS* bra, const BFS* ket);
    ~ArrayBraket() throw();

    /// Comparison function
    bool operator==(const this_type&) const;
    /// Returns pointer to the i-th function for particle p
    const BFS& member(unsigned int p, unsigned int i=0) const;
    /// Returns pointer to the SubIterator for i-th BFS of particle p
    SubIterator* member_subiter(unsigned int p, unsigned int i=0) const;
    /// Sets i-th function for particle p
    void set_member(const BFS&, unsigned int p, unsigned int i=0);
    /// Sets i-th function for particle p (does a dynamic cast inside)
    void set_member(const ConstructablePolymorphically&, unsigned int p, unsigned int i=0);
    /// Returns the number of BFS for particle p
    const unsigned int num_members(unsigned int p) const { return 1; }
    /// Returns the number of particles
    const unsigned int num_part() const { return NP; }
    

  private:
    BFS bra_[NP];
    BFS ket_[NP];
    
  };
  
};

#endif
