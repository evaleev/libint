/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint compiler.
 *
 *  Libint compiler is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint compiler is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Libint compiler.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <limits.h>
#include <smart_ptr.h>

#include <cstdint>
#include <list>

#ifndef _libint2_src_bin_libint_memory_h_
#define _libint2_src_bin_libint_memory_h_

namespace libint2 {

/**
   MemoryBlock<Address,Size> describes a block of raw memory addressed via
   Address and size described by Size
 */
template <typename A, typename S>
class MemoryBlock {
 public:
  typedef A Address;
  typedef S Size;

  MemoryBlock(const Address& address, const Size& size, bool free,
              const std::shared_ptr<MemoryBlock>& left,
              const std::shared_ptr<MemoryBlock>& right)
      : address_(address),
        size_(size),
        free_(free),
        left_(left),
        right_(right) {}
  MemoryBlock(const MemoryBlock& other)
      : address_(other.address_),
        size_(other.size_),
        free_(other.free_),
        left_(other.left_),
        right_(other.right_) {}

  ~MemoryBlock() {}

  /// copy A to this
  const MemoryBlock& operator=(const MemoryBlock& other) {
    address_ = other.address_;
    size_ = other.size_;
    free_ = other.free_;
    left_ = other.left_;
    right_ = other.right_;
    return *this;
  }

  /// Returns address
  Address address() const { return address_; }
  /// Returns size
  Size size() const { return size_; }
  /// Returns true if the block is free
  bool free() const { return free_; }
  /// Returns the left adjacent block
  std::shared_ptr<MemoryBlock> left() const { return left_; }
  /// Returns the right adjacent block
  std::shared_ptr<MemoryBlock> right() const { return right_; }
  /// Sets the left adjacent block
  void left(const std::shared_ptr<MemoryBlock>& l) { left_ = l; }
  /// Sets the right adjacent block
  void right(const std::shared_ptr<MemoryBlock>& r) { right_ = r; }

  /// Sets the address
  void set_address(const Address& address) { address_ = address; }
  /// Sets the size
  void set_size(const Size& size) { size_ = size; }
  /// Sets block's free status
  void set_free(bool free) { free_ = free; }

  /// Returns true if the size of *i is less than the size of *j
  static bool size_less_than(const std::shared_ptr<MemoryBlock>& i,
                             const std::shared_ptr<MemoryBlock>& j) {
    return i->size() < j->size();
  }
  /** Returns true if the size of *i equals sz. Note that the arguments are
      not passed by reference since this function is designed to be converted
      to std::pointer_to_binary_function, which adds references to the arguments
   */
  static bool size_eq(std::shared_ptr<MemoryBlock> i, Size sz) {
    return i->size() == sz;
  }
  /** Returns true if the size of *i greater or equal than sz. Note that the
     arguments are not passed by reference since this function is designed to be
     converted to std::pointer_to_binary_function, which adds references to the
     arguments */
  static bool size_geq(std::shared_ptr<MemoryBlock> i, Size sz) {
    return i->size() >= sz;
  }
  /// Returns true if the address of *i is less than the address of *j
  static bool address_less_than(const std::shared_ptr<MemoryBlock>& i,
                                const std::shared_ptr<MemoryBlock>& j) {
    return i->address() < j->address();
  }
  /** Returns true if the address of *i equals a. Note that the arguments are
      not passed by reference since this function is designed to be converted
      to std::pointer_to_binary_function, which adds references to the arguments
   */
  static bool address_eq(std::shared_ptr<MemoryBlock> i, Address a) {
    return i->address() == a;
  }
  /// Returns true if *i is free
  static bool is_free(const std::shared_ptr<MemoryBlock>& i) {
    return i->free();
  }

  /// Merge A to this (does not check if merge can happen -- can_merge(*this,*A)
  /// must be already satisfied). The left/right pointers are not changed
  const MemoryBlock& merge(const MemoryBlock& other) {
    if (address() > other.address()) {
      address_ = other.address_;
    }
    size_ += other.size();
    return *this;
  }

 private:
  Address address_;
  Size size_;
  bool free_;
  typedef MemoryBlock<Address, Size> this_type;
  std::shared_ptr<this_type> left_;
  std::shared_ptr<MemoryBlock> right_;

  MemoryBlock();
};

/**
   Class MemoryManager handles allocation and deallocation of raw
   memory (stack) provided at runtime of the library.
*/

class MemoryManager {
 public:
  /// Negative Address is used to denote an invalid address -- hence signed
  /// integer
  typedef intptr_t Address;
  typedef size_t Size;
  typedef MemoryBlock<Address, Size> MemBlock;

  static const Address InvalidAddress = -1;

 protected:
  typedef std::list<std::shared_ptr<MemBlock> > memblkset;

 private:
  /// Upper limit on the amount of memory this MemoryManager can handle
  Size maxmem_;
  /// manages MemBlocks
  memblkset blks_;
  /// This block is guaranteed to be free until all memory is exhausted
  std::shared_ptr<MemBlock> superblock_;
  /// Max amount of memory used
  Size max_memory_used_;

  std::shared_ptr<MemBlock> merge_blocks(
      const std::shared_ptr<MemBlock>& left,
      const std::shared_ptr<MemBlock>& right);
  std::shared_ptr<MemBlock> merge_to_superblock(
      const std::shared_ptr<MemBlock>& blk);
  void update_max_memory();

 public:
  virtual ~MemoryManager();

  /// Reserve a block and return its address
  virtual Address alloc(const Size& size) = 0;
  /// Release a block previously reserved using alloc
  virtual void free(const Address& address);
  /// Returns the max amount of memory used up to this moment
  Size max_memory_used() const { return max_memory_used_; }

  /// resets the state of MemoryManager; does not invalidate stats, however
  void reset();

 protected:
  MemoryManager(const Size& maxmem);

  /// Returns maxmem
  Size maxmem() const { return maxmem_; }
  /// Returns blocks
  memblkset& blocks() { return blks_; }
  /// Returns the superblock
  std::shared_ptr<MemBlock> superblock() const { return superblock_; }
  /// steals size memory from block blk and returns the new block
  std::shared_ptr<MemBlock> steal_from_block(
      const std::shared_ptr<MemBlock>& blk, const Size& size);
  /// finds the block at Address a
  std::shared_ptr<MemBlock> find_block(const Address& a);
};

/**
   WorstFitMemoryManager allocates memory by trying to find the largest-possible
   free block. If search_exact == true -- exact fit is sought first.
*/
class WorstFitMemoryManager : public MemoryManager {
 public:
  WorstFitMemoryManager(bool search_exact = true,
                        const Size& maxsize = ULONG_MAX);
  virtual ~WorstFitMemoryManager();

  /// Implementation of MemoryManager::alloc()
  Address alloc(const Size& size) override;

 private:
  /// If seacrh_exact_ == true -- look for exact fit first
  bool search_exact_;
};

/**
   BestFitMemoryManager allocates memory by trying to find
   a suitable free block, which is is larger than the requested
   amount by at least tight_fit.
   If search_exact == true -- exact fit is sought first.
*/
class BestFitMemoryManager : public MemoryManager {
 public:
  BestFitMemoryManager(bool search_exact = true, const Size& tight_fit = 0,
                       const Size& maxsize = ULONG_MAX);
  virtual ~BestFitMemoryManager();

  /// Implementation of MemoryManager::alloc()
  Address alloc(const Size& size) override;

 private:
  /// If seacrh_exact_ == true -- look for exact fit first
  bool search_exact_;
  /// If size of a block - requested size < tight_fit_ it is rejected
  Size tight_fit_;
};

/**
   FirstFitMemoryManager allocates memory by finding first suitable free block.
   If search_exact == true -- exact fit is sought first.
*/
class FirstFitMemoryManager : public MemoryManager {
 public:
  FirstFitMemoryManager(bool search_exact = true,
                        const Size& maxsize = ULONG_MAX);
  virtual ~FirstFitMemoryManager();

  /// Implementation of MemoryManager::alloc()
  Address alloc(const Size& size) override;

 private:
  /// If seacrh_exact_ == true -- look for exact fit first
  bool search_exact_;
};

/**
   LastFitMemoryManager allocates memory by finding last suitable free block.
   If search_exact == true -- exact fit is sought first (from the back of the
   list).
*/
class LastFitMemoryManager : public MemoryManager {
 public:
  LastFitMemoryManager(bool search_exact = true,
                       const Size& maxsize = ULONG_MAX);
  virtual ~LastFitMemoryManager();

  /// Implementation of MemoryManager::alloc()
  Address alloc(const Size& size) override;

 private:
  /// If seacrh_exact_ == true -- look for exact fit first
  bool search_exact_;
};

/**
   MemoryManagerFactory is a very dumb factory for MemoryManagers
*/
class MemoryManagerFactory {
 public:
  static const unsigned int ntypes = 8;
  std::shared_ptr<MemoryManager> memman(unsigned int type) const;
  std::string label(unsigned int type) const;
};

/**
   Very useful nonmember functions to operate on MemBlocks and their containers
*/
typedef MemoryManager::MemBlock MemBlock;
typedef std::list<MemoryManager::MemBlock> MemBlockSet;
bool size_lessthan(const MemoryManager::MemBlock& A,
                   const MemoryManager::MemBlock& B);
bool address_lessthan(const MemoryManager::MemBlock& A,
                      const MemoryManager::MemBlock& B);
/// True if can merge blocks
bool can_merge(const MemoryManager::MemBlock& A,
               const MemoryManager::MemBlock& B);
/// Merge blocks, if possible
void merge(MemBlockSet& blocks);

};  // namespace libint2

#endif
