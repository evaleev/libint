
#include <limits.h>
#include <vector>
#include <smart_ptr.h>

#ifndef _libint2_src_bin_libint_memory_h_
#define _libint2_src_bin_libint_memory_h_

using namespace std;

namespace libint2 {

    /**
       MemoryBlock<Address,Size> describes a block of raw memory addressed via Address and size described by Size
     */
  template <typename Address, typename Size>
    class MemoryBlock
    {
    public:
      MemoryBlock(Address address, Size size, bool free,
                  const SafePtr<MemoryBlock>& prev,
                  const SafePtr<MemoryBlock>& next) :
        address_(address), size_(size), free_(free),
        prev_(prev), next_(next)
        {
        }

      ~MemoryBlock() {}

      /// Returns address
      Address address() const { return address_; }
      /// Returns size
      Size size() const { return size_; }
      /// Returns true if the block is free
      bool free() const { return free_; }

      /// Sets the address
      void set_address(Address address) { address_ = address; }
      /// Sets the size
      void set_size(Size size) { size_ = size; }
      /// Sets block's free status
      void set_free(bool free) { free_ = free; }

      /// Returns true if the size of *i is less than the size of *j
      static bool size_less_than(const SafePtr<MemoryBlock>& i,
                                 const SafePtr<MemoryBlock>& j) {
        return i->size() < j->size();
      }
      /// Returns true if the address of *i is less than the address of *j
      static bool address_less_than(const SafePtr<MemoryBlock>& i,
                                    const SafePtr<MemoryBlock>& j) {
        return i->address() < j->address();
      }
      /// Returns true if *i is free
      static bool is_free(const SafePtr<MemoryBlock>& i) {
        return i->free();
      }

    private:
      Address address_;
      Size size_;
      bool free_;
      SafePtr<MemoryBlock> prev_;
      SafePtr<MemoryBlock> next_;

      MemoryBlock();
    };

  /**
     Class MemoryManager handles allocation and deallocation of raw
     memory (stack) provided at runtime of the library.
  */

  class MemoryManager {
  public:
    typedef unsigned long int Address;
    typedef unsigned long int Size;
    typedef MemoryBlock<Address,Size> MemBlock;

    virtual ~MemoryManager();

    /// Reserve a block and return its address
    virtual Address alloc(Size size) =0;
    /// Release a block previously reserved using alloc
    virtual void free(Address address) =0;

  protected:
    MemoryManager(Size maxmem);

    typedef std::vector< SafePtr<MemBlock> > blkstore;
    // manages MemBlocks
    blkstore blks_;
    unsigned int last_freed_;

    /// Returns maxmem
    Size maxmem() const { return maxmem_; }
    /// steals size memory from block blk and returns the new block
    SafePtr<MemBlock> steal_from_block(const SafePtr<MemBlock>& blk, Size size);
    /// finds the block at Address a
    SafePtr<MemBlock> find_block(Address a);

  private:
    Size maxmem_;

  };


  /**
     WorstFitMemoryManager is a worst-fit memory manager
  */
  class WorstFitMemoryManager : public MemoryManager {
  public:
    WorstFitMemoryManager(Size maxsize = ULONG_MAX);
    ~WorstFitMemoryManager();

    /// Implementation of MemoryManager::alloc()
    Address alloc(Size size);
    /// Implementation of MemoryManager::free()
    void free(Address address);

  };

};

#endif
