
#include <algorithm>
#include <list>
#include <stdexcept>
#include <iostream>
#include <memory.h>

using namespace std;
using namespace libint2;

MemoryManager::MemoryManager(const Size& maxmem) :
  maxmem_(maxmem), blks_()
{
  SafePtr<MemBlock> null_ptr;
  SafePtr<MemBlock> default_block(new MemBlock(Address(0),maxmem_,true,null_ptr,null_ptr));
  blks_.push_back(default_block);
  last_freed_ = 0;
}

MemoryManager::~MemoryManager()
{
}

SafePtr<MemoryManager::MemBlock>
MemoryManager::steal_from_block(const SafePtr<MemBlock>& blk, const Size& size)
{
  if (!blk->free())
    throw std::runtime_error("MemoryManager::steal_from_block() -- block is not free");

  Size old_size = blk->size();
  if (old_size <= size)
    throw std::runtime_error("MemoryManager::steal_from_block() -- block is too small");

  Size new_size = old_size - size;
  Address address = blk->address();
  blk->set_size(new_size);
  blk->set_address(address+size);
  SafePtr<MemBlock> null_ptr;
  SafePtr<MemBlock> newblk(new MemBlock(address,size,false,null_ptr,null_ptr));
  blks_.push_back(newblk);
  return newblk;
}

SafePtr<MemoryManager::MemBlock>
MemoryManager::find_block(const Address& address)
{
  typedef blkstore::iterator iter;
  iter begin = blks_.begin();
  iter end = blks_.end();
  for(iter b = begin; b != end; b++)
    if ( (*b)->address() == address)
      return *b;
  
  throw std::runtime_error("MemoryManager::find_block() -- didn't find a block at this address");
}

///////////////

WorstFitMemoryManager::WorstFitMemoryManager(const Size& maxsize) :
  MemoryManager(maxsize)
{
}

WorstFitMemoryManager::~WorstFitMemoryManager()
{
}

MemoryManager::Address
WorstFitMemoryManager::alloc(const Size& size)
{
  if (size > maxmem())
    throw std::runtime_error("WorstFitMemoryManager::alloc() -- requested more memory than available");

  // starting at block 1, try to find the exact match
  typedef blkstore::iterator iter;
  iter begin = blks_.begin();
  iter end = blks_.end();
  begin++;
  for(iter b=begin; b!=end; b++) {
    if((*b)->size() == size && (*b)->free()) {
      (*b)->set_free(false);
      return (*b)->address();
    }
  }

  // find all free_blocks
  std::list< SafePtr<MemBlock> > free_blks;
  for(iter b=begin; b!=end; b++) {
    b = find_if(b,end,&MemBlock::is_free);
    if (b != end)
      free_blks.push_back(*b);
    else // No more blocks left
      break;
  }

  // if no exact match found -- find the largest free block and grab memory from it
  std::list< SafePtr<MemBlock> >::iterator largest_free_block = max_element(free_blks.begin(),free_blks.end(),&MemBlock::size_less_than);
  if (largest_free_block != free_blks.end() &&
      (*largest_free_block)->size() > size) {
    SafePtr<MemBlock> result = steal_from_block(*largest_free_block,size);
    return result->address();
  }

  // lastly, if all failed -- steal from the master block
  SafePtr<MemBlock> result = steal_from_block(blks_[0],size);
  return result->address();
}

void
WorstFitMemoryManager::free(const Address& address)
{
  SafePtr<MemBlock> blk = find_block(address);
  if (!blk->free())
    blk->set_free(true);
  else
    throw std::runtime_error("WorstFitMemoryManager::free() tried to free a free block");
}
