
#define HAVE_STD_BINARY_COMPOSE 1

#include <algorithm>
#include <functional>
#include <list>
#include <stdexcept>
#include <iostream>
#include <memory.h>

using namespace std;
using namespace libint2;

#if HAVE_STD_BINARY_COMPOSE
#include <ext/functional>
using namespace __gnu_cxx;
#endif

MemoryManager::MemoryManager(const Size& maxmem) :
  maxmem_(maxmem), blks_(), superblock_(new MemBlock(Address(0),maxmem,true,SafePtr<MemBlock>(),SafePtr<MemBlock>())),
  max_memory_used_(0)
{
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
  SafePtr<MemBlock> left = blk->left();
  SafePtr<MemBlock> newblk(new MemBlock(address,size,false,left,blk));
  if (left)
    left->right(newblk);
  blk->left(newblk);
  blks_.push_back(newblk);

  update_max_memory();

  return newblk;
}

SafePtr<MemoryManager::MemBlock>
MemoryManager::find_block(const Address& address)
{
  typedef blkstore::iterator iter;
  iter blk = find_if(blks_.begin(),blks_.end(),bind2nd(ptr_fun(MemBlock::address_eq),address));
  if (blk == blks_.end())
    throw std::runtime_error("MemoryManager::find_block() -- didn't find a block at this address");
  else
    return *blk;
}

void
MemoryManager::free(const Address& address)
{
  SafePtr<MemBlock> blk = find_block(address);
  if (!blk->free())
    blk->set_free(true);
  else
    throw std::runtime_error("WorstFitMemoryManager::free() tried to free a free block");

  // Find blocks adjacent to this one and, if they are free, merge them
  SafePtr<MemBlock> left = blk->left();
  SafePtr<MemBlock> right = blk->right();
  if (left && left->free())
    blk = merge_blocks(left,blk);
  if (right && right->free())
    merge_blocks(blk,right);
}

SafePtr<MemoryManager::MemBlock>
MemoryManager::merge_blocks(const SafePtr<MemBlock>& left, const SafePtr<MemBlock>& right)
{
  if (left->free() != right->free())
    throw std::runtime_error("MemoryManager::merge_block() -- both blocks must be occupied or free");
  bool free = left->free();
  Address address = left->address();
  if (right->address() <= address)
    throw std::runtime_error("MemoryManager::merge_block() -- address of left block >= address of right block");
  if (left->address() + left->size() != right->address())
    throw std::runtime_error("MemoryManager::merge_block() -- address of left block + size of left block != address of right block");
  Size size = left->size() + right->size();

  if (right == superblock())
    return merge_to_superblock(left);
  else {
    SafePtr<MemBlock> lleft = left->left();
    SafePtr<MemBlock> rright = right->right();
    
    typedef blkstore::iterator iter;
    iter liter = find(blks_.begin(),blks_.end(),left);
    if (liter != blks_.end())
      blks_.erase(liter);
    else
      throw std::runtime_error("MemoryManager::merge_block() -- left block is not found");
    iter riter = find(blks_.begin(),blks_.end(),right);
    if (riter != blks_.end())
      blks_.erase(riter);
    else
      throw std::runtime_error("MemoryManager::merge_block() -- right block is not found");
  
    SafePtr<MemBlock> newblk(new MemBlock(address,size,free,lleft,rright));
    blks_.push_back(newblk);
    if (lleft) {
      lleft->right(newblk);
    }
    if (rright) {
      rright->left(newblk);
    }

    return newblk;
  }
}

SafePtr<MemoryManager::MemBlock>
MemoryManager::merge_to_superblock(const SafePtr<MemBlock>& blk)
{
  SafePtr<MemBlock> sblk = superblock();
  typedef blkstore::iterator iter;
  iter biter = find(blks_.begin(),blks_.end(),blk);
  if (biter != blks_.end())
    blks_.erase(biter);
  else
    throw std::runtime_error("MemoryManager::merge_to_superblock(blk) --  blk is not found");
  sblk->set_address(blk->address());
  sblk->set_size(sblk->size() + blk->size());
  SafePtr<MemBlock> left = blk->left();
  if (left)
    left->right(sblk);
  sblk->left(left);
  return sblk;
}

void
MemoryManager::update_max_memory()
{
  Address saddr =  superblock()->address();
  if (saddr > max_memory_used_)
    max_memory_used_  = saddr;
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
  if (size == 0)
    throw std::runtime_error("WorstFitMemoryManager::alloc(size) -- size is 0");

  // try to find the exact match first
  typedef blkstore::iterator iter;
  blkstore& blks = blocks();

#if HAVE_STD_BINARY_COMPOSE
  iter blk;
  blk = find_if(blks.begin(),blks.end(),
                compose2(logical_and<bool>(),
                         bind2nd(ptr_fun(MemBlock::size_eq),size),
                         &MemBlock::is_free));
  if (blk != blks.end()) {
    (*blk)->set_free(false);
    return (*blk)->address();
  }
#else
  iter begin = blks.begin();
  iter end = blks.end();
  for(iter b=begin; b!=end; b++) {
    if((*b)->size() == size && (*b)->free()) {
      (*b)->set_free(false);
      return (*b)->address();
    }
  }
#endif

  // find all free_blocks
  std::list< SafePtr<MemBlock> > free_blks;
  for(iter b=blks.begin(); b!=blks.end(); b++) {
    b = find_if(b,blks.end(),&MemBlock::is_free);
    if (b != blks.end())
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

  // lastly, if all failed -- steal from the super block
  SafePtr<MemBlock> result = steal_from_block(superblock(),size);
  return result->address();
}

///////////////

BestFitMemoryManager::BestFitMemoryManager(const Size& maxsize) :
  MemoryManager(maxsize)
{
}

BestFitMemoryManager::~BestFitMemoryManager()
{
}

MemoryManager::Address
BestFitMemoryManager::alloc(const Size& size)
{
  if (size > maxmem())
    throw std::runtime_error("BestFitMemoryManager::alloc() -- requested more memory than available");
  if (size == 0)
    throw std::runtime_error("BestFitMemoryManager::alloc(size) -- size is 0");

  // try to find the exact match first
  typedef blkstore::iterator iter;
  blkstore& blks = blocks();

#if HAVE_STD_BINARY_COMPOSE
  iter blk;
  blk = find_if(blks.begin(),blks.end(),
                compose2(logical_and<bool>(),
                         bind2nd(ptr_fun(MemBlock::size_eq),size),
                         &MemBlock::is_free));
  if (blk != blks.end()) {
    (*blk)->set_free(false);
    return (*blk)->address();
  }
#else
  iter begin = blks.begin();
  iter end = blks.end();
  for(iter b=begin; b!=end; b++) {
    if((*b)->size() == size && (*b)->free()) {
      (*b)->set_free(false);
      return (*b)->address();
    }
  }
#endif

  // find all free_blocks
  std::list< SafePtr<MemBlock> > free_blks;
  typedef std::list< SafePtr<MemBlock> >::iterator fiter;
  for(iter b=blks.begin(); b!=blks.end(); b++) {
    b = find_if(b,blks.end(),&MemBlock::is_free);
    if (b != blks.end())
      free_blks.push_back(*b);
    else // No more blocks left
      break;
  }

  // if there are no free blocks left -- steal from the super block
  if (free_blks.empty()) {
    SafePtr<MemBlock> result = steal_from_block(superblock(),size);
    return result->address();
  }

  // else find the smallest free block and grab memory from it
  fiter smallest_free_block = min_element(free_blks.begin(),free_blks.end(),&MemBlock::size_less_than);

  do {
    
    if ((*smallest_free_block)->size() > size) {
      SafePtr<MemBlock> result = steal_from_block(*smallest_free_block,size);
      return result->address();
    }
    else {
      free_blks.erase(smallest_free_block);
    }
    
    smallest_free_block = min_element(free_blks.begin(),free_blks.end(),&MemBlock::size_less_than);
    
  } while (smallest_free_block != free_blks.end());

  // Steal from superblock as a last resort
  SafePtr<MemBlock> result = steal_from_block(superblock(),size);
  return result->address();

}

