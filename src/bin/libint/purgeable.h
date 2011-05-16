
#ifndef _libint2_src_bin_libint_purgeable_h_
#define _libint2_src_bin_libint_purgeable_h_

#include <vector>

namespace libint2 {

  /** Determines whether an object should be purged from a stack. The default policy is to
      purge DGVector if it doesn't belong to any DirectedGraph.
    */
  template <typename T>
  struct DefaultPurgingPolicy {
    /// returns true if objects of this type can be purged
    static bool purgeable();
    /// returns true if obj should be purged
    static bool purge(const T* obj);
  };


  /**
     * PurgeableStack is a container that can be purged by calling purge() method.
     */
  class AbstractPurgeableStack {
    public:
      virtual ~AbstractPurgeableStack();
      virtual void purge() =0;
  };

  /**
   * PurgeableStack is an AbstractPurgeableStack that contains objects of type T.
   * Whether to purge an object is determined by calling Policy::purge(const T*)
   */
  template <typename T, typename Policy = DefaultPurgingPolicy<T> >
    class PurgeableStack : public AbstractPurgeableStack
    {
      protected:
        typedef Policy PurgingPolicy;
    };

  /// Collection of AbstractPurgeableStack objects
  class PurgeableStacks {
    public:
      typedef PurgeableStacks this_type;
      typedef AbstractPurgeableStack stack_type;

      static this_type* Instance();

      void purge();
      void register_stack(stack_type* stack);

    private:
      PurgeableStacks() {}
      static this_type* instance_;
      std::vector<stack_type*> stacks_;
  };

};

#endif
