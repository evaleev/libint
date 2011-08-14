
#ifndef _libint2_src_bin_libint_purgeable_h_
#define _libint2_src_bin_libint_purgeable_h_

#include <dgvertex.h>
#include <vector>
#include <boost/type_traits.hpp>

namespace libint2 {

  /** Determines whether an object should be purged from a stack. The default policy is to
      purge DGVector if it doesn't belong to any DirectedGraph.
    */
  template <typename T>
  struct DefaultPurgingPolicy {
    /// returns true if objects of this type can be purged
    static bool purgeable() {

      bool result = false;

      if (boost::is_base_of<DGVertex,T>::value == true) { // can only purge DGVertex objects
        result = true;
      }

      return result;
    }

    /// returns true if obj should be purged
    static bool purge(const T* ref) {

      bool result = false;

      try {
        const DGVertex* dgv_ptr = dynamic_cast<const DGVertex*>(ref);
        if (dgv_ptr->dg() == 0)
          result = true;
      }
      catch(...) {
      }

      return result;
    }
  };


  /**
     * PurgeableStack is a container that can be purged by calling purge() method.
     */
  class AbstractPurgeableStack {
    public:
      virtual ~AbstractPurgeableStack() {}
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

        virtual ~PurgeableStack() {}
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
