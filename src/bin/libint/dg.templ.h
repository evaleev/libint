
#include <dg.h>
#include <dgvertex.h>

#ifndef _libint2_src_bin_libint_dgtempl_h_
#define _libint2_src_bin_libint_dgtempl_h_

using namespace std;


namespace libint2 {

  /// Apply RR to target
  template <class I, class RR>
    void
    DirectedGraph::append_target(const SafePtr<I>& target)
    {
      target->make_a_target();
      recurse<I,RR>(target);
    };

  /// Apply RR to target
  template <class I, class RR>
    void
    DirectedGraph::recurse(const SafePtr<I>& vertex)
    {
      SafePtr<DGVertex> dgvertex = add_vertex(vertex);
      if (dgvertex != vertex)
	return;
      
      SafePtr<RR> rr0(new RR(vertex));
      const int num_children = rr0->num_children();
      
      for(int c=0; c<num_children; c++) {
        
        SafePtr<DGVertex> child = rr0->child(c);
        SafePtr<DGArc> arc(new DGArcRel<RR>(vertex,child,rr0));
        vertex->add_exit_arc(arc);
        
	SafePtr<I> child_cast = dynamic_pointer_cast<I,DGVertex>(child);
	if (child_cast == 0)
	  throw std::runtime_error("DirectedGraph::recurse(const SafePtr<I>& vertex) -- dynamic cast failed, most probably this is a logic error!");
        recurse<I,RR>(child_cast);
        
      }
    };

  /// Apply RR recursively starting with vertex
  template <class RR>
    void
    DirectedGraph::recurse(const SafePtr<DGVertex>& vertex)
    {
      SafePtr<DGVertex> dgvertex = add_vertex(vertex);
      if (dgvertex != vertex)
	return;
      
      typedef typename RR::TargetType TT;
      SafePtr<TT> tptr = dynamic_pointer_cast<TT,DGVertex>(vertex);
      if (tptr == 0)
        return;
      
      SafePtr<RR> rr0(new RR(tptr));
      const int num_children = rr0->num_children();
      
      for(int c=0; c<num_children; c++) {
        
        SafePtr<DGVertex> child = rr0->child(c);
        SafePtr<DGArc> arc(new DGArcRel<RR>(vertex,child,rr0));
        vertex->add_exit_arc(arc);
        
        recurse<RR>(child);
      }
    };

  /// Apply RR to all classes already on the graph
  template <class RR>
    void
    DirectedGraph::apply_to_all()
    {
      typedef typename RR::TargetType TT;
      typedef vertices::const_iterator citer;
      typedef vertices::iterator iter;
      for(citer v=stack_.begin(); v!=stack_.end(); ++v) {
	ver_ptr& vptr = vertex_ptr(*v);
        if ((vptr)->num_exit_arcs() != 0)
          continue;
        SafePtr<TT> tptr = dynamic_pointer_cast<TT,DGVertex>(v);
        if (tptr == 0)
          continue;
      
        SafePtr<RR> rr0(new RR(tptr));
        const int num_children = rr0->num_children();
      
        for(int c=0; c<num_children; c++) {
        
          SafePtr<DGVertex> child = rr0->child(c);
          SafePtr<DGArc> arc(new DGArcRel<RR>(tptr,child,rr0));
          tptr->add_exit_arc(arc);
        
          recurse<RR>(child);
        
        }
      }
    }

  template <class RR>
    unsigned int
    DirectedGraph::num_children_on(const SafePtr<RR>& rr) const
    {
      unsigned int nchildren = rr->num_children();
      unsigned int nchildren_on_stack = 0;
      for(int c=0; c<nchildren; c++) {
	if (!vertex_is_on(rr->rr_child(c)))
          continue;
	else
	  nchildren_on_stack++;
      }

      return nchildren_on_stack;
    }
    
  template <DirectedGraph::DGVertexMethodPtr method>
    void
    DirectedGraph::apply_at(const SafePtr<DGVertex>& vertex) const
    {
      ((vertex.get())->*method)();
      typedef DGVertex::ArcSetType::const_iterator aciter;
      const aciter abegin = vertex->first_exit_arc();
      const aciter aend = vertex->plast_exit_arc();
      for(aciter a=abegin; a!=aend; ++a)
        apply_at<method>((*a)->dest());
    }

  template <class Method>
    void
    DirectedGraph::foreach(Method& m)
    {
      typedef vertices::const_iterator citer;
      typedef vertices::iterator iter;
      for(iter v=stack_.begin(); v!=stack_.end(); ++v) {
	ver_ptr& vptr = vertex_ptr(*v);
	m(vptr);
      }
    }

  template <class Method>
    void
    DirectedGraph::foreach(Method& m) const
    {
      typedef vertices::const_iterator citer;
      typedef vertices::iterator iter;
      for(citer v=stack_.begin(); v!=stack_.end(); ++v) {
	const ver_ptr& vptr = vertex_ptr(*v);
	m(vptr);
      }
    }

  template <class Method>
    void
    DirectedGraph::rforeach(Method& m)
    {
      typedef vertices::const_reverse_iterator criter;
      typedef vertices::reverse_iterator riter;
      for(riter v=stack_.rbegin(); v!=stack_.rend(); ++v) {
	ver_ptr& vptr = vertex_ptr(*v);
	m(vptr);
      }
    }

  template <class Method>
    void
    DirectedGraph::rforeach(Method& m) const
    {
      typedef vertices::const_reverse_iterator criter;
      typedef vertices::reverse_iterator riter;
      for(criter v=stack_.rbegin(); v!=stack_.rend(); ++v) {
	const ver_ptr& vptr = vertex_ptr(*v);
	m(vptr);
      }
    }

};


#endif

