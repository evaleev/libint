
#include <functional>
#include <fstream>
#include <dg.h>
#include <dg.templ.h>
#include <rr.h>
#include <strategy.h>
#include <prefactors.h>
#include <codeblock.h>
#include <default_params.h>
#include <graph_registry.h>
#include <global_macros.h>

using namespace std;
using namespace libint2;

#define ONLY_CLONE_IF_DIFF 1

DirectedGraph::DirectedGraph() :
  stack_(default_size_,SafePtr<DGVertex>()), targets_(0), func_names_(),
  first_free_(0), first_to_compute_(), registry_(SafePtr<GraphRegistry>(new GraphRegistry)),
  iregistry_(SafePtr<InternalGraphRegistry>(new InternalGraphRegistry))
{
}

DirectedGraph::~DirectedGraph()
{
  reset();
}

void
DirectedGraph::append_target(const SafePtr<DGVertex>& target)
{
  target->make_a_target();
  append_vertex(target);
  targets_.push_back(target);
}

void
DirectedGraph::append_vertex(const SafePtr<DGVertex>& vertex) throw(VertexAlreadyOnStack)
{
  try {
    add_vertex(vertex);
    // If this is a new vertex -- tell the vertex who its owner is now
    vertex->dg(this);
  }
  // Handle the trivial case when this exact object is already on graph -- there's no reason to throw the exception then
  catch (VertexAlreadyOnStack& e) {
    //    if (e.vertex() != vertex)
      throw;
  }
}

void
DirectedGraph::add_vertex(const SafePtr<DGVertex>& vertex) throw(VertexAlreadyOnStack)
{
  vertex_is_on(vertex);
  add_new_vertex(vertex);
  return;
}

void
DirectedGraph::add_new_vertex(const SafePtr<DGVertex>& vertex)
{
  if (first_free_ == stack_.size()) {
    stack_.resize( stack_.size() + default_size_ );
#if DEBUG
    cout << "Increased size of DirectedGraph's stack to "
         << stack_.size() << endl;
#endif
  }
  char label[80];  sprintf(label,"vertex%d",first_free_);
  vertex->set_graph_label(label);
  stack_[first_free_++] = vertex;
#if DEBUG
  cout << "Added vertex " << vertex->description() << endl;
#endif
  return;
}

void
DirectedGraph::vertex_is_on(const SafePtr<DGVertex>& vertex) const throw(VertexAlreadyOnStack)
{
  for(int i=0; i<first_free_; i++)
    if(vertex->equiv(stack_[i]))
      throw VertexAlreadyOnStack(stack_[i]);
}

void
DirectedGraph::del_vertex(const SafePtr<DGVertex>& v) throw(CannotPerformOperation)
{
  // Cannot delete targets. Should I be able to? Probably not
  if (v->is_a_target())
    throw CannotPerformOperation("DirectedGraph::del_vertex() cannot delete targets");

  vector< SafePtr<DGVertex> >::iterator pos = find(stack_.begin(),stack_.end(),v);
  if (pos == stack_.end())
    throw CannotPerformOperation("DirectedGraph::del_vertex() cannot delete vertex");

  if (v->num_exit_arcs() == 0 && v->num_entry_arcs() == 0)
    stack_.erase(pos);
  else
    throw CannotPerformOperation("DirectedGraph::del_vertex() cannot delete vertex");
  --first_free_;
}

void
DirectedGraph::prepare_to_traverse()
{
  for(int i=0; i<first_free_; i++)
    stack_[i]->prepare_to_traverse();
}

void
DirectedGraph::traverse()
{
  // Initialization
  prepare_to_traverse();

  // Start at the targets which don't have parents
  for(int i=0; i<first_free_; i++) {
    if (stack_[i]->is_a_target() && stack_[i]->num_entry_arcs() == 0) {
      SafePtr<DGVertex> vertex_ptr = stack_[i];
      // First, since this target doesn't have parents we can schedule its computation
      schedule_computation(vertex_ptr);
      int nchildren = vertex_ptr->num_exit_arcs();
      for(int c=0; c<nchildren; c++)
        traverse_from(vertex_ptr->exit_arc(c));
    }
  }
}

void
DirectedGraph::traverse_from(const SafePtr<DGArc>& arc)
{
  SafePtr<DGVertex> orig = arc->orig();
  SafePtr<DGVertex> dest = arc->dest();
  // no need to compute if precomputed
  if (dest->precomputed())
    return;
  // if has been hit by all parents ...
  const unsigned int num_tags = dest->tag();
  const unsigned int num_parents = dest->num_entry_arcs();
  if (num_tags == num_parents) {

    // ... check if it is on a subtree ...
    if (SafePtr<DRTree> stree = dest->subtree()) {
      // ... if yes, schedule only if it is a root of a subtree
      if (stree->root() == dest)
        schedule_computation(dest);
    }
    // else schedule
    else
      schedule_computation(dest);
    
    int nchildren = dest->num_exit_arcs();
    for(int c=0; c<nchildren; c++)
      traverse_from(dest->exit_arc(c));
  }
}

void
DirectedGraph::schedule_computation(const SafePtr<DGVertex>& vertex)
{
  vertex->set_postcalc(first_to_compute_);
  first_to_compute_ = vertex;
}


void
DirectedGraph::debug_print_traversal(std::ostream& os) const
{
  SafePtr<DGVertex> current_vertex = first_to_compute_;

  os << "Debug print of traversal order" << endl;

  do {
    current_vertex->print(os);
    current_vertex = current_vertex->postcalc();
  } while (current_vertex != 0);
}

void
DirectedGraph::print_to_dot(bool symbols, std::ostream& os) const
{
  os << "digraph G {" << endl
     << "  size = \"8,8\"" << endl;
  for(int i=0; i<first_free_; i++) {
    SafePtr<DGVertex> vertex = stack_[i];
    os << "  " << vertex->graph_label()
       << " [ label = \"";
    if (symbols && vertex->symbol_set())
      os << vertex->symbol();
    else
      os << vertex->label();
    os << "\"]" << endl;
  }
  for(int i=0; i<first_free_; i++) {
    SafePtr<DGVertex> vertex = stack_[i];
    unsigned int narcs = vertex->num_exit_arcs();
    for(int a=0; a<narcs; a++) {
      SafePtr<DGVertex> dest = vertex->exit_arc(a)->dest();
      os << "  " << vertex->graph_label() << " -> "
         << dest->graph_label() << endl;
    }
  }

  // Print traversal order using dotted lines
  SafePtr<DGVertex> current_vertex = first_to_compute_;
  if (current_vertex != 0) {
    do {
      SafePtr<DGVertex> next = current_vertex->postcalc();
      if (current_vertex && next) {
        os << "  " << current_vertex->graph_label() << " -> "
           << next->graph_label() << " [ style = dotted ]";
      }
      current_vertex = next;
    } while (current_vertex != 0);
  }

  os << "}" << endl;
}

void
DirectedGraph::reset()
{
  // Reset each vertex, releasing all arcs
  for(int i=0; i<first_free_; i++) {
#if DEBUG
    std::cout << "DGStack::reset: will unregister " << stack_[i]->label() << std::endl;
#endif
    // remove this vertex from its SingletonManager
    stack_[i]->unregister();
    stack_[i]->reset();
  }
  for(int i=0; i<first_free_; i++)
    stack_[i].reset();
  // if everything went OK then resize stack_ to default and targets_ to 0
  stack_.resize(default_size_);
  targets_.resize(0);
  first_free_ = 0;
  first_to_compute_.reset();
}


/// Apply a strategy to all vertices not yet computed (i.e. which do not have exit arcs)
void
DirectedGraph::apply(const SafePtr<Strategy>& strategy,
                     const SafePtr<Tactic>& tactic)
{
  const int num_vertices_on_graph = first_free_;
  for(int v=0; v<num_vertices_on_graph; v++) {
    if (stack_[v]->num_exit_arcs() != 0 || stack_[v]->precomputed() || !stack_[v]->need_to_compute())
      continue;

    SafePtr<DirectedGraph> this_ptr = SafePtr_from_this();
    SafePtr<RecurrenceRelation> rr0 = strategy->optimal_rr(this_ptr,stack_[v],tactic);
    if (rr0 == 0)
      continue;

    // add children to the graph
    SafePtr<DGVertex> target = rr0->rr_target();
    const int num_children = rr0->num_children();
    for(int c=0; c<num_children; c++) {
      SafePtr<DGVertex> child = rr0->rr_child(c);
      bool new_vertex = true;
      try { append_vertex(child); }
      catch (VertexAlreadyOnStack& e) { child = e.vertex(); new_vertex = false; }
      SafePtr<DGArc> arc(new DGArcRel<RecurrenceRelation>(target,child,rr0));
      target->add_exit_arc(arc);
      if (new_vertex)
        apply_to(child,strategy,tactic);
    }
  }
}

/// Add vertex to graph and apply a strategy to vertex recursively
void
DirectedGraph::apply_to(const SafePtr<DGVertex>& vertex,
                        const SafePtr<Strategy>& strategy,
                        const SafePtr<Tactic>& tactic)
{
  if (vertex->precomputed() || !vertex->need_to_compute())
    return;
  SafePtr<RecurrenceRelation> rr0 = strategy->optimal_rr(SafePtr_from_this(),vertex,tactic);
  if (rr0 == 0)
    return;

  SafePtr<DGVertex> target = rr0->rr_target();
  const int num_children = rr0->num_children();
  for(int c=0; c<num_children; c++) {
    SafePtr<DGVertex> child = rr0->rr_child(c);
    bool new_vertex = true;
    try { append_vertex(child); }
    catch (VertexAlreadyOnStack& e) { child = e.vertex(); new_vertex = false; }
    SafePtr<DGArc> arc(new DGArcRel<RecurrenceRelation>(target,child,rr0));
    target->add_exit_arc(arc);
    if (new_vertex)
      apply_to(child,strategy,tactic);
  }
}

// Optimize out simple recurrence relations
void
DirectedGraph::optimize_rr_out()
{
  replace_rr_with_expr();
  remove_trivial_arithmetics();
  handle_trivial_nodes();
  remove_disconnected_vertices();
  find_subtrees();
}

// Replace recurrence relations with expressions
void
DirectedGraph::replace_rr_with_expr()
{
  const unsigned int nvertices = first_free_;
  for(int v=0; v<nvertices; v++) {

    SafePtr<DGVertex> vertex = stack_[v];
    if (vertex->num_exit_arcs()) {
      SafePtr<DGArc> arc0 = vertex->exit_arc(0);
      SafePtr<DGArcRR> arc0_cast = dynamic_pointer_cast<DGArcRR,DGArc>(arc0);
      if (arc0_cast == 0)
        continue;
      SafePtr<RecurrenceRelation> rr = arc0_cast->rr();

      // Optimize if the recurrence relation is simple and the target and
      // children are of the same type
      if (rr->is_simple() && rr->invariant_type()) {

        unsigned int nchildren = rr->num_children();

        // Remove arcs connecting this vertex to children
        vertex->del_exit_arcs();

        // and instead insert the numerical expression
        SafePtr<RecurrenceRelation::ExprType> rr_expr = rr->rr_expr();
        SafePtr<DGVertex> expr_vertex = static_pointer_cast<RecurrenceRelation::ExprType,DGVertex>(rr_expr);
	try {
          insert_expr_at(vertex,rr_expr);
	}
	catch (VertexAlreadyOnStack& e) { expr_vertex = e.vertex(); }
        SafePtr<DGArc> arc(new DGArcDirect(vertex,expr_vertex));
        vertex->add_exit_arc(arc);

      }
    }
  }
}


//
// This function is very tricky at the moment. The operands have to be added before the operator
// such that operands are guaranteed to be on graph before the operator. This way ExprType::equiv
// can simply compare operand pointers and the operator type (very cheap operations compared
// to the fully recursive explicit comparison).
//
void
DirectedGraph::insert_expr_at(const SafePtr<DGVertex>& where, const SafePtr<RecurrenceRelation::ExprType>& expr) throw(VertexAlreadyOnStack)
{
#if DEBUG
  cout << "insert_expr_at: " << expr->description() << endl;
#endif
  // If it's already on then throw an exception
  if (expr->dg() == this)
    throw VertexAlreadyOnStack(expr);

  typedef RecurrenceRelation::ExprType ExprType;
  SafePtr<DGVertex> expr_vertex = static_pointer_cast<DGVertex,ExprType>(expr);

  bool new_left = true;
  bool new_right = true;
  bool need_to_clone = false;
  SafePtr<DGVertex> left_oper = expr->left();
  SafePtr<DGVertex> right_oper = expr->right();

  // See if left operand is also an operator
  SafePtr<ExprType> left_cast = dynamic_pointer_cast<ExprType,DGVertex>(left_oper);
  // if yes -- add it to the graph recursively
  if (left_cast) {
    try { insert_expr_at(expr_vertex,left_cast); }
#if ONLY_CLONE_IF_DIFF
    catch (VertexAlreadyOnStack& e) {
      if (left_oper != e.vertex()) {
        need_to_clone = true;
        left_oper = e.vertex();
      }
      new_left = false;
    }
#else
    catch (VertexAlreadyOnStack& e) { left_oper = e.vertex(); new_left = false; need_to_clone = true;}
#endif 
  }
  // else add it directly
  else {
    try { append_vertex(left_oper); }
#if ONLY_CLONE_IF_DIFF
    catch (VertexAlreadyOnStack& e) {
      if (left_oper != e.vertex()) {
        need_to_clone = true;
        left_oper = e.vertex();
      }
      new_left = false;
    }
#else
    catch (VertexAlreadyOnStack& e) { left_oper = e.vertex(); new_left = false; need_to_clone = true; }
#endif
  }

  // See if right operand is also an operator
  SafePtr<ExprType> right_cast = dynamic_pointer_cast<ExprType,DGVertex>(right_oper);
  // if yes -- add it to the graph recursively
  if (right_cast) {
    try { insert_expr_at(expr_vertex,right_cast); }
#if ONLY_CLONE_IF_DIFF
    catch (VertexAlreadyOnStack& e) {
      if (right_oper != e.vertex()) {
        need_to_clone = true;
        right_oper = e.vertex();
      }
      new_right = false;
    }
#else
    catch (VertexAlreadyOnStack& e) { right_oper = e.vertex(); new_right = false; need_to_clone = true; }
#endif
  }
  // else add it directly
  else {
    try { append_vertex(right_oper); }
#if ONLY_CLONE_IF_DIFF
    catch (VertexAlreadyOnStack& e) {
      if (right_oper != e.vertex()) {
        need_to_clone = true;
        right_oper = e.vertex();
      }
      new_right = false;
    }
#else
    catch (VertexAlreadyOnStack& e) { right_oper = e.vertex(); new_right = false; need_to_clone = true; }
#endif 
  }

  if (need_to_clone) {
    SafePtr<ExprType> expr_new(new ExprType(expr,left_oper,right_oper));
    expr_vertex = static_pointer_cast<DGVertex,ExprType>(expr_new);
    int nchildren = expr->num_exit_arcs();
#if DEBUG
    cout << "Cloned AlgebraicOperator with " << expr->num_exit_arcs() << " children" << endl;
    if (nchildren) {
      cout << "Left:  " << expr->left()->description() << endl;
      cout << "Right: " << expr->right()->description() << endl;
    }
#endif
  }

  bool new_vertex = true;
  if (new_left || new_right)
    add_new_vertex(expr_vertex);
  else {
    const bool do_cse = registry()->do_cse();
    try {
      if (do_cse) {
        append_vertex(expr_vertex);
      }
      else {
        add_new_vertex(expr_vertex);
      }
    }
    catch (VertexAlreadyOnStack& e) {
      if (new_left || new_right) {
        cout << "Problem detected: AlgebraicOperator is found on the stack but one of its operands was new" << endl;
        cout << expr_vertex->description() << endl;
        cout << e.vertex()->description() << endl;
        throw std::runtime_error("DirectedGraph::insert_expr_at() -- vertex is not new but one of the operands is");
      }
      expr_vertex = e.vertex(); new_vertex = false;
      throw;
    }
  }
  SafePtr<DGArc> left_arc(new DGArcDirect(expr_vertex,left_oper));
  expr_vertex->add_exit_arc(left_arc);
  SafePtr<DGArc> right_arc(new DGArcDirect(expr_vertex,right_oper));
  expr_vertex->add_exit_arc(right_arc);
  //SafePtr<DGArc> arc(new DGArcDirect(where,expr_vertex));
  //where->add_exit_arc(arc);
#if DEBUG
  cout << "insert_expr_at: added arc between " << where << " and " << expr_vertex << endl;
#endif

  // throw if the added expression was a cloned copy of the original
  if (need_to_clone)
    throw VertexAlreadyOnStack(expr_vertex);
}

// Replace recurrence relations with expressions
void
DirectedGraph::remove_trivial_arithmetics()
{
  const unsigned int nvertices = first_free_;
  for(int v=0; v<nvertices; v++) {

    SafePtr<DGVertex> vertex = stack_[v];
    SafePtr< AlgebraicOperator<DGVertex> > oper_cast = dynamic_pointer_cast<AlgebraicOperator<DGVertex>,DGVertex>(vertex);
    if (oper_cast) {

      SafePtr<DGVertex> left = oper_cast->exit_arc(0)->dest();
      SafePtr<DGVertex> right = oper_cast->exit_arc(1)->dest();

      // 1.0 * x = x
      if (left->equiv(prefactors.N_i[1])) {
        remove_vertex_at(vertex,right);
#if DEBUG
        cout << "Removed vertex " << vertex->description() << endl;
#endif
      }

      // x * 1.0 = x
      if (right->equiv(prefactors.N_i[1])) {
        remove_vertex_at(vertex,left);
#if DEBUG
        cout << "Removed vertex " << vertex->description() << endl;
#endif
      }

      // NOTE : more cases to come
    }
  }
}

//
// Handles "trivial" nodes. A node is trivial is it satisfies the following conditions:
// 1) it has only one child
// 2) the exit arc is of a trivial type (DGArvDirect or IntegralSet_to_Integral applied to node of size 1)
//
// By "handling" I mean either removing the node from the graph or making a node refer to another node so that
// no code is generated for it.
//
void
DirectedGraph::handle_trivial_nodes()
{
  const unsigned int nvertices = first_free_;
  for(int v=0; v<nvertices; v++) {

    SafePtr<DGVertex> vertex = stack_[v];
    if (vertex->num_exit_arcs() != 1)
      continue;
    SafePtr<DGArc> arc = vertex->exit_arc(0);

    // Is the exit arc DGArcDirect?
    {
      SafePtr<DGArcDirect> arc_cast = dynamic_pointer_cast<DGArcDirect,DGArc>(arc);
      if (arc_cast) {
        // remove the vertex, if possible
        try { remove_vertex_at(vertex,arc->dest()); }
        catch (CannotPerformOperation& c) {
        }
      }
    }

    // Is the exit arc DGArcRel<IntegralSet_to_Integrals> and vertex->size() == 1?
    {
      if (vertex->size() == 1) {
        SafePtr<DGArcRR> arc_cast = dynamic_pointer_cast<DGArcRR,DGArc>(arc);
        if (arc_cast) {
          SafePtr<RecurrenceRelation> rr = arc_cast->rr();
          SafePtr<IntegralSet_to_Integrals_base> rr_cast = dynamic_pointer_cast<IntegralSet_to_Integrals_base,RecurrenceRelation>(rr);
          if (rr_cast)
            vertex->refer_this_to(arc->dest());
        }
      }
    }
    
      // NOTE : more cases to come
  }
}


// If v1 and v2 are connected by DGArcDirect and all entry arcs to v1 are of the DGArcDirect type as well,
// this function will reattach all arcs extering v1 to v2 and remove v1 from the graph alltogether.
void
DirectedGraph::remove_vertex_at(const SafePtr<DGVertex>& v1, const SafePtr<DGVertex>& v2) throw(CannotPerformOperation)
{
  typedef vector<SafePtr<DGArc> > arcvec;
  typedef arcvec::iterator aiter;
  // Collect all entry arcs in a container
  arcvec v1_entry;
  // Verify that all entry arcs are DGArcDirect
  for(int i=0; i<v1->num_entry_arcs(); i++) {
    // See if this is a direct arc -- otherwise cannot do this
    SafePtr<DGArc> arc = v1->entry_arc(i);
    SafePtr<DGArcDirect> arc_cast = dynamic_pointer_cast<DGArcDirect,DGArc>(arc);
    if (arc_cast == 0)
      throw CannotPerformOperation("DirectedGraph::remove_vertex_at() -- cannot remove vertex");
    v1_entry.push_back(v1->entry_arc(i));
  }

  // Verify that v1 and v2 are connected by an arc and it is the only arc exiting v1
  /*if (v1->num_exit_arcs() != 1 || v1->exit_arc(0)->dest() != v2)
    throw CannotPerformOperation("DirectedGraph::remove_vertex_at() -- cannot remove vertex");*/

  // See if this is a direct arc -- otherwise cannot do this
  SafePtr<DGArc> arc = v1->exit_arc(0);
  SafePtr<DGArcDirect> arc_cast = dynamic_pointer_cast<DGArcDirect,DGArc>(arc);
  if (arc_cast == 0)
    throw CannotPerformOperation("DirectedGraph::remove_vertex_at() -- cannot remove vertex");

  //
  // OK, now do work!
  //

  // Reconnect each of v1's entry arcs to v2
  for(arcvec::iterator i=v1_entry.begin(); i != v1_entry.end(); i++) {
    SafePtr<DGVertex> parent = (*i)->orig();
    SafePtr<DGArcDirect> new_arc(new DGArcDirect(parent,v2));
    parent->replace_exit_arc(*i,new_arc);
#if DEBUG
    cout << "Replaced arcs: parent " << parent->description() << " now connected to " << new_arc->dest()->description() << endl;
    cout << "                ptr = " << parent << endl;
    cout << "               parent has " << parent->num_exit_arcs() << " children" << endl;
    cout << "               child 0 " << parent->exit_arc(0)->dest()->description() << endl;
    cout << "               child 1 " << parent->exit_arc(1)->dest()->description() << endl;
#endif
  }

  // and fully disconnect this vertex
  v1->detach();
}

void
DirectedGraph::remove_disconnected_vertices()
{
  for(int i=0; i<first_free_; i++) {
    SafePtr<DGVertex> vertex = stack_[i];
    if (vertex->num_entry_arcs() == 0 && vertex->num_exit_arcs() == 0) {
#if DEBUG
      cout << "Trying to erase disconnected vertex " << vertex->description() << " first_free_ = " << first_free_ << endl;
#endif
      try { del_vertex(vertex); }
      catch (CannotPerformOperation) {
#if DEBUG
        cout << "But couldn't!!!" << endl;
#endif
        continue;
      }
      // first_free_ was decremented, so need to decrease i as well
      --i;
    }
  }
}


//
//
//
void
DirectedGraph::generate_code(const SafePtr<CodeContext>& context, const SafePtr<MemoryManager>& memman,
                             const SafePtr<ImplicitDimensions>& dims, const SafePtr<CodeSymbols>& args,
                             const std::string& label,
                             std::ostream& decl, std::ostream& def)
{
  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  const std::string tlabel = taskmgr.current().label();

  decl << context->std_header();
  std::string comment("This code computes "); comment += label; comment += "\n";
  if (context->comments_on())
    decl << context->comment(comment) << endl;

  std::string function_name = label_to_funcname(label);
  function_name = context->label_to_name(function_name);

  decl << context->code_prefix();
  std::string func_decl;
  ostringstream oss;
  oss << context->type_name<void>() << " "
      << function_name << "(" << context->inteval_type_name(tlabel) << "* inteval";
  unsigned int nargs = args->n();
  for(unsigned int a=0; a<nargs; a++) {
    oss << ", " << context->type_name<double*>() << " "
        << args->symbol(a);
  }
  if (!dims->high_is_static()) {
    oss << ", " << context->type_name<int>() << " "
        << dims->high()->id();
  }
  if (!dims->low_is_static()) {
    oss << ", " << context->type_name<int>() << " "
        <<dims->low()->id();
  }
  oss << ")";
  func_decl = oss.str();
  
  decl << func_decl << context->end_of_stat() << endl;
  decl << context->code_postfix();

  //
  // Generate function's definition
  //

  // include standard headers
  def << context->std_header();
  // include declarations for all function calls:
  // 1) update func_names_
  // 2) include their headers into the current definition file
  update_func_names();
  for(FuncNameContainer::const_iterator fn=func_names_.begin(); fn!=func_names_.end(); fn++) {
    string function_name = (*fn).first;
    def << "#include <"
        << context->label_to_name(context->cparams()->api_prefix() + function_name)
        << ".h>" << endl;
  }
  def << endl;
  
  def << context->code_prefix();
  def << func_decl << context->open_block() << endl;
  def << context->std_function_header();
  
  context->reset();
  // if we vectorize by-line then all data is allocated on Libint's stack
  if (context->cparams()->vectorize_by_line())
    allocate_mem(memman,dims,0);
  // otherwise only arrays go on Libint's stack (scalars are handled by the compiler)
  else
    allocate_mem(memman,dims,1);
  assign_symbols(context,dims);
  print_def(context,def,dims);
  def << context->close_block() << endl;
  def << context->code_postfix();
}

void
DirectedGraph::allocate_mem(const SafePtr<MemoryManager>& memman,
                            const SafePtr<ImplicitDimensions>& dims,
                            unsigned int min_size_to_alloc)
{
  // NOTE does this belong here?
  // First, reset tag counters
  for(int i=0; i<first_free_; i++)
    stack_[i]->prepare_to_traverse();

  //
  // If need to accumulate targets, special events must happen here.
  //
  // NOTES on how to handle accumulation
  // 1) if all targets are unrolled then need to identify which integrals are part of target sets and use += instead of =
  // 2) if no targets are unrolled then allocate extra space for the target quartets and, after code has been generated,
  //    accumulate target sets into those
  // 3) if some targets are not unrolled then still need the extra space. The targets which were not unrolled should be handled
  //    as usual, i.e. not accumulated -- accumulation happens at the end
  if (registry()->accumulate_targets()) {
    // need extra buffers for targets if some are unrolled
    const bool need_copies_of_targets = nonunrolled_targets(targets_);
    iregistry()->accumulate_targets_directly(!need_copies_of_targets);

    if (need_copies_of_targets) {

      // compute the aggregate size of all targets and allocate all space at once
      ver_citer end = targets_.end();
      size size_of_targets(0);
      for(ver_iter t=targets_.begin(); t!=end; ++t) {
	size_of_targets += (*t)->size();
      }
      const address targets_buffer = memman->alloc(size_of_targets);
      iregistry()->size_of_target_accum(size_of_targets);

      // make sure that the space is at the beginning of stack
      if (targets_buffer != 0)
	throw ProgrammingError("DirectedGraph::allocate_mem() -- buffer for targets accumulation is not at the beginning of stack, need MemoryManager::alloc_at");

      // allocate every target accumulator manually
      address curr_ptr(0);
      for(ver_iter t=targets_.begin(); t!=end; ++t) {
	target_accums_.push_back(curr_ptr);
	curr_ptr += (*t)->size();
      }

    } // need copies of targets
  } // need to accumulate targets

  // Second, MUST allocate space for all targets whose symbols are not set explicitly
  // If a symbol is set means the object is not on stack (e.g. if location of target
  // is passed as an argument to set-level function)
  // This code ensures that target quartets are persistent, i.e. never overwritten, and can be accumulated into
  for(int i=0; i<first_free_; i++) {
    SafePtr<DGVertex> vertex = stack_[i];
    if (vertex->is_a_target() && !vertex->symbol_set()) {
      vertex->set_address(memman->alloc(vertex->size()));
    }
  }
  
  //
  // Go through the traversal order and at each step tag every child
  // Once a child receives same number of tags as the number of parents,
  // it can be deallocated
  //
  SafePtr<DGVertex> vertex = first_to_compute_;
  do {
    // If symbol is set then the object is not on stack
    if (!vertex->symbol_set() &&
        // address may already by set
        !vertex->address_set() &&
        // precomputed objects don't go on stack
        !vertex->precomputed() &&
        // if don't need to compute ..
        vertex->need_to_compute() &&
        // don't put on stack if smaller than min_size_to_alloc, unless it's a target
        (vertex->is_a_target() || vertex->size() > min_size_to_alloc)) {
      MemoryManager::Address addr = memman->alloc(vertex->size());
      vertex->set_address(addr);
      const unsigned int nchildren = vertex->num_exit_arcs();
      for(int c=0; c<nchildren; c++) {
        SafePtr<DGVertex> child = vertex->exit_arc(c)->dest();
        const unsigned int ntags = child->tag();
        // Do NOT deallocate if it's a target!
        if (ntags == child->num_entry_arcs() && child->address_set() && !child->is_a_target()) {
          memman->free(child->address());
        }
      }
    }
    vertex = vertex->postcalc();
  } while (vertex != 0);
}

namespace {
  std::string stack_symbol(const SafePtr<CodeContext>& ctext, const DGVertex::Address& address, const DGVertex::Size& size,
                           const std::string& low_rank, const std::string& veclen,
                           const std::string& prefix)
  {
    ostringstream oss;
    oss << prefix << "[((hsi*" << size << "+"
        << ctext->stack_address(address) << ")*" << low_rank << "+lsi)*"
        << veclen << "]";
    return oss.str();
  }
  
  /// Returns a "vector" form of stack symbol, e.g. converts libint->stack[x] to libint->stack[x+vi]
  inline std::string to_vector_symbol(const SafePtr<DGVertex>& v)
  {
    int current_pos = 0;
    std::string symb = v->symbol();
    // replace repeatedly until the string is exhausted
    while(current_pos != std::string::npos) {

      // find "[" first
      const std::string left_braket("[");
      int where = symb.find(left_braket,current_pos);
      current_pos = where;
      // if the prefix indicating a stack symbol found:
      // 1) make sure vi doesn't appear between the brakets
      // 2) replace "]" with "+vi]"
      if (where != std::string::npos) {
        const std::string right_braket("]");
        int where = symb.find(right_braket,current_pos);
        if (where == std::string::npos)
          throw logic_error("to_vector_symbol() -- address is set but no right braket found");
        
        const std::string forbidden("vi");
        int pos = symb.find(forbidden,current_pos);
        if (pos == std::string::npos || pos > where) {
          const std::string what_to_add("+vi");
          symb.insert(where,what_to_add);
          current_pos = where + 4;
        }
        else {
          current_pos = where + 1;
        }
      }
    } // end of while
    return symb;
  }
};

void
DirectedGraph::assign_symbols(const SafePtr<CodeContext>& context, const SafePtr<ImplicitDimensions>& dims)
{
  std::ostringstream os;
  const std::string null_str("");
  const std::string& stack_name = registry()->stack_name();
  
  // Generate the label for the rank of the low dimension
  std::string low_rank = dims->low_label();
  std::string veclen = dims->vecdim_label();

  // First, set symbols for all vertices which have address assigned
  for(int i=0; i<first_free_; i++) {
    SafePtr<DGVertex> vertex = stack_[i];
    if (!vertex->symbol_set() && vertex->address_set()) {
      vertex->set_symbol(stack_symbol(context,vertex->address(),vertex->size(),low_rank,veclen,stack_name));
    }
  }

  // Second, find all nodes which were unrolled using IntegralSet_to_Integrals:
  // 1) such nodes do not need symbols generated since they never appear in the code expicitly
  // 2) children of such nodes have symbols that depend on the parent's address
  for(int i=0; i<first_free_; i++) {
    SafePtr<DGVertex> vertex = stack_[i];
    if (vertex->num_exit_arcs() == 0)
      continue;
    SafePtr<DGArc> arc = vertex->exit_arc(0);
    SafePtr<DGArcRR> arc_rr = dynamic_pointer_cast<DGArcRR,DGArc>(arc);
    if (arc_rr == 0)
      continue;
    SafePtr<RecurrenceRelation> rr = arc_rr->rr();
    SafePtr<IntegralSet_to_Integrals_base> iset_to_i = dynamic_pointer_cast<IntegralSet_to_Integrals_base,RecurrenceRelation>(rr);
    if (iset_to_i == 0) {
      continue;
    }
    else {
      unsigned int nchildren = vertex->num_exit_arcs();
      for(int c=0; c<nchildren; c++) {
        SafePtr<DGVertex> child = vertex->exit_arc(c)->dest();
        // If a child is precomputed and it's parent symbol is not set -- its symbol will be set as usual
        if (!child->precomputed() || vertex->symbol_set()) {
          if (vertex->address_set()) {
            child->set_symbol(stack_symbol(context,vertex->address()+c,vertex->size(),low_rank,veclen,stack_name));
          }
          else {
            child->set_symbol(stack_symbol(context,c,vertex->size(),low_rank,veclen,vertex->symbol()));
          }
        }
      }
      vertex->refer_this_to(vertex->exit_arc(0)->dest());
      vertex->reset_symbol();
    }
  }

  // then process all other symbols, EXCEPT operators
  for(int i=0; i<first_free_; i++) {

    SafePtr<DGVertex> vertex = stack_[i];
#if DEBUG
    cout << "Trying to assign symbol to " << vertex->description() << endl;
#endif
    if (vertex->symbol_set()) {
      continue;
    }

#if 0
    // test if the vertex is an operator
    {
      typedef AlgebraicOperator<DGVertex> oper;
      SafePtr<oper> ptr_cast = dynamic_pointer_cast<oper,DGVertex>(vertex);
      if (ptr_cast) {
        vertex->set_symbol(context->unique_name<EntityTypes::FP>());
        continue;
      }
    }
#endif

    // test if the vertex is a static quantity, like a constant
    {
      typedef CTimeEntity<double> cdouble;
      SafePtr<cdouble> ptr_cast = dynamic_pointer_cast<cdouble,DGVertex>(vertex);
      if (ptr_cast) {
        vertex->set_symbol(ptr_cast->label());
        continue;
      }
    }

    // test if the vertex is precomputed runtime quantity, like a geometric parameter
    if (vertex->precomputed()) {
      std::string symbol("inteval->");
      symbol += context->label_to_name(vertex->label());
      symbol += "[vi]";
      vertex->set_symbol(symbol);
      continue;
    }

    // test if the vertex is other runtime quantity
    {
      typedef RTimeEntity<double> cdouble;
      SafePtr<cdouble> ptr_cast = dynamic_pointer_cast<cdouble,DGVertex>(vertex);
      if (ptr_cast) {
        vertex->set_symbol(ptr_cast->label());
        continue;
      }
    }

  } // done with everything BUT operators
  
  // finally, process all operators (start with most recently added vertices since those are
  // much more likely to be on the bottom of the graph).
  for(int i=first_free_-1; i>=0; --i) {

    SafePtr<DGVertex> vertex = stack_[i];
#if DEBUG
    cout << "Trying to assign symbol to " << vertex->description() << endl;
#endif
    assign_oper_symbol(context,vertex);
  }

}

void
DirectedGraph::assign_oper_symbol(const SafePtr<CodeContext>& context, SafePtr<DGVertex>& vertex)
{
  // do nothing if the vertex has a symbol or is not an operator
  if (vertex->symbol_set())
    return;
  
  {
    typedef AlgebraicOperator<DGVertex> oper;
    SafePtr<oper> ptr_cast = dynamic_pointer_cast<oper,DGVertex>(vertex);
    if (ptr_cast) {
      // is it in a subtree?
      const bool on_a_subtree = (vertex->subtree());
      
      // If no -- it will be an automatic variable
      if (!on_a_subtree)
        vertex->set_symbol(context->unique_name<EntityTypes::FP>());
      // else assign symbols to left and right arguments
      else {
        SafePtr<DGVertex> left = ptr_cast->exit_arc(0)->dest();
        SafePtr<DGVertex> right = ptr_cast->exit_arc(1)->dest();
        assign_oper_symbol(context,left);
        assign_oper_symbol(context,right);
        
        std::ostringstream oss;
        oss << "( " << left->symbol() << " ) "
            << ptr_cast->label()
            << " ( " << right->symbol() << " )";
        vertex->set_symbol(oss.str());
      }
    }
  }
}


namespace {
  /// returns how many times token appears in str
  unsigned int nfind(const std::string& str, const std::string& token)
  {
    unsigned int nfinds = 0;
    typedef std::string::size_type size_type;
    size_type current_pos = 0;
    while(1) {
      current_pos = str.find(token,current_pos);
      if (current_pos != std::string::npos) {
        ++nfinds;
        ++current_pos;
      }
      else
        return nfinds;
    }
  }
  
#define DO_NOT_COUNT_DIV 1
  /// Returns the number of FLOPs in an expression
  unsigned int nflops(const std::string& expr)
  {
    static const std::string mul(" * ");
    static const std::string div(" / ");
    static const std::string plus(" + ");
    static const std::string minus(" - ");
    const unsigned int nflops =
      nfind(expr,mul) +
      nfind(expr,plus) +
      nfind(expr,minus)
#if !DO_NOT_COUNT_DIV
      + nfind(expr,div)
#endif
      ;
    return nflops;
  }
}

namespace {
  SafePtr<MemBlockSet>
  to_memoryblks(DirectedGraph::vertices& vertices) {
    SafePtr<MemBlockSet> result(new MemBlockSet);
    typedef DirectedGraph::ver_citer citer;
    typedef DirectedGraph::ver_iter iter;
    citer end(vertices.end());
    for(iter v=vertices.begin(); v!=end; ++v) {
      result->push_back(MemBlock((*v)->address(),(*v)->size(),false,SafePtr<MemBlock>(),SafePtr<MemBlock>()));
    }
    return result;
  }
};

void
DirectedGraph::print_def(const SafePtr<CodeContext>& context, std::ostream& os,
                        const SafePtr<ImplicitDimensions>& dims)
{
  std::ostringstream oss;
  const std::string null_str("");
  SafePtr<Entity > ctimeconst_zero(new CTimeEntity<int>(0));

  const bool accumulate_targets_directly = registry()->accumulate_targets() && iregistry()->accumulate_targets_directly();
  const bool accumulate_targets_indirectly = registry()->accumulate_targets() && !accumulate_targets_directly;

  //
  // To optimize accumulation and setting blocks to zero, need to merge maximally the targets (the accumulation area is one block)
  //
  SafePtr<MemBlockSet> targetblks;
  if (registry()->accumulate_targets()) {
    if (accumulate_targets_directly) {
      targetblks = to_memoryblks(targets_);
      merge(*targetblks);
    }
    else {
      targetblks = SafePtr<MemBlockSet>(new MemBlockSet);
      targetblks->push_back(MemBlock(0,iregistry()->size_of_target_accum(),false,SafePtr<MemBlock>(),SafePtr<MemBlock>()));
    }
  }

  //
  // If accumulating integrals, check inteval's zero_out_targets. If set to 1 -- zero out accumulated targets
  //
  if (registry()->accumulate_targets()) {

    const bool vecdim_is_static = dims->vecdim_is_static();

    os << "if (inteval->zero_out_targets) {" << std::endl;

    typedef MemBlockSet::const_iterator citer;
    typedef MemBlockSet::iterator iter;
    const citer end = targetblks->end();
    for(iter b=targetblks->begin(); b!=end; ++b) {

      size s = b->size();
      SafePtr<CTimeEntity<int> > bdim(new CTimeEntity<int>(s));

      SafePtr<Entity> bvecdim;
      if (!vecdim_is_static) {
	SafePtr< RTimeEntity<EntityTypes::Int> > vecdim = dynamic_pointer_cast<RTimeEntity<EntityTypes::Int>,Entity>(dims->vecdim());
	bvecdim = vecdim * bdim;
      }
      else {
	SafePtr< CTimeEntity<int> > vecdim = dynamic_pointer_cast<CTimeEntity<int>,Entity>(dims->vecdim());
	bvecdim = vecdim * bdim;
      }

#if 0
      std::string loopvar("i");
      ForLoop loop(context,loopvar,bvecdim,ctimeconst_zero);
      os << loop.open();
      {
	ostringstream oss;
	oss << registry()->stack_name() << "[" << loopvar << "]";
	const std::string zero("0");
	os << context->assign(oss.str(),zero);
      }
      os << loop.close();
#endif
      os << "_libint2_static_api_bzero_short_(" << registry()->stack_name() << "+"
	 << b->address() << "*" << dims->vecdim()->id() << "," << bvecdim->id() << ")" << endl;

    }

    os << "inteval->zero_out_targets = 0;" << std::endl << "}" << std::endl;
  }

  std::string varname("hsi");
  SafePtr<ForLoop> hsi_loop(new ForLoop(context,varname,dims->high(),SafePtr<Entity>(new CTimeEntity<int>(0))));
  os << hsi_loop->open();
  
  varname = "lsi";
  SafePtr<ForLoop> lsi_loop(new ForLoop(context,varname,dims->low(),SafePtr<Entity>(new CTimeEntity<int>(0))));
  os << lsi_loop->open();
  
  // the vector loop is created outside of the body of the function if
  // 1) blockwise vectorization is requested
  // and
  // 2) this is a purely int-unit code, i.e. there are no explicit RRs on sets in the body
  // Otherwise, I create a dummy vector loop with the vector loop index set to 0
  const unsigned int max_vector_length = context->cparams()->max_vector_length();
  const bool vectorize = (max_vector_length != 1);
  const bool vectorize_by_line = context->cparams()->vectorize_by_line();
  const bool create_outer_vector_loop = !vectorize_by_line && !contains_nontrivial_rr();
  varname = "vi";
  // outer vector loop
  SafePtr<ForLoop> outer_vloop;
  // vector loop for each code line
  SafePtr<ForLoop> line_vloop;
  if (create_outer_vector_loop) {
    SafePtr<ForLoop> tmp_vi_loop(new ForLoop(context,varname,dims->vecdim(),SafePtr<Entity>(new CTimeEntity<int>(0))));
    outer_vloop = tmp_vi_loop;
  }
  else {
    SafePtr<Entity> unit_dim(new CTimeEntity<int>(1));
    SafePtr<ForLoop> tmp_vi_loop(new ForLoop(context,varname,unit_dim,SafePtr<Entity>(new CTimeEntity<int>(0))));
    outer_vloop = tmp_vi_loop;
    SafePtr<ForLoop> tmp2_vi_loop(new ForLoop(context,varname,dims->vecdim(),SafePtr<Entity>(new CTimeEntity<int>(0))));
    line_vloop = tmp2_vi_loop;
    // note that both loops use same variable name -- standard C++ scoping rules allow it -- but the outer
    // loop will become a declaration of a constant variable
  }
  os << outer_vloop->open();

  //
  // generate code for vertices
  //
  unsigned int nflops_total = 0;
  SafePtr<DGVertex> current_vertex = first_to_compute_;
  do {

    // for every vertex that has a defined symbol, hence must be defined in code
    if (current_vertex->symbol_set()) {

      const bool address_set = current_vertex->address_set();

      // print algebraic expression
      {
        typedef AlgebraicOperator<DGVertex> oper_type;
        SafePtr<oper_type> oper_ptr = dynamic_pointer_cast<oper_type,DGVertex>(current_vertex);
        if (oper_ptr) {

          // Type declaration if this is an automatic variable (i.e. not on Libint's stack)
          if (!address_set)
            os << context->declare(context->type_name<double>(),
                                   current_vertex->symbol());
          
	  // If this is an Integral in a target IntegralSet AND
	  // can accumulate targets directly -- use '+=' instead of '='
	  const bool accumulate_not_assign = accumulate_targets_directly && IntegralInTargetIntegralSet()(current_vertex);

          if (context->comments_on()) {
            oss.str(null_str);
            oss << current_vertex->label() << (accumulate_not_assign ? " += " : " = ")
                << oper_ptr->exit_arc(0)->dest()->label()
                << oper_ptr->label()
                << oper_ptr->exit_arc(1)->dest()->label();
            os << context->comment(oss.str()) << endl;
          }

          // expression
          SafePtr<DGVertex> left_arg = oper_ptr->exit_arc(0)->dest();
          SafePtr<DGVertex> right_arg = oper_ptr->exit_arc(1)->dest();
#if DEBUG
          cout << "Generating code for " << current_vertex->description() << endl;
          cout << "              ptr = " << current_vertex << endl;
          cout << "         left_arg = " << left_arg->description() << endl;
          cout << "        right_arg = " << right_arg->description() << endl;
#endif
          
          // convert symbols to their vector form if needed
          std::string curr_symbol = current_vertex->symbol();
          std::string left_symbol = left_arg->symbol();
          std::string right_symbol = right_arg->symbol();
          if (vectorize) {
            curr_symbol = to_vector_symbol(current_vertex);
            left_symbol = to_vector_symbol(left_arg);
            right_symbol = to_vector_symbol(right_arg);
          }
          
          if (vectorize_by_line)
            os << line_vloop->open();
	  // the statement that does the work
	  {
	    if (accumulate_not_assign) {
	      os << context->accumulate_binary_expr(curr_symbol,left_symbol,oper_ptr->label(),right_symbol);
	      nflops_total += (1 + nflops(left_symbol) + nflops(right_symbol)) + 1;
	    }
	    else {
	      os << context->assign_binary_expr(curr_symbol,left_symbol,oper_ptr->label(),right_symbol);
	      nflops_total += (1 + nflops(left_symbol) + nflops(right_symbol));
	    }
	  }
          if (vectorize_by_line)
            os << line_vloop->close();


          goto next;
        }
      }
      
      // print simple assignment statement
      if (current_vertex->num_exit_arcs() == 1) {
        typedef DGArcDirect arc_type;
        SafePtr<arc_type> arc_ptr = dynamic_pointer_cast<arc_type,DGArc>(current_vertex->exit_arc(0));
        if (arc_ptr) {

	  // If this is an Integral in a target IntegralSet AND
	  // can accumulate targets directly -- use '+=' instead of '='
	  const bool accumulate_not_assign = accumulate_targets_directly && IntegralInTargetIntegralSet()(current_vertex);

          if (context->comments_on()) {
            oss.str(null_str);
            oss << current_vertex->label() << (accumulate_not_assign ? " += " : " = ")
                << arc_ptr->dest()->label();
            os << context->comment(oss.str()) << endl;
          }

          // convert symbols to their vector form if needed
          std::string curr_symbol = current_vertex->symbol();
          std::string rhs_symbol = arc_ptr->dest()->symbol();
          if (vectorize) {
            curr_symbol = to_vector_symbol(current_vertex);
            rhs_symbol = to_vector_symbol(arc_ptr->dest());
          }
          
          if (vectorize_by_line)
            os << line_vloop->open();
	  if (accumulate_not_assign) {
	    os << context->accumulate(curr_symbol,
				      rhs_symbol);
	    nflops_total += nflops(rhs_symbol) + 1;   // +1 due to +=
	  }
	  else {
	    os << context->assign(curr_symbol,
				  rhs_symbol);
	    nflops_total += nflops(rhs_symbol);
	  }
          if (vectorize_by_line)
            os << line_vloop->close();
          
          
          goto next;
        }
      }

      // print out a recurrence relation
      {
        typedef DGArcRR arc_type;
        SafePtr<arc_type> arc_ptr = dynamic_pointer_cast<arc_type,DGArc>(current_vertex->exit_arc(0));
        if (arc_ptr) {
          
          SafePtr<RecurrenceRelation> rr = arc_ptr->rr();
          os << rr->spfunction_call(context,dims);
          
          goto next;
        }
      }

      {
        throw std::runtime_error("DirectedGraph::print_def() -- cannot handle this vertex yet");
      }
    }

    next:
    current_vertex = current_vertex->postcalc();
  } while (current_vertex != 0);

  os << outer_vloop->close();
  os << lsi_loop->close();
  os << hsi_loop->close();

  //
  // Accumulate targets
  //
  if (accumulate_targets_indirectly) {
    os << context->comment("Accumulate target integral sets") << std::endl;
    const std::string& stack_name = registry()->stack_name();
    const bool vecdim_is_static = dims->vecdim_is_static();
#if 0
    unsigned int vecdim_rank;
    if (vecdim_is_static) {
      SafePtr<CTimeEntity<int> > cptr = dynamic_pointer_cast<CTimeEntity<int> ,Entity >(dims->vecdim());
      vecdim_rank = cptr->value();
    }
    const std::string times_vecdim("*" + dims->vecdim_label());
#endif

    // Loop over targets
    unsigned int curr_target = 0;
    const ver_iter end = targets_.end();
    for(ver_iter t=targets_.begin(); t!=targets_.end(); ++t, ++curr_target) {

      size s = (*t)->size();
      SafePtr<CTimeEntity<int> > bdim(new CTimeEntity<int>(s));

      SafePtr<Entity> bvecdim;
      if (!vecdim_is_static) {
	SafePtr< RTimeEntity<EntityTypes::Int> > vecdim = dynamic_pointer_cast<RTimeEntity<EntityTypes::Int>,Entity>(dims->vecdim());
	bvecdim = vecdim * bdim;
      }
      else {
	SafePtr< CTimeEntity<int> > vecdim = dynamic_pointer_cast<CTimeEntity<int>,Entity>(dims->vecdim());
	bvecdim = vecdim * bdim;
      }

      // For now write an explicit loop for each target. In the future should:
      // 1) check if all computed targets are adjacent (accumulated targets are adjacent by the virtue of the allocation mechanism)
      // 2) check the sizes and insert optimized calls, if possible

      // form a single loop over integrals and vector dimension
      // NOTE: single loop suffices because if outer/inner strides are not 1 this block of code should not be executed

#if 0
      SafePtr<Entity> loopmax;
      if (vecdim_is_static) {
	const int loopmax_value = s*vecdim_rank;
	loopmax = SafePtr<Entity>(new CTimeEntity<int>(loopmax_value));
      }
      else {
	ostringstream oss;  oss << s << times_vecdim;
	loopmax = SafePtr<Entity>(new RTimeEntity<EntityTypes::Int>(oss.str()));
      }

      std::string loopvar("i");
      ForLoop loop(context,loopvar,loopmax,ctimeconst_zero);
      os << loop.open();
      std::string acctarget;
      {
	ostringstream oss;
	oss << stack_name << "[" << target_accums_[curr_target] << times_vecdim << "+" << loopvar << "]";
	acctarget = oss.str();
      }
      std::string target;
      {
	ostringstream oss;
	oss << stack_name << "[" << (*t)->address() << times_vecdim << "+" << loopvar << "]";
	target = oss.str();
      }
      os << context->accumulate(acctarget,target);
      os << loop.close();
#endif
      os << "_libint2_static_api_inc_short_("
	 << registry()->stack_name() << "+" << target_accums_[curr_target] << "*" << dims->vecdim()->id() << ","
	 << registry()->stack_name() << "+" << (*t)->address() << "*" << dims->vecdim()->id() << ","
	 << bvecdim->id() << ","
	 << "1.0)" << endl;

      nflops_total += s;
    }
  }

  // Outside of loops stack symbols don't make sense, so we must define loop variables hsi, lsi, and vi to 0
  os << context->decldef(context->type_name<const int>(), "hsi", "0");
  os << context->decldef(context->type_name<const int>(), "lsi", "0");
  os << context->decldef(context->type_name<const int>(), "vi", "0");
  
  //
  // Now pass back all targets through the inteval object, if needed.
  //
  if (registry()->return_targets()) {
    unsigned int curr_target = 0;
    const ver_iter end = targets_.end();
    for(ver_iter t=targets_.begin(); t!=targets_.end(); ++t, ++curr_target) {
      const std::string& symbol = (accumulate_targets_indirectly
				   //                                                                    is this correct?         ???
				   ? stack_symbol(context,target_accums_[curr_target],(*t)->size(),dims->low_label(),dims->vecdim_label(),registry()->stack_name())
				   : (*t)->symbol());
      os << "inteval->targets[" << curr_target << "] = "
	 << context->value_to_pointer(symbol) << context->end_of_stat() << endl;
    }
  }
  
  // Print out the number of flops
  oss.str(null_str);
  oss << "Number of flops = " << nflops_total;
  os << context->comment(oss.str()) << endl;
  
  if (context->cparams()->count_flops()) {
    oss.str(null_str);
    oss << nflops_total << " * " << dims->high_label() << " * "
        << dims->low_label() << " * "
        << dims->vecdim_label();
    os << context->assign_binary_expr("inteval->nflops","inteval->nflops","+",oss.str());
  }

}

void
DirectedGraph::update_func_names()
{
  // Loop over all vertices
  for(int i=0; i<first_free_; i++) {
    ver_ptr v = stack_[i];
    // for every vertex with children
    if (v->num_exit_arcs() > 0) {
      // if it must be computed using a RR
      SafePtr<DGArc> arc = v->exit_arc(0);
      SafePtr<DGArcRR> arcrr = dynamic_pointer_cast<DGArcRR,DGArc>(arc);
      if (arcrr != 0) {
        SafePtr<RecurrenceRelation> rr = arcrr->rr();
        // and the RR is complex (i.e. likely to result in a function call)
        if (!rr->is_simple())
          // add it to the RRStack
          func_names_[rr->label()] = true;
      }
    }
  }
}

bool
DirectedGraph::contains_nontrivial_rr() const
{
  SafePtr<DGVertex> current_vertex = first_to_compute_;
  do {
    const int nchildren = current_vertex->num_exit_arcs();
    if (nchildren > 0) {
      arc_ptr aptr = current_vertex->exit_arc(0);
      typedef RecurrenceRelation RR;
      SafePtr<DGArcRR> aptr_cast = dynamic_pointer_cast<DGArcRR,arc>(aptr);
      // if this is a RR
      if (aptr_cast != 0) {
        // and a non-trivial one
        if (!aptr_cast->rr()->is_simple())
          return true;
      }
    }
    current_vertex = current_vertex->postcalc();
  } while (current_vertex != 0);
  
  return false;
}

void
DirectedGraph::find_subtrees()
{
  // Find subtrees by starting from the targets and moving down ...
  for(int i=0; i<first_free_; i++) {
    SafePtr<DGVertex> v = stack_[i];
    
    if (v->is_a_target() && v->num_entry_arcs() == 0) {
      find_subtrees_from(v);
    }
  }
}

void
DirectedGraph::find_subtrees_from(const SafePtr<DGVertex>& v)
{
  // is not on a subtree already
  if (!v->subtree()) {
    
    bool useless_subtree = false;
    
    //
    // Subtrees are useless in the following cases:
    // 1) root is computed via RRs, not explicitly
    //
    {
      if (v->num_exit_arcs() > 0) {
        SafePtr<DGArc> arc = v->exit_arc(0);
        SafePtr<DGArcRR> arc_rr = dynamic_pointer_cast<DGArcRR,DGArc>(arc);
        if (arc_rr)
          useless_subtree = true;
      }
    }
    
    // create subtree
    if (!useless_subtree) {
      SafePtr<DRTree> stree = DRTree::CreateRootedAt(v);
      
      // Remove all trivial subtrees
      if (stree) {
#if DISABLE_SUBTREES
        if (stree->nvertices() >= 0) {
          stree->detach();
        }
#else
        if (stree->nvertices() < 3) {
          stree->detach();
        }
#endif
      }
    }
    
    // move on to children
    const unsigned int nchildren = v->num_exit_arcs();
    for(unsigned int c=0; c<nchildren; c++)
      find_subtrees_from(v->exit_arc(c)->dest());
  }
}

////

namespace libint2 {

  bool
  nonunrolled_targets(const DirectedGraph::vertices& targets) {
    typedef DirectedGraph::ver_citer citer;
    citer end = targets.end();
    if (end != find_if(targets.begin(),end,NotUnrolledIntegralSet()))
      return true;
    else
      return false;
  }

};
