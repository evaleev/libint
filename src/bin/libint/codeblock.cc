
#include <codeblock.h>
#include <context.h>

using namespace std;
using namespace libint2;

ForLoop::ForLoop(const SafePtr<CodeContext>& context, std::string& varname,
                 const SafePtr<Entity>& less_than, const SafePtr<Entity>& start_at) :
  CodeBlock(context), varname_(varname), less_than_(less_than), start_at_(start_at)
{
  init_();
}

ForLoop::~ForLoop()
{
}

void
ForLoop::init_()
{
  SafePtr<CodeContext> ctext = context();
  /// Why in the hell is dynamic_pointer_cast broken?
  SafePtr< CTimeEntity<int> > lt_cptr = boost::shared_ptr< CTimeEntity<int> >(less_than_,boost::detail::dynamic_cast_tag());
  SafePtr< CTimeEntity<int> > sa_cptr = boost::shared_ptr< CTimeEntity<int> >(start_at_,boost::detail::dynamic_cast_tag());
  SafePtr< RTimeEntity<EntityTypes::Int> > lt_rptr = boost::shared_ptr< RTimeEntity<EntityTypes::Int> >(less_than_,boost::detail::dynamic_cast_tag());
  SafePtr< RTimeEntity<EntityTypes::Int> > sa_rptr = boost::shared_ptr< RTimeEntity<EntityTypes::Int> >(start_at_,boost::detail::dynamic_cast_tag());
  
  if (lt_cptr != 0) {
    ostringstream oss;
    oss << lt_cptr->value();
    lt_expr_ = oss.str();
  }
  else if (lt_rptr != 0) {
    lt_expr_ = ctext->label_to_name(lt_rptr->label());
  }
  else
    throw std::logic_error("ForLoop::open -- less_than does not have one of desired types");
  if (sa_cptr != 0) {
    ostringstream oss;
    oss << sa_cptr->value();
    sa_expr_ = oss.str();
  }
  else if (sa_rptr != 0) {
    sa_expr_ = ctext->label_to_name(sa_rptr->label());
  }
  else
    throw std::logic_error("ForLoop::open -- less_than does not have one of desired types");
  
  if (lt_cptr !=0 && sa_cptr !=0 &&
      lt_cptr->value() == sa_cptr->value()+1 )
    dummy_loop_ = true;
  else
    dummy_loop_ = false;
}

std::string
ForLoop::open()
{
  SafePtr<CodeContext> ctext = context();
  if (dummy_loop_) {
      return ctext->decldef(ctext->type_name<const int>(), varname_, sa_expr_);
  }
  else {
    ostringstream oss;
    oss << "for(" << ctext->type_name<int>() << " " << varname_ << " = " << sa_expr_ << ctext->end_of_stat()
    << " " << varname_ << "<" << lt_expr_ << ctext->end_of_stat() << " " << varname_ << "++) {" << endl;
    return oss.str();
  }
}

std::string
ForLoop::close()
{
  ostringstream oss;
  if (dummy_loop_)
    oss << "";
  else
    oss << "}";
  oss << endl;
  return oss.str();
}

