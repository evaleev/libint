/*
 *  Copyright (C) 2004-2021 Edward F. Valeev
 *
 *  This file is part of Libint.
 *
 *  Libint is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Libint.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

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
  SafePtr< CTimeEntity<int> > lt_cptr = dynamic_pointer_cast<CTimeEntity<int>,Entity>(less_than_);
  SafePtr< CTimeEntity<int> > sa_cptr = dynamic_pointer_cast<CTimeEntity<int>,Entity>(start_at_);
  SafePtr< RTimeEntity<EntityTypes::Int> > lt_rptr = dynamic_pointer_cast<RTimeEntity<EntityTypes::Int>,Entity>(less_than_);
  SafePtr< RTimeEntity<EntityTypes::Int> > sa_rptr = dynamic_pointer_cast<RTimeEntity<EntityTypes::Int>,Entity>(start_at_);
  
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
  ostringstream oss;

  if (dummy_loop_) {
    oss << "{" << endl
        << ctext->decldef(ctext->type_name<const int>(), varname_, sa_expr_);
  }
  else {
    oss << "#ifdef __INTEL_COMPILER\n#pragma ivdep\n#endif\n"
        << "for(" << ctext->type_name<int>() << " " << varname_ << " = " << sa_expr_ << ctext->end_of_stat()
    << " " << varname_ << "<" << lt_expr_ << ctext->end_of_stat() << " " << varname_ << "++) {" << endl;
  }
  return oss.str();
}

std::string
ForLoop::close()
{
  ostringstream oss;
  oss << "}" << endl;
  return oss.str();
}

