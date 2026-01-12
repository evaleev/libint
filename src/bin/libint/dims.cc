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

#include <default_params.h>
#include <dims.h>

using namespace std;
using namespace libint2;

std::shared_ptr<ImplicitDimensions> ImplicitDimensions::default_dims_ =
    std::shared_ptr<ImplicitDimensions>();

void ImplicitDimensions::set_default_dims(
    const std::shared_ptr<CompilationParameters>& cparams) {
  std::shared_ptr<ImplicitDimensions> new_default(
      new ImplicitDimensions(1, 1, cparams->max_vector_length()));
  default_dims_ = new_default;
}

std::shared_ptr<ImplicitDimensions> ImplicitDimensions::default_dims() {
  if (default_dims_ == 0)
    throw std::logic_error(
        "ImplicitDimensions::default_dims() -- set_default_dims() has not been "
        "called yet");
  return default_dims_;
}

ImplicitDimensions::ImplicitDimensions(const std::shared_ptr<Entity>& high,
                                       const std::shared_ptr<Entity>& low,
                                       const std::shared_ptr<Entity>& vecdim)
    : high_(high), low_(low), vecdim_(vecdim) {
  init_();
}

ImplicitDimensions::ImplicitDimensions()
    : high_(std::shared_ptr<Entity>(
          new RTimeEntity<EntityTypes::Int>("highdim"))),
      low_(
          std::shared_ptr<Entity>(new RTimeEntity<EntityTypes::Int>("lowdim"))),
      vecdim_(std::shared_ptr<Entity>(
          new RTimeEntity<EntityTypes::Int>("libint->veclen"))) {
  init_();
}

ImplicitDimensions::ImplicitDimensions(int high, int low, int vec)
    : high_(std::shared_ptr<Entity>(new CTimeEntity<int>(high))),
      low_(std::shared_ptr<Entity>(new CTimeEntity<int>(low))),
      vecdim_(std::shared_ptr<Entity>(new CTimeEntity<int>(vec))) {
  init_();
}

void ImplicitDimensions::init_() {
  std::shared_ptr<CTimeEntity<int> > cptr =
      std::dynamic_pointer_cast<CTimeEntity<int>, Entity>(high_);
  if (cptr != 0) {
    high_is_static_ = true;
    ostringstream oss;
    oss << cptr->value();
    high_label_ = oss.str();
  } else {
    high_is_static_ = false;
    std::shared_ptr<DGVertex> dptr =
        std::dynamic_pointer_cast<DGVertex, Entity>(high_);
    high_label_ = dptr->label();
  }
  cptr = std::dynamic_pointer_cast<CTimeEntity<int>, Entity>(low_);
  if (cptr != 0) {
    low_is_static_ = true;
    ostringstream oss;
    oss << cptr->value();
    low_label_ = oss.str();
  } else {
    low_is_static_ = false;
    std::shared_ptr<DGVertex> dptr =
        std::dynamic_pointer_cast<DGVertex, Entity>(low_);
    low_label_ = dptr->label();
  }
  cptr = std::dynamic_pointer_cast<CTimeEntity<int>, Entity>(vecdim_);
  if (cptr != 0) {
    vecdim_is_static_ = true;
    ostringstream oss;
    oss << cptr->value();
    vecdim_label_ = oss.str();
  } else {
    vecdim_is_static_ = false;
    std::shared_ptr<DGVertex> dptr =
        std::dynamic_pointer_cast<DGVertex, Entity>(vecdim_);
    vecdim_label_ = dptr->label();
  }
}
