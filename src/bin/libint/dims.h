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
#include <entity.h>
#include <smart_ptr.h>

#ifndef _libint2_src_bin_libint_dims_h_
#define _libint2_src_bin_libint_dims_h_

namespace libint2 {

using namespace EntityTypes;

/** ImplicitDimensions describes basis functions or other "degrees of freedom"
    not actively engaged in a recurrence relation. For example, horizontal
    AM transfer to obtain a (dd|ps) integral from (fp|ps) and (dp|ps) does
    not involve the |ps) part. Thus there's may be no reason to generate
   transfer routine specific to this integral, only a routine specific to the
   (dd| part. Such function will require the information about the rank of the
   |ps) part. This information is encoded in ImplicitDimensions.

    Another special dimension is the vector length...
*/

class ImplicitDimensions {
 public:
  /// Explicitly initialize both quantities. Their exact type is not known.
  ImplicitDimensions(const std::shared_ptr<Entity>& high,
                     const std::shared_ptr<Entity>& low,
                     const std::shared_ptr<Entity>& vecdim);
  /// Default assumes runtime (dynamical) quantities
  ImplicitDimensions();
  /// Handy constructor to initialize dimensions as compile-time (static)
  /// quatities
  ImplicitDimensions(int high, int low, int vec);
  ~ImplicitDimensions() {}

  /// Returns the high dimension
  std::shared_ptr<Entity> high() const { return high_; }
  /// Returns the low dimension
  std::shared_ptr<Entity> low() const { return low_; }
  /// Returns the vector dimension
  std::shared_ptr<Entity> vecdim() const { return vecdim_; }
  /// Returns true if the rank of high dimension is known
  bool high_is_static() const { return high_is_static_; }
  /// Returns true if the rank of low dimension is known
  bool low_is_static() const { return low_is_static_; }
  /// Returns true if the rank of vector dimension is known
  bool vecdim_is_static() const { return vecdim_is_static_; }
  /// Returns the label of the high dimension
  const std::string& high_label() const { return high_label_; }
  /// Returns the label of the low dimension
  const std::string& low_label() const { return low_label_; }
  /// Returns the label of the vector dimension
  const std::string& vecdim_label() const { return vecdim_label_; }

  /// Sets default ImplicitDimension object
  static void set_default_dims(
      const std::shared_ptr<CompilationParameters>& cparams);
  /// Default ImplicitDimension object
  static std::shared_ptr<ImplicitDimensions> default_dims();

 private:
  // Dimensions can be runtime or compile-time quantities
  const std::shared_ptr<Entity> high_;
  const std::shared_ptr<Entity> low_;
  const std::shared_ptr<Entity> vecdim_;

  // checks if the dimensions are CTImeEntities
  void init_();
  bool high_is_static_;
  bool low_is_static_;
  bool vecdim_is_static_;
  // Cached labels for
  std::string high_label_;
  std::string low_label_;
  std::string vecdim_label_;

  /// Default dimension
  static std::shared_ptr<ImplicitDimensions> default_dims_;
};

};  // namespace libint2

#endif
