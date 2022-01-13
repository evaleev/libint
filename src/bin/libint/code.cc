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

#include <code.h>

using namespace std;
using namespace libint2;

CodeSymbols::CodeSymbols(): symbols_() {}

CodeSymbols::~CodeSymbols() {}

void
CodeSymbols::append_symbol(const std::string& s)
{
  symbols_.push_back(s);
}

unsigned int
CodeSymbols::n() const
{
  return symbols_.size();
}

const std::string&
CodeSymbols::symbol(unsigned int i) const
{
  return symbols_.at(i);
}

