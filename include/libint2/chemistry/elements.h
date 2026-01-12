/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint library.
 *
 *  Libint library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef INCLUDE_LIBINT2_CHEMISTRY_ELEMENTS_H_
#define INCLUDE_LIBINT2_CHEMISTRY_ELEMENTS_H_

#include <string>
#include <vector>

namespace libint2 {
namespace chemistry {

struct element {
  unsigned short Z;
  std::string name;
  std::string symbol;
};

static const std::vector<element>& get_element_info() {
  static std::vector<element> element_info{
      {1, "hydrogen", "H"},       {2, "helium", "He"},
      {3, "lithium", "Li"},       {4, "beryllium", "Be"},
      {5, "boron", "B"},          {6, "carbon", "C"},
      {7, "nitrogen", "N"},       {8, "oxygen", "O"},
      {9, "fluorine", "F"},       {10, "neon", "Ne"},
      {11, "sodium", "Na"},       {12, "magnesium", "Mg"},
      {13, "aluminum", "Al"},     {14, "silicon", "Si"},
      {15, "phosphorus", "P"},    {16, "sulfur", "S"},
      {17, "chlorine", "Cl"},     {18, "argon", "Ar"},
      {19, "potassium", "K"},     {20, "calcium", "Ca"},
      {21, "scandium", "Sc"},     {22, "titanium", "Ti"},
      {23, "vanadium", "V"},      {24, "chromium", "Cr"},
      {25, "manganese", "Mn"},    {26, "iron", "Fe"},
      {27, "cobalt", "Co"},       {28, "nickel", "Ni"},
      {29, "copper", "Cu"},       {30, "zinc", "Zn"},
      {31, "gallium", "Ga"},      {32, "germanium", "Ge"},
      {33, "arsenic", "As"},      {34, "selenium", "Se"},
      {35, "bromine", "Br"},      {36, "krypton", "Kr"},
      {37, "rubidium", "Rb"},     {38, "strontium", "Sr"},
      {39, "yttrium", "Y"},       {40, "zirconium", "Zr"},
      {41, "niobium", "Nb"},      {42, "molybdenum", "Mo"},
      {43, "technetium", "Tc"},   {44, "ruthenium", "Ru"},
      {45, "rhodium", "Rh"},      {46, "palladium", "Pd"},
      {47, "silver", "Ag"},       {48, "cadminium", "Cd"},
      {49, "indium", "In"},       {50, "tin", "Sn"},
      {51, "antimony", "Sb"},     {52, "tellurium", "Te"},
      {53, "iodine", "I"},        {54, "xenon", "Xe"},
      {55, "cesium", "Cs"},       {56, "barium", "Ba"},
      {57, "lanthanium", "La"},   {58, "cerium", "Ce"},
      {59, "praseodymium", "Pr"}, {60, "neodymium", "Nd"},
      {61, "promethium", "Pm"},   {62, "samarium", "Sm"},
      {63, "europium", "Eu"},     {64, "gadolinium", "Gd"},
      {65, "terbium", "Tb"},      {66, "dysprosium", "Dy"},
      {67, "holmium", "Ho"},      {68, "erbium", "Er"},
      {69, "thulium", "Tm"},      {70, "ytterbium", "Yb"},
      {71, "lutetium", "Lu"},     {72, "hafnium", "Hf"},
      {73, "tantalum", "Ta"},     {74, "tungsten", "W"},
      {75, "rhenium", "Re"},      {76, "osmium", "Os"},
      {77, "iridium", "Ir"},      {78, "platinum", "Pt"},
      {79, "gold", "Au"},         {80, "mercury", "Hg"},
      {81, "thallium", "Tl"},     {82, "lead", "Pb"},
      {83, "bismuth", "Bi"},      {84, "polonium", "Po"},
      {85, "astatine", "At"},     {86, "radon", "Rn"},
      {87, "francium", "Fr"},     {88, "radium", "Ra"},
      {89, "actinium", "Ac"},     {90, "thorium", "Th"},
      {91, "protactinium", "Pa"}, {92, "uranium", "U"},
      {93, "neptunium", "Np"},    {94, "plutonium", "Pu"},
      {95, "americium", "Am"},    {96, "curium", "Cm"},
      {97, "berkelium", "Bk"},    {98, "californium", "Cf"},
      {99, "einsteinum", "Es"},   {100, "fermium", "Fm"},
      {101, "mendelevium", "Md"}, {102, "nobelium", "No"},
      {103, "lawrencium", "Lr"},  {104, "rutherfordium", "Rf"},
      {105, "dubnium", "Db"},     {106, "seaborgium", "Sg"},
      {107, "bohrium", "Bh"},     {108, "hassium", "Hs"},
      {109, "meitnerium", "Mt"},  {110, "darmstadtium", "Ds"},
      {111, "roentgenium", "Rg"}, {112, "copernicium", "Cn"},
      {113, "nihonium", "Nh"},    {114, "flerovium", "Fl"},
      {115, "moscovium", "Mc"},   {116, "livermorium", "Lv"},
      {117, "tennessine", "Ts"},  {118, "oganesson", "Og"}};
  return element_info;
}
}  // namespace chemistry
}  // namespace libint2

#endif  // INCLUDE_LIBINT2_CHEMISTRY_ELEMENTS_H_
