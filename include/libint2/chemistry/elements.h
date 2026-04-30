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

#include <cmath>
#include <stdexcept>
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
/// Most-abundant isotope mass number A for elements Z=1..118.
/// Used by gaussian_nuclear_exponent() to compute nuclear charge exponents.
/// Source: Visscher & Dyall, At. Data Nucl. Data Tables 67, 207 (1997),
/// Table I, extended to Z=118 using most-abundant/stable isotopes.
// clang-format off
static constexpr unsigned short most_abundant_isotope_A[] = {
  /*  0 */   0,
  /*  1 H  */   1,  /*  2 He */   4,  /*  3 Li */   7,  /*  4 Be */   9,
  /*  5 B  */  11,  /*  6 C  */  12,  /*  7 N  */  14,  /*  8 O  */  16,
  /*  9 F  */  19,  /* 10 Ne */  20,  /* 11 Na */  23,  /* 12 Mg */  24,
  /* 13 Al */  27,  /* 14 Si */  28,  /* 15 P  */  31,  /* 16 S  */  32,
  /* 17 Cl */  35,  /* 18 Ar */  40,  /* 19 K  */  39,  /* 20 Ca */  40,
  /* 21 Sc */  45,  /* 22 Ti */  48,  /* 23 V  */  51,  /* 24 Cr */  52,
  /* 25 Mn */  55,  /* 26 Fe */  56,  /* 27 Co */  59,  /* 28 Ni */  58,
  /* 29 Cu */  63,  /* 30 Zn */  64,  /* 31 Ga */  69,  /* 32 Ge */  74,
  /* 33 As */  75,  /* 34 Se */  80,  /* 35 Br */  79,  /* 36 Kr */  84,
  /* 37 Rb */  85,  /* 38 Sr */  88,  /* 39 Y  */  89,  /* 40 Zr */  90,
  /* 41 Nb */  93,  /* 42 Mo */  98,  /* 43 Tc */  98,  /* 44 Ru */ 102,
  /* 45 Rh */ 103,  /* 46 Pd */ 106,  /* 47 Ag */ 107,  /* 48 Cd */ 114,
  /* 49 In */ 115,  /* 50 Sn */ 120,  /* 51 Sb */ 121,  /* 52 Te */ 130,
  /* 53 I  */ 127,  /* 54 Xe */ 132,  /* 55 Cs */ 133,  /* 56 Ba */ 138,
  /* 57 La */ 139,  /* 58 Ce */ 140,  /* 59 Pr */ 141,  /* 60 Nd */ 144,
  /* 61 Pm */ 145,  /* 62 Sm */ 152,  /* 63 Eu */ 153,  /* 64 Gd */ 158,
  /* 65 Tb */ 159,  /* 66 Dy */ 162,  /* 67 Ho */ 165,  /* 68 Er */ 168,
  /* 69 Tm */ 169,  /* 70 Yb */ 174,  /* 71 Lu */ 175,  /* 72 Hf */ 180,
  /* 73 Ta */ 181,  /* 74 W  */ 184,  /* 75 Re */ 187,  /* 76 Os */ 192,
  /* 77 Ir */ 193,  /* 78 Pt */ 195,  /* 79 Au */ 197,  /* 80 Hg */ 202,
  /* 81 Tl */ 205,  /* 82 Pb */ 208,  /* 83 Bi */ 209,  /* 84 Po */ 209,
  /* 85 At */ 210,  /* 86 Rn */ 222,  /* 87 Fr */ 223,  /* 88 Ra */ 226,
  /* 89 Ac */ 227,  /* 90 Th */ 232,  /* 91 Pa */ 231,  /* 92 U  */ 238,
  /* 93 Np */ 237,  /* 94 Pu */ 244,  /* 95 Am */ 243,  /* 96 Cm */ 247,
  /* 97 Bk */ 247,  /* 98 Cf */ 251,  /* 99 Es */ 252, /* 100 Fm */ 257,
  /*101 Md */ 258, /* 102 No */ 259, /* 103 Lr */ 262, /* 104 Rf */ 261,
  /*105 Db */ 262, /* 106 Sg */ 263, /* 107 Bh */ 262, /* 108 Hs */ 265,
  /*109 Mt */ 266, /* 110 Ds */ 281, /* 111 Rg */ 282, /* 112 Cn */ 285,
  /*113 Nh */ 286, /* 114 Fl */ 289, /* 115 Mc */ 290, /* 116 Lv */ 293,
  /*117 Ts */ 294, /* 118 Og */ 294
};
// clang-format on

/// Gaussian nuclear model exponent ξ for element Z, in bohr^{-2}.
/// From Visscher & Dyall, At. Data Nucl. Data Tables 67, 207 (1997):
///   √⟨R²⟩ = (0.836·A^(1/3) + 0.570) fm,  ξ = 3/(2⟨R²⟩)
/// where A is the most-abundant isotope mass number.
/// @param Z atomic number (1..118)
/// @return ξ in bohr^{-2}
/// @throw std::out_of_range if Z is out of range
inline double gaussian_nuclear_exponent(int Z) {
  if (Z < 1 || Z > 118)
    throw std::out_of_range("gaussian_nuclear_exponent: Z=" +
                            std::to_string(Z) + " out of range [1,118]");
  const double A = most_abundant_isotope_A[Z];
  // RMS radius in fm: sqrt(<R^2>) = 0.836 * A^(1/3) + 0.570
  const double rms_fm = 0.836 * std::cbrt(A) + 0.570;
  // Convert fm to bohr: 1 fm = 1e-15 m, 1 bohr = 0.529177210903e-10 m
  constexpr double fm_per_bohr = 52917.72109030;
  const double rms_bohr = rms_fm / fm_per_bohr;
  // xi = 3 / (2 * <R^2>)
  return 1.5 / (rms_bohr * rms_bohr);
}

/// Gaussian nuclear model exponent ξ for a given mass number A, in bohr^{-2}.
/// @sa gaussian_nuclear_exponent(int Z)
/// @param A mass number (must be positive and finite)
/// @throw std::invalid_argument if A is not positive or is not finite
inline double gaussian_nuclear_exponent_from_A(double A) {
  if (!std::isfinite(A) || A <= 0.0)
    throw std::invalid_argument(
        "gaussian_nuclear_exponent_from_A: A=" + std::to_string(A) +
        " is invalid; must be positive and finite");
  const double rms_fm = 0.836 * std::cbrt(A) + 0.570;
  constexpr double fm_per_bohr = 52917.72109030;
  const double rms_bohr = rms_fm / fm_per_bohr;
  return 1.5 / (rms_bohr * rms_bohr);
}

}  // namespace chemistry
}  // namespace libint2

#endif  // INCLUDE_LIBINT2_CHEMISTRY_ELEMENTS_H_
