/*
 *  Copyright (C) 1996-2017 Edward F. Valeev and Justin T. Fermann
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

#ifndef _libint_copyright_h_
#define _libint_copyright_h_

#include <stdio.h>

static inline void copyright(FILE* ofile)
{
  fprintf(ofile,"/*\n");
  fprintf(ofile," *  Copyright (C) 1996-2017 Edward F. Valeev and Justin T. Fermann\n");
  fprintf(ofile," *\n");
  fprintf(ofile," *  This file is part of Libint.\n");
  fprintf(ofile," *\n");
  fprintf(ofile," *  Libint is free software: you can redistribute it and/or modify\n");
  fprintf(ofile," *  it under the terms of the GNU Lesser General Public License as published by\n");
  fprintf(ofile," *  the Free Software Foundation, either version 3 of the License, or\n");
  fprintf(ofile," *  (at your option) any later version.\n");
  fprintf(ofile," *\n");
  fprintf(ofile," *  Libint is distributed in the hope that it will be useful,\n");
  fprintf(ofile," *  but WITHOUT ANY WARRANTY; without even the implied warranty of\n");
  fprintf(ofile," *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n");
  fprintf(ofile," *  GNU Lesser General Public License for more details.\n");
  fprintf(ofile," *\n");
  fprintf(ofile," *  You should have received a copy of the GNU Lesser General Public License\n");
  fprintf(ofile," *  along with Libint.  If not, see <http://www.gnu.org/licenses/>.\n");
  fprintf(ofile," *\n");
  fprintf(ofile," */\n\n");
}

#endif
