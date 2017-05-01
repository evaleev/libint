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

/*-------------------------------------------------------------------------
  This is not real DERIV_LVL, it is just a constant that reflects the fact
  that t-integrals require incremented angular momentum (effectively, 
  derivatives of ERIs).
 -------------------------------------------------------------------------*/
#define MAX_AM 20
#define DERIV_LVL 1
#define DEFAULT_NEW_AM 6
#define DEFAULT_OPT_AM 6
#define DEFAULT_MAX_CLASS_SIZE 785

typedef struct {

  /* Twice the maximum AM for which manager routines need to be generated */
  int new_am;

  /* Twice the maximum AM for which workers need to be generated */
  int opt_am;

  /* Twice the AM for which manager routines are already
     generated and stored in an existing library */
  int old_am;

  /* Max number of integrals to be processed in one worker. If a class
   is larger than this then split the worker into several. */ 
  int max_class_size;

} Libr12Params_t;
