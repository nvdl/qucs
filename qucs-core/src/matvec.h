/*
 * matvec.h - matrix vector class definitions
 *
 * Copyright (C) 2004, 2005, 2006, 2007, 2009 Stefan Jahn <stefan@lkcc.org>
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this package; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * $Id$
 *
 */

#ifndef __MATVEC_H__
#define __MATVEC_H__

#include <string>
#include <vector>

#include "matrix.h"


namespace qucs {

class matvec;

// forward declarations
matvec transpose (const matvec&);
matvec conj      (const matvec&);
vector det       (const matvec&);
matvec inverse   (const matvec&);
matvec sqr       (const matvec&);
matvec pow       (const matvec&, int);
matvec pow       (const matvec&, const qucs::vector &);
matvec twoport   (const matvec&, char, char);
matvec real      (const matvec&);
matvec imag      (const matvec&);
matvec abs       (const matvec&);
matvec dB        (const matvec&);
matvec arg       (const matvec&);
matvec adjoint   (const matvec&);
vector rollet    (const matvec&);
vector b1        (const matvec&);
matvec rad2deg       (const matvec&);
matvec deg2rad       (const matvec&);
matvec stos (const matvec&, nr_complex_t, nr_complex_t z0 = 50.0);
matvec stos (const matvec&, nr_double_t, nr_double_t z0 = 50.0);
matvec stos (const matvec&, const qucs::vector &, nr_complex_t z0 = 50.0);
matvec stos (const matvec&, nr_complex_t, const qucs::vector &);
matvec stos (const matvec&, const qucs::vector &, const qucs::vector &);
matvec stoz (const matvec&, nr_complex_t z0 = 50.0);
matvec stoz (const matvec&, const qucs::vector &);
matvec ztos (const matvec&, nr_complex_t z0 = 50.0);
matvec ztos (const matvec&, const qucs::vector &);
matvec ztoy (const matvec&);
matvec stoy (const matvec&, nr_complex_t z0 = 50.0);
matvec stoy (const matvec&, const qucs::vector &);
matvec ytos (const matvec&, nr_complex_t z0 = 50.0);
matvec ytos (const matvec&, const qucs::vector &);
matvec ytoz (const matvec&);
 
class matvec : public std::vector<matrix>
{
 private:
  typedef decltype(matrix().rows()) index;
  index r = 0;
  index c = 0;
  std::string name = "";
 public:
  using std::vector<matrix>::vector;
  matvec(int s,int rr, int cc) : vector(s,matrix(rr,cc)), r(rr), c(cc) {
  }
  index cols (void) const { return c; }
  index rows (void) const { return r; }
  void setName (const std::string &);
  std::string getName (void) const;
  
  qucs::vector get (int, int);
  static char * createMatrixString (const char *, int, int);
  static char * createMatrixString (char, int, int);
  static char * isMatrixVector (const char *, int&, int&);
  static matvec * getMatrixVector (qucs::vector *, char *);
  static void getMatrixVectorSize (qucs::vector *, char *, int&, int&, int&);

  // intrinsic operator functions
  matvec operator  - ();
  matvec operator += (const matvec&);
  matvec operator -= (const matvec&);
};



  
} // namespace qucs

#endif /* __MATVEC_H__ */
