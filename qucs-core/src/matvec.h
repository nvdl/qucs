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

class matvec
{
 private:
  typedef decltype(matrix().rows()) index;
  typedef std::vector<matrix> data_type;
  index rows;
  index cols;
  std::string name;
  data_type data;
 public:
  typedef data_type::value_type value_type;
  typedef data_type::allocator_type allocator_type;
  typedef data_type::reference reference;
  typedef data_type::const_reference const_reference;
  typedef data_type::iterator iterator;
  typedef data_type::const_iterator const_iterator;
  typedef data_type::size_type size_type;

  /* begin iterator */
  iterator begin() noexcept {
    return data.begin();
  }
  const_iterator begin() const noexcept {
    return data.begin();
  }

  /* end iterator */
  iterator end() noexcept {
    return data.begin();
  }
  const_iterator end() const noexcept {
    return data.begin();
  }

  reference operator[] (const size_type n) {
    return data[n];
  };
  const_reference operator[] (const size_type n) const {
    return data[n];
  }
  
  /*! default constructor: empty */
  matvec (): rows(0), cols(0), name(), data() {};
  /*! Constructor creates an unnamed instance of the matvec class with a
   certain number of empty matrices. */
  matvec (decltype(data.size()) size, index r, index c):
	  rows(r),
	  cols(c),
	  name(),
	  data(size,matrix(r,c))
	  {};
  matvec (const matvec &) = default;
  ~matvec () = default;
  std::size_t size (void) const { return data.size(); }
  int getCols (void) const { return cols; }
  int getRows (void) const { return rows; }
  void setName (const std::string &);
  std::string getName (void) const;
  void set (qucs::vector, int, int);
  void set (matrix, int);
  qucs::vector get (int, int);
  static char * createMatrixString (const char *, int, int);
  static char * createMatrixString (char, int, int);
  static char * isMatrixVector (const char *, int&, int&);
  static matvec * getMatrixVector (qucs::vector *, char *);
  static void getMatrixVectorSize (qucs::vector *, char *, int&, int&, int&);

  // operator functions
  friend matvec operator + (const matvec&, const matvec&);
  friend matvec operator + (const matvec&, const matrix&);
  friend matvec operator + (const matrix&, const matvec&);
  friend matvec operator + (const matvec&, const nr_complex_t);
  friend matvec operator + (const nr_complex_t, const matvec&);
  friend matvec operator + (const matvec&, const nr_double_t);
  friend matvec operator + (const nr_double_t, const matvec&);
  friend matvec operator + (const matvec&, const qucs::vector&);
  friend matvec operator + (const qucs::vector&, const matvec&);
  friend matvec operator - (const matvec&, const matvec&);
  friend matvec operator - (const matvec&, const matrix&);
  friend matvec operator - (const matrix&, const matvec&);
  friend matvec operator - (const matvec&, nr_complex_t);
  friend matvec operator - (nr_complex_t, const matvec&);
  friend matvec operator - (const matvec&, nr_double_t);
  friend matvec operator - (nr_double_t, const matvec&);
  friend matvec operator - (const matvec&, const qucs::vector&);
  friend matvec operator - (const qucs::vector&, const matvec&);
  friend matvec operator / (const matvec&, const qucs::vector&);
  friend matvec operator * (const matvec&, const qucs::vector&);
  friend matvec operator * (const qucs::vector&, const matvec&);
  friend matvec operator * (const matvec&, nr_complex_t);
  friend matvec operator * (nr_complex_t, const matvec&);
  friend matvec operator * (const matvec&, nr_double_t);
  friend matvec operator * (nr_double_t, const matvec&);
  friend matvec operator * (const matvec&, const matvec&);
  friend matvec operator * (const matvec&, const matrix&);
  friend matvec operator * (const matrix&, const matvec&);

  // intrinsic operator functions
  matvec operator  - ();
  matvec operator += (const matvec&);
  matvec operator -= (const matvec&);

  // other operations
  friend matvec transpose (const matvec&);
  friend matvec conj      (const matvec&);
  friend qucs::vector det       (const matvec&);
  friend matvec inverse   (const matvec&);
  friend matvec sqr       (const matvec&);
  friend matvec pow       (const matvec&, int);
  friend matvec pow       (const matvec&, const qucs::vector&);
  friend matvec twoport   (const matvec&, char, char);
  friend matvec real      (const matvec&);
  friend matvec imag      (const matvec&);
  friend matvec abs       (const matvec&);
  friend matvec dB        (const matvec&);
  friend matvec arg       (const matvec&);
  friend matvec adjoint   (const matvec&);
  friend qucs::vector rollet    (const matvec&);
  friend qucs::vector b1        (const matvec&);
  friend matvec rad2deg       (const matvec&);
  friend matvec deg2rad       (const matvec&);

  friend matvec stos (const matvec&, nr_complex_t, nr_complex_t);
  friend matvec stos (const matvec&, nr_double_t, nr_double_t);
  friend matvec stos (const matvec&, const qucs::vector&, nr_complex_t);
  friend matvec stos (const matvec&, nr_complex_t, const qucs::vector&);
  friend matvec stos (const matvec&, const qucs::vector&, const qucs::vector&);
  friend matvec stoz (const matvec&, nr_complex_t);
  friend matvec stoz (const matvec&, const qucs::vector&);
  friend matvec ztos (const matvec&, nr_complex_t);
  friend matvec ztos (const matvec&, const qucs::vector&);
  friend matvec ztoy (const matvec&);
  friend matvec stoy (const matvec&, nr_complex_t);
  friend matvec stoy (const matvec&, const qucs::vector&);
  friend matvec ytos (const matvec&, nr_complex_t);
  friend matvec ytos (const matvec&, const qucs::vector&);
  friend matvec ytoz (const matvec&);

};



  
} // namespace qucs

#endif /* __MATVEC_H__ */
