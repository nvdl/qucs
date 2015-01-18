/*
 * tvector.h - simple vector template class definitions
 *
 * Copyright (C) 2004, 2005, 2006 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __TVECTOR_H__
#define __TVECTOR_H__

#include <vector>
#include <cassert>

#include <limits>

#include "precision.h"

namespace qucs {

template <class T>
class tvector;

// Forward declarations of friend functions.
template <class T>
T   scalar (tvector<T>, tvector<T>);
template <class T>
nr_double_t maxnorm (tvector<T>);
template <class T>
nr_double_t norm (tvector<T>);
template <class T>
T   sum (tvector<T>);
template <class T>
tvector<T> conj (tvector<T>);
template <class T>
tvector<T> operator + (tvector<T>, tvector<T>);
template <class T>
tvector<T> operator + (tvector<T>, T);
template <class T>
tvector<T> operator + (T, tvector<T>);
template <class T>
tvector<T> operator - (tvector<T>, tvector<T>);
template <class T>
tvector<T> operator * (tvector<T>, nr_double_t);
template <class T>
tvector<T> operator * (nr_double_t, tvector<T>);
template <class T>
tvector<T> operator * (tvector<T>, tvector<T>);
template <class T>
tvector<T> operator - (tvector<T>);
template <class T>
bool operator < (tvector<T>, tvector<T>);
template <class T>
bool operator > (tvector<T>, tvector<T>);

template <class T>
class tvector
{
 private:
  std::vector<T> data;
 public:
  // Constructor creates an unnamed instance of the tvector class.
  tvector () = default;
  tvector (const tvector &c) = default;
  ~tvector () = default;
  /* Constructor creates an unnamed instance of the tvector class with a
     certain length */
  explicit tvector (const std::size_t n) : data(n) {};
  const tvector& operator = (const tvector &);
  /* Returns the tvector element at the given position. */
  T get (std::size_t i) {
    assert (i >= 0 && i < data.size ());
    return data[i];
  }
  void set (int, T);
  void set (T);
  void set (T, int, int);
  void set (tvector, int, int);
  int  getSize (void) { return (int)data.size (); }
  const std::vector<T> * getData (void) const { return &data; }
  void add (T);
  void clear (void);
  void drop (int);
  void truncate (int);
  void exchangeRows (int, int);
  int  isFinite (void);
  void print (void);
  void reorder (int *);
  int  contains (T, nr_double_t eps = std::numeric_limits<nr_double_t>::epsilon());

  // some basic vector operations
#ifndef _MSC_VER
  friend tvector operator +<> (tvector, tvector);
  friend tvector operator -<> (tvector, tvector);
  friend tvector operator *<> (tvector, nr_double_t);
  friend tvector operator *<> (nr_double_t, tvector);
  friend tvector operator *<> (tvector, tvector);
  friend tvector operator -<> (tvector);
  friend tvector operator +<> (tvector, T);
  friend tvector operator +<> (T, tvector);
#endif

  // other operations
#ifndef _MSC_VER
  friend nr_double_t norm<> (tvector);
  friend nr_double_t maxnorm<> (tvector);
  friend T   sum<> (tvector);
  friend T   scalar<> (tvector, tvector);
  friend tvector     conj<> (tvector);
#endif

  // comparisons
#ifndef _MSC_VER
  friend bool operator < <> (tvector, tvector);
  friend bool operator > <> (tvector, tvector);
#endif

  // intrinsic operators
  tvector operator += (tvector);
  tvector operator -= (tvector);
  tvector operator *= (nr_double_t);
  tvector operator /= (nr_double_t);

  // assignment operators
  tvector operator = (const T);

  // easy accessor operators
  T  operator () (int i) const {
    assert (i >= 0 && i < (int)data.size ()); return (data)[i]; }
  T& operator () (int i) {
    assert (i >= 0 && i < (int)data.size ()); return (data)[i]; }
};

} // namespace qucs

#include "tvector.cpp"

#endif /* __TVECTOR_H__ */
