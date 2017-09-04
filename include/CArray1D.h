#ifndef CARRAY_1D_H
#define CARRAY_1D_H

#include <iostream>
#include <cstddef>

// 1D Array of specified type, data is optionally owned
template<typename T, bool OWNER>
class CArray2D;

template<typename T, bool OWNER=true>
class CArray1D {
 public:
  typedef T                 value_type;
  typedef value_type       &reference;
  typedef value_type const &const_reference;
  typedef value_type       *pointer;
  typedef value_type const *const_pointer;

  typedef size_t            size_type;
  typedef size_t            index_type;
  typedef ptrdiff_t         difference_type;

  typedef pointer       iterator;
  typedef const_pointer const_iterator;

  typedef pointer       reverse_iterator;
  typedef const_pointer const_reverse_iterator;

  typedef reference slice_type;

 public:
  CArray1D(value_type *data, index_type d1) :
   data_(data), d1_(d1), size_(d1) {
    if (OWNER) assert(false && "Owner");
  }

  explicit CArray1D(index_type d1) :
   d1_(d1), size_(d1) {
    if (! OWNER) assert(false && "Not Owner");

    allocate(size_);
  }

  CArray1D(index_type d1, const_reference idata) :
   d1_(d1), size_(d1) {
    if (! OWNER) assert(false && "Not Owner");

    allocate(size_, idata);
  }

 private:
  CArray1D(const CArray1D &array);

  const CArray1D &operator=(const CArray1D &array);

 public:
 ~CArray1D() {
    if (OWNER)
      deallocate();
  }

 private:
  void allocate(size_type size_) {
    data_ = new T [size_];
  }

  void allocate(size_type size_, value_type idata) {
    data_ = new T [size_];

    for (size_type i = 0; i < size_; ++i)
      data_[i] = idata;
  }

  void deallocate() {
    delete [] data_;
  }

 public:
  size_type size() { return size_; }
  bool      empty() { return d1_ == 0; }

  iterator begin() { return &data_[0]; }
  iterator end  () { return &data_[size_]; }

  const_iterator begin() const { return &data_[0]; }
  const_iterator end  () const { return &data_[size_]; }

  slice_type operator[](index_type ind) {
    return slice_type(data_[ind]);
  }

  const slice_type operator[](index_type ind) const {
    return slice_type(data_[ind]);
  }

  value_type at(index_type i1) {
    return data_[i1];
  }

  void print(std::ostream &os) const {
    os << "[";

    for (index_type i = 0; i < d1_; ++i)
      os << " " << data_[i];

    os << " ]";
  }

  friend std::ostream &operator<<(std::ostream &os, const CArray1D &array) {
    array.print(os);

    return os;
  }

  friend class CArray2D<T, false>;
  friend class CArray2D<T, true>;

 private:
  T          *data_ { nullptr };
  index_type  d1_ { 0 };
  size_type   size_ { 0 };
};

#endif
