#ifndef CARRAY_2D_H
#define CARRAY_2D_H

#include <CArray1D.h>

// 2D Array of specified type, data is optionally owned
template<typename T, bool OWNER=true>
class CArray2D {
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

  typedef CArray1D<T, false> slice_type;

 private:
  T          *data_;
  index_type  d1_;
  index_type  d2_;
  size_type   size_;

 public:
  CArray2D(value_type *data, index_type d1, index_type d2) :
   data_(data), d1_(d1), d2_(d2), size_(d1*d2) {
    if (OWNER) assert(false && "Not Owner");
  }

  CArray2D(index_type d1, index_type d2) :
   d1_(d1), d2_(d2), size_(d1*d2) {
    if (! OWNER) assert(false && "Owner");

    allocate(size_);
  }

  CArray2D(index_type d1, index_type d2, const_reference idata) :
   d1_(d1), d2_(d2), size_(d1*d2) {
    if (! OWNER) assert(false && "Owner");

    allocate(size_, idata);
  }

 private:
  CArray2D(const CArray2D &array);

  const CArray2D &operator=(const CArray2D &array);

 public:
 ~CArray2D() {
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
  bool      empty() { return d1_ == 0 || d2_ == 0; }

  iterator begin() { return &data_[0]; }
  iterator end  () { return &data_[size_]; }

  const_iterator begin() const { return &data_[0]; }
  const_iterator end  () const { return &data_[size_]; }

  slice_type operator[](index_type ind) {
    return slice_type(&data_[ind*d2_], d2_);
  }

  const slice_type operator[](index_type ind) const {
    return slice_type(&data_[ind*d2_], d2_);
  }

  void print(std::ostream &os) const {
    os << "[";

    index_type k = 0;

    for (index_type i = 0; i < d1_; ++i) {
      os << " [";

      for (index_type j = 0; j < d2_; ++j, k++) {
        os << " " << data_[k];
      }

      os << " ]";
    }

    os << " ]";
  }

  friend std::ostream &operator<<(std::ostream &os, const CArray2D &array) {
    array.print(os);

    return os;
  }

  value_type at(index_type i1, index_type i2) {
    return data_[i1*d2_ + i2];
  }
};

#endif
