#ifndef __HOPPET___BASE_TYPES__
#define __HOPPET___BASE_TYPES__

#include <tuple>

/// an empty class as a default for template parameters    
struct Empty {};

//-----------------------------------------------------------------------------
/// @brief Generic object view class, to be used as a base for non-owning wrappers
/// of Fortran objects that do not provide a view of the underlying data
///
/// @tparam T  the pointer type to the underlying Fortran object
/// @tparam E  an extra type to hold additional information
template<typename T, typename E = Empty>
class obj_view {
public:
  using ptr_t = T;
  using extras_t = E;
  using view_t = obj_view<T,E>;

  //obj_view() noexcept {}
  obj_view() noexcept = default; 
  obj_view(T * ptr) noexcept : _ptr(ptr) {}
  obj_view(T * ptr, const E & extra) noexcept : _ptr(ptr), _extra(extra) {}
  obj_view(const E & extra) noexcept : _ptr(nullptr), _extra(extra) {}
  obj_view(const obj_view & other) noexcept : _ptr(other._ptr), _extra(other._extra) {}
  obj_view & operator= (const obj_view & other) noexcept {
    _ptr   = other._ptr;
    _extra = other._extra;
    return *this;
  }

  /// @brief  make this object a view of the "other" object
  /// @param other the object of which to take a view
  void take_view(const obj_view & other) noexcept {
    _ptr   = other._ptr;
    _extra = other._extra;
  }

  const E & extra() const { return _extra; }
  E & extra() { return _extra; }

  /// @brief return a const pointer to the underlying Fortran object
  const T * ptr() const { return _ptr; }
  /// @brief return a non-const pointer to the underlying Fortran object
  T * ptr() { return _ptr; }

  /// @brief  return a validated const pointer to the underlying Fortran object
  ///
  /// If the pointer is null, an exception is thrown
  const T * valid_ptr() const { ensure_valid(); return _ptr; }

  /// @brief  return a validated non-const pointer to the underlying Fortran object
  ///
  /// If the pointer is null, an exception is thrown
  T * valid_ptr() { ensure_valid(); return _ptr; }

  void ensure_valid() const {if (!_ptr) throw std::runtime_error("hoppet::obj_view::ensure_valid: object pointer is null");}
protected:
  /// @brief pointer to the underlying Fortran grid_def object
  T * _ptr = nullptr;
  E _extra;
};

//-----------------------------------------------------------------------------
/// @brief Generic object-owning class, to be used as a base for owning wrappers
/// of Fortran objects that do not provide a view of the underlying data
///
/// @tparam V  the obj_view type to be wrapped
template<typename V>
//requires (requires { typename V::ptr_t; } && std::derived_from<V, obj_view<typename V::ptr_t, typename V::extras_t>>)
class obj_owner : public V {
public:
  using view_t = V;
  using ptr_t    = typename V::ptr_t;
  using extras_t = typename V::extras_t;

  obj_owner() noexcept = default;
  obj_owner(ptr_t * ptr) noexcept : V(ptr) {}
  obj_owner(ptr_t * ptr, extras_t extras) noexcept : view_t(ptr, extras) {}
  obj_owner(const view_t & obj) : view_t(generic_copy(obj.ptr()),obj.extra()) {}
  obj_owner(const obj_owner & obj) : view_t(generic_copy(obj.ptr()),obj.extra()) {}
  obj_owner & operator= (const view_t & obj) {
    if (this != &obj) {
      generic_delete(this->_ptr);
      this->_ptr = generic_copy(obj.ptr());
      this->_extra = obj._extra;
    }
    return *this;
  }
  obj_owner & operator= (const obj_owner & obj) {
    if (this != &obj) {
      generic_delete(this->_ptr);
      this->_ptr = generic_copy(obj.ptr());
      this->_extra = obj._extra;
    }
    return *this;
  }

  obj_owner(obj_owner && other) noexcept : view_t(other.ptr(), other.extra()) {
    // null out the other pointer to avoid double deletion
    other._ptr = nullptr;
  }
  obj_owner & operator= (obj_owner && other) noexcept {
    if (this != &other) {
      generic_delete(this->_ptr);
      this->_ptr = other.ptr();
      this->_extra = other._extra;
      other._ptr = nullptr;
    }
    return *this;
  }
  virtual ~obj_owner() {generic_delete(V::ptr());}

  void take_view() noexcept = delete; // obj_owner cannot safely take_view(), so disable it

};


//-----------------------------------------------------------------------------
/// @brief Generic object view class, to be used as a base for non-owning wrappers
/// of Fortran objects that provide a view of the underlying double-precision data
///
/// @tparam E  a type for holding extra information
template<typename E, typename D = double>
class data_view {
public:

  using extras_t = E;
  using data_t = D;
  using view_t = data_view<E,D>;

  data_view() noexcept {}

  data_view(D * data_ptr, std::size_t size, const extras_t & extras) noexcept
    : _data(data_ptr), _size(size), _extras(extras) {}

  explicit data_view(const data_view<E, D> & other) noexcept {
    take_view(other);
  }

  data_view<E, D> & operator=(const data_view<E, D> & other) {
    this->copy_data(other);
    return *this;
  }

  void take_view(const data_view<E, D> & other) noexcept {
    _data = other._data;
    _size = other._size;
    _extras = other._extras;
  }

  /// @brief  return pointer to the underlying data
  D       * data()       {return _data;}

  /// @brief  return a const pointer to the underlying data
  const D * data() const {return _data;}

  /// @brief  return the size of the data array
  std::size_t    size() const {return _size;}

  data_view & set_data_ptr(D * data_in) noexcept {_data = data_in; return *this;}
  data_view & set_size(std::size_t size_in) noexcept {_size = size_in; return *this;}
  data_view & set_extras(const E & extras_in) noexcept {_extras = extras_in; return *this;}
  /// compound assignment arithmetic operators
  ///@{
  data_view<E, D> & operator+=(const data_view<E, D> & other) {
    auto [sz, this_data, other_data] = prepare_compound(other);
    for (std::size_t iy=0; iy<sz; ++iy) this_data[iy] += other_data[iy];
    return *this;
  }

  data_view<E, D> & operator-=(const data_view<E, D> & other) {
    auto [sz, this_data, other_data] = prepare_compound(other);
    for (std::size_t iy=0; iy<sz; ++iy) this_data[iy] -= other_data[iy];
    return *this;
  }

  data_view<E, D> & operator*=(double val) {
    std::size_t sz = size();
    D * this_data = data();
    for (std::size_t iy=0; iy<sz; ++iy) this_data[iy] *= val;
    return *this;
  }
  data_view<E, D> & operator/=(double val) {
    std::size_t sz = size();
    D * this_data = data();
    for (std::size_t iy=0; iy<sz; ++iy) this_data[iy] /= val;
    return *this;
  }
  ///@}

  /// copy the data from other, assuming *this is initialised and of the correct size, etc.
  void copy_data(const data_view<E, D> & other) {
    auto [sz, this_data, other_data] = prepare_compound(other);
    std::copy(other_data, other_data + sz, this_data);
  }

  const E & extras() const {return _extras;}
  E & extras() {return _extras;}

  void reset() {
    _data = nullptr;
    _size = 0;
    _extras = E();
  }

protected:
  D *           _data = nullptr;
  std::size_t   _size = 0;
  E _extras = E();

  /// @brief  prepare for a compound assignment operation
  /// @param b  the other data_view to operate on
  /// @return  a tuple containing the size, this data pointer, and the other data pointer
  inline std::tuple<std::size_t, double *, const double *> prepare_compound(const data_view<E> & b ) {
    extras().ensure_compatible(b.extras());
    return std::make_tuple(size(), data(), b.data());
  }
};


//-----------------------------------------------------------------------------
/// @brief Generic data_owner class, to be used as a base for owning wrappers
/// of Fortran objects that provide a view of the underlying data (typically double-precision)
///
/// @tparam V  the data_view type from which this derives
/// @tparam P  the pointer type to the underlying Fortran object that owns the data
///
template<typename V, typename P>
class data_owner : public V {
public:

  using base_t = data_owner<V,P>;
  using view_t = V;
  using ptr_t  = P;

  data_owner() {}
  virtual ~data_owner() {del();}
  void del() {generic_delete(_ptr); reset_ptrs();}

  virtual void alloc_virtual(std::size_t sz, const typename V::extras_t & extras) = 0;

  data_owner & operator=(const data_owner & other) {
    if (_ptr == other._ptr) return *this; // self-assignment check
    copy(other);
    return *this;
  }

  /// copy constructor: perform a deep copy using the `copy` helper
  data_owner(const data_owner & other) {
    _ptr = nullptr;
    copy(other);
  }

  /// @brief move constructor: take ownership of other's data
  /// @param other the other data_owner to move from
  data_owner(data_owner && other) noexcept  {move_no_del(other);}

  data_owner & operator=(data_owner && other) noexcept {
    if (this != &other) {
      del();
      move_no_del(other);
    }
    return *this;
  }

  /// @brief  make a copy of other, including allocating new storage if needed
  /// @param  other the other grid_quant to copy from
  void copy(const V & other) {
    if (_ptr && this->extras() == other.extras()) {
      //std::cout << " reusing existing storage in grid_quant::copy\n";
      this->copy_data(other);
    } else {
      //std::cout << " allocating new storage in grid_quant::copy, size = " << other.size() << "\n";
      del();
      alloc_virtual(other.size(), other.extras());
      this->copy_data(other);
      //std::cout << " data[50] = " << data()[50] << " v " << other.data()[50] << "\n";
    }
  }

  const P* ptr() const { return _ptr; }
  P* ptr()  { return _ptr; }

protected:

  /// @brief  resets the internal pointers to null and the extras to their default, without deleting any data
  void reset_ptrs() {
    V::reset();
    _ptr = nullptr;
  }

  /// @brief  move the contents from other into this, without deleting existing data
  /// @param  other: the other grid_quant to move from
  ///
  /// Note that this does not delete any existing data in *this, so be careful to call del() 
  /// first if needed.
  void move_no_del(data_owner & other) noexcept {
    // move the semantics
    this->take_view(other);
    _ptr = other._ptr;
    other.reset_ptrs();
  }

  using view_t::take_view; //< data_owner cannot safely take_view(), so disable it


  P * _ptr = nullptr;
};


#endif // __HOPPET___BASE_TYPES__