#pragma once
/**
 * @file core/views.h
 * @brief Non-owning views of 3D and 4D arrays, the parameter type of the field kernels
 *
 * A kernel receives arrays from three kinds of owner: carray owners, ctensor, ctetrad, the
 * per-material slices of a material_field, and null placeholders for arrays a material does not
 * use. A view is a plain pointer to the row-pointer table: it converts from all three and indexes
 * like the raw pointer (`v[i][j][k]`, `v[m][i][j][k]`). `row(i, j)` is the contiguous k-row that a
 * vectorised, DAMM, kernel takes with its length.
 *
 * @verbatim
      view3<T> v           v[i]      -> T**   plane i
                           v[i][j]   -> T*    row along k, contiguous,  == v.row(i, j)
                           v[i][j][k]-> T&
   @endverbatim
 */
#include <cstddef>
#include <type_traits>

template <class T>
class view3
{
	T*** _p = nullptr;

public:
	view3() = default;
	view3(T*** p) : _p(p) {}
	/// anything that converts to T*** a ctensor, a material_field slice
	template <class C>
		requires(!std::is_same_v<std::remove_cvref_t<C>, view3> && !std::is_pointer_v<std::remove_cvref_t<C>> && std::is_convertible_v<C&, T***>)
	view3(C& c) : _p(static_cast<T***>(c)) {}

	T** operator[](std::ptrdiff_t i) const { return _p[i]; }
	T* row(std::ptrdiff_t i, std::ptrdiff_t j) const { return _p[i][j]; }
	T*** data() const { return _p; }
	explicit operator bool() const { return _p != nullptr; }
	operator T***() const { return _p; }
};

template <class T>
class view4
{
	T**** _p = nullptr;

public:
	view4() = default;
	view4(T**** p) : _p(p) {}
	/// anything that converts to T**** a ctetrad, a material_field
	template <class C>
		requires(!std::is_same_v<std::remove_cvref_t<C>, view4> && !std::is_pointer_v<std::remove_cvref_t<C>> && std::is_convertible_v<C&, T****>)
	view4(C& c) : _p(static_cast<T****>(c)) {}

	view3<T> operator[](std::ptrdiff_t m) const { return view3<T>(_p[m]); }
	T**** data() const { return _p; }
	explicit operator bool() const { return _p != nullptr; }
	operator T****() const { return _p; }
};
