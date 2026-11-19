/** \file   core/array.h
 *  \brief  The solver's arrays: carray, extern/carray, with an empty state, boxed storage and a
 *          memory charge.
 *
 *  carray allocates a full array and its row-pointer tables. On top of it an array here
 *  - can be empty: declared now, assigned later, and released with reset();
 *  - can be boxed: only the rows of a box are stored, the others share one zero row;
 *  - knows its shape and charges its blocks to the memory tally, core/memory.h.
 *  The shape, data and view are kept here, so every layout indexes the same way. The workarounds
 *  for what carray does not offer are listed in extern/README.md.
 */
#ifndef FDTD_CORE_ARRAY_H
#define FDTD_CORE_ARRAY_H

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <functional>
#include <memory>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>

#include <carray.h>

#include "core/memory.h"

namespace fdtd::core
{
template <class T, std::size_t N, std::size_t A = 64>
class array : public ::carray<T, N, A>
{
	using base = ::carray<T, N, A>;

public:
	using view = std::conditional_t<N == 1, T*, std::conditional_t<N == 2, T**, std::conditional_t<N == 3, T***, T****>>>;
	using shape_t = std::array<std::size_t, N>;

	// creating, copying and releasing an array are cold: one out-of-line copy per instantiation
	// instead of the shared-pointer bookkeeping inlined into every owner

	/// empty: holds no memory
	[[gnu::noinline]] array() : array(std::make_index_sequence<N>{}) {}
	[[gnu::noinline]] ~array() {}

	/// full array of the given extents, memory not initialized
	template <std::integral... I>
		requires(sizeof...(I) == N)
	explicit array(I... n)
		: base(n...), _shape{static_cast<std::size_t>(n)...}, _view(base::get()), _data(base::begin()),
		  _count(product(_shape)), _owned(std::make_shared<memory::Charge>(full_bytes(_shape)))
	{
	}

	// carray's copy constructor takes a non-const reference; it only copies shared pointers
	[[gnu::noinline]] array(const array& a)
		: base(static_cast<base&>(const_cast<array&>(a))), _shape(a._shape), _view(a._view), _data(a._data),
		  _count(a._count), _boxed(a._boxed), _owned(a._owned)
	{
	}
	[[gnu::noinline]] array(array&& a) noexcept
		: base(static_cast<base&&>(a)), _shape(std::exchange(a._shape, {})), _view(std::exchange(a._view, nullptr)),
		  _data(std::exchange(a._data, nullptr)), _count(std::exchange(a._count, 0)), _boxed(std::exchange(a._boxed, false)),
		  _owned(std::move(a._owned))
	{
	}
	[[gnu::noinline]] array& operator=(const array& a)
	{
		base::operator=(a);
		_shape = a._shape;
		_view = a._view;
		_data = a._data;
		_count = a._count;
		_boxed = a._boxed;
		_owned = a._owned;
		return *this;
	}
	[[gnu::noinline]] array& operator=(array&& a) noexcept
	{
		base::operator=(static_cast<base&&>(a));
		_shape = std::exchange(a._shape, {});
		_view = std::exchange(a._view, nullptr);
		_data = std::exchange(a._data, nullptr);
		_count = std::exchange(a._count, 0);
		_boxed = std::exchange(a._boxed, false);
		_owned = std::move(a._owned);
		return *this;
	}

	/// rank 3 or 4: only rows (i, j) in [i0, i1] x [j0, j1] of dimensions N-3, N-2 are stored, every k;
	/// the other rows share one zero row. Rank 4: the same box in every slice.
	static array boxed(shape_t shape, std::size_t i0, std::size_t i1, std::size_t j0, std::size_t j1);

	/// row i, as carray indexes it; any integer type, as carray takes it
	template <std::integral I>
	decltype(auto) operator[](I i) const
	{
		return _view[i];
	}

	/// element at (i, j, ...), read only
	template <std::integral... I>
		requires(sizeof...(I) == N)
	const T& operator()(I... ijk) const
	{
		return element(_view, static_cast<std::size_t>(ijk)...);
	}

	view get() const { return _view; }
	operator view() const { return _view; }

	T* begin() const { return _data; }
	T* end() const { return _data + _count; }

	std::size_t size() const { return _view ? product(_shape) : 0; } ///< elements indexed
	std::size_t extent(std::size_t d) const { return _shape[d]; }
	std::size_t stored() const { return _count; } ///< elements stored, fewer than size() when boxed
	T* storage() const { return _data; }

	array& zeros()
	{
		memory::zero_fill(_data, _count);
		return *this;
	}

	/// release this reference to the memory; the array becomes empty
	[[gnu::noinline]] void reset() { *this = array{}; }

	/// a boxed array's shared zero row is still zero; always true otherwise
	bool padding_zero() const
	{
		return !_boxed || std::all_of(end() - _shape[N - 1], end(), [](const T& v) { return v == T{}; });
	}

private:
	// carray has no empty state: every empty array shares one placeholder of one element
	template <std::size_t... I>
	explicit array(std::index_sequence<I...>) : base(placeholder<I...>())
	{
	}
	template <std::size_t... I>
	static base& placeholder()
	{
		static base p(((void)I, std::size_t{1})...);
		return p;
	}

	static std::size_t product(const shape_t& s)
	{
		return std::accumulate(s.begin(), s.end(), std::size_t{1}, std::multiplies<>());
	}

	/// the blocks carray allocates: one pointer table per leading dimension, then the data
	static std::size_t full_bytes(const shape_t& s)
	{
		std::size_t bytes = memory::block(product(s) * sizeof(T), A);
		std::size_t rows = 1;
		for (std::size_t d = 0; d + 1 < N; ++d)
			bytes += memory::block((rows *= s[d]) * sizeof(void*), A);
		return bytes;
	}

	template <class V, class... I>
	static const T& element(V v, std::size_t i, I... rest)
	{
		if constexpr (sizeof...(I) == 0)
			return v[i];
		else
			return element(v[i], rest...);
	}

	shape_t _shape{};
	view _view = nullptr;
	T* _data = nullptr;
	std::size_t _count = 0;
	bool _boxed = false;
	std::shared_ptr<void> _owned; ///< the memory charge; for a boxed array also its storage
};

template <class T, std::size_t N, std::size_t A>
array<T, N, A> array<T, N, A>::boxed(shape_t shape, std::size_t i0, std::size_t i1, std::size_t j0, std::size_t j1)
{
	static_assert(N == 3 || N == 4, "boxed arrays are rank 3 or 4");
	using tables = std::vector<unique_ptr_aligned<void>>;
	const std::size_t nt = N == 4 ? shape[0] : 1;
	const std::size_t ni = shape[N - 3];
	const std::size_t nj = shape[N - 2];
	const std::size_t nk = shape[N - 1];
	const std::size_t bi = i1 - i0 + 1;
	const std::size_t bj = j1 - j0 + 1;
	const std::size_t count = nt * bi * bj * nk + nk; // the box rows, then the zero row
	auto table = [&](tables& t, auto* p) { t.emplace_back(p, &pdelete); return p; };

	struct Storage
	{
		explicit Storage(std::size_t bytes) : charge(bytes) {}
		memory::Charge charge;
		std::shared_ptr<T[]> data;
		tables pointers;
	};
	const std::size_t bytes = memory::block(count * sizeof(T), A) + memory::block((nt * bi + 1) * nj * sizeof(T*), A) +
							  memory::block(nt * ni * sizeof(T**), A) + (N == 4 ? memory::block(nt * sizeof(T***), A) : 0);
	auto s = std::make_shared<Storage>(bytes);
	s->data = make_shared_aarray<T>(A, count);
	T* data = s->data.get();
	memory::zero_fill(data, count);
	T* zero_row = data + count - nk;

	T** rows = table(s->pointers, palign<T*>(A, (nt * bi + 1) * nj));
	T** zero_plane = rows + nt * bi * nj;
	std::fill_n(zero_plane, nj, zero_row);
	T*** planes = table(s->pointers, palign<T**>(A, nt * ni));
	for (std::size_t t = 0; t < nt; ++t)
		for (std::size_t i = 0; i < ni; ++i)
		{
			if (i < i0 || i > i1)
			{
				planes[t * ni + i] = zero_plane;
				continue;
			}
			T** r = rows + (t * bi + (i - i0)) * nj;
			for (std::size_t j = 0; j < nj; ++j)
				r[j] = (j < j0 || j > j1) ? zero_row : data + ((t * bi + (i - i0)) * bj + (j - j0)) * nk;
			planes[t * ni + i] = r;
		}

	array a;
	if constexpr (N == 3)
		a._view = planes;
	else
	{
		T**** tetrad = table(s->pointers, palign<T***>(A, nt));
		for (std::size_t t = 0; t < nt; ++t)
			tetrad[t] = planes + t * ni;
		a._view = tetrad;
	}
	a._shape = shape;
	a._data = data;
	a._count = count;
	a._boxed = true;
	a._owned = std::move(s);
	return a;
}
} // namespace fdtd::core

namespace fdtd
{
template <class T>
using array1 = core::array<T, 1>; ///< vector
template <class T>
using array2 = core::array<T, 2>; ///< matrix
template <class T>
using array3 = core::array<T, 3>; ///< tensor
template <class T>
using array4 = core::array<T, 4>; ///< tetrad

/// zero-filled (nx+1) x (ny+1) x (nz+1): one ghost cell per axis, the solver's field shape
template <class T>
array3<T> field3(int nx, int ny, int nz)
{
	array3<T> a(nx + 1, ny + 1, nz + 1);
	a.zeros();
	return a;
}

/// n fields3 in one block, [n][nx+1][ny+1][nz+1], zero-filled; at least one
template <class T>
array4<T> field4(int n, int nx, int ny, int nz)
{
	array4<T> a(n > 0 ? n : 1, nx + 1, ny + 1, nz + 1);
	a.zeros();
	return a;
}

/// zero-filled vector of n elements, at least one
template <class T>
array1<T> vector1(std::size_t n)
{
	array1<T> v(n > 0 ? n : 1);
	v.zeros();
	return v;
}
} // namespace fdtd

#endif // FDTD_CORE_ARRAY_H
