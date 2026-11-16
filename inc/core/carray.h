#ifndef __C_ARRAY_H__
#define __C_ARRAY_H__
// Copyright 2025  Constantine Papakonstantinou
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files, the "Software", to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
/**
 * \file core/carray.h header only support for dynamically allocated memory aligned arrays
 * \author cpapakonstantinou
 * \date 2025
 **/
#include <algorithm>
#include <array>
#include <cstdint>
#include <numeric>
#include <memory>
#include <type_traits>
#include <assert.h>
#include "core/aligned_memory.h"

/**
* 	\brief	dynamically allocated and aligned contiguous memory arrays.
*	Row major contiguous memory allocation.
*	Memory alignment can be specified by template instantiation.
*	T - Type
*	N - Rank
*	A - Alignment
*
*	\note Memory is NOT initialized; call zeros(), or use field3(), when the
*	      legacy MemAllocat*Array zero-fill semantics are required.
*/
template<class T, size_t N, size_t A>
class carray
{
	static_assert(N >= 1 && N <= 4, "carray supports rank <= 4");

	public:
		/** \brief Empty array: owns no memory, declare now, assign a sized carray later. */
		carray() = default;

		/**
		 *	\brief constructor.
		 *	Variadic constructor which accepts a parameter pack.
		 *	Parameter pack must expand to N elements.
		 *
		 *	\param 	ijk	parameter pack
		 */
		template<class... IJK>
		explicit
		carray(IJK&&... ijk):
		_shape(new size_t[N]{ static_cast<size_t>(ijk)... })
		{
			static_assert(sizeof...(IJK) == N, "number of indices must match rank of array");
			allocate_memory();
		}

		/**
		*	\brief copy constructor.
		*	Internally carray uses a shared_ptr to manage the buffer.
		*	Copy construction is safely supported
		*	\param 	an initialized carray class
		*/
		carray(carray<T, N, A>& a):
		_shape(a._shape),
		_buffer(a._buffer),
		_ptr(a._ptr),
		_count(a._count),
		_boxed(a._boxed)
		{}

		/**
		*	\brief	move constructor.
		*	Internally carray uses a shared_ptr to manage the buffer.
		*	Move construction is safely supported.
		*	\param 	i	number of rows
		*/
		carray(carray<T, N, A>&& a) noexcept:
		_shape(std::move(a._shape)),
		_buffer(std::move(a._buffer)),
		_ptr(std::move(a._ptr)),
		_count(a._count),
		_boxed(a._boxed)
		{}

		/**	\brief	Copy assignment operator. */
		carray<T, N, A>&
		operator=(const carray<T, N, A>& a)
		{
			if (&a != this)
			{
				_shape=a._shape;
				_buffer=a._buffer;
				_ptr = a._ptr;
				_count = a._count;
				_boxed = a._boxed;
			}
			return *this;
		}

		/**	\brief Move assignment operator. */
		carray<T, N, A>&
		operator=(carray<T, N, A>&& a) noexcept
		{
			if (&a != this)
			{
				_shape = std::move(a._shape);
				_buffer = std::move(a._buffer);
				_ptr = std::move(a._ptr);
				_count = a._count;
				_boxed = a._boxed;
			}
			return *this;
		}

	private:

		using vector_t = T*; ///< vector datatype
		using matrix_t = T**; ///< matrix datatype
		using tensor_t = T***; ///< tensor datatype
		using tetrad_t = T****; ///< tetrad datatype
		using shape_t = std::shared_ptr<size_t[]>;///< shared array shape datatype
		using buffer_t = std::shared_ptr<T[]>; ///< shared contiguous array datatype underlying the memory of the carray
		using ptr_t = std::shared_ptr<void>; ///< Hierarchical pointer datatype

		shape_t _shape; ///< shape of the carray
		buffer_t _buffer; ///< contiguous memory of the carray
		ptr_t _ptr; ///< Hierarchical pointer for structured memory access to buffer
		size_t _count = 0; ///< elements stored in _buffer, size(), or a box and its zero row
		bool _boxed = false; ///< rows outside a box point to one shared zero row, see boxed()

		inline void
		allocate_memory()
		{
			// size_t accumulator: an int initial value would make std::accumulate
			// compute in int and overflow beyond 2^31 elements, large 3D/4D fields.
			const size_t n = size();
			_count = n;
			_buffer = make_shared_aarray<T>(A, n);

			if constexpr (N == 1)
			{
				 _ptr = std::shared_ptr<void>(_buffer.get(), [](void*){});
			}
			else if constexpr (N == 2)
			{
				matrix_t matrix = static_cast<matrix_t>(palign<vector_t>(A, _shape[0]));
				for (size_t i = 0; i < _shape[0]; ++i)
					matrix[i] = &_buffer[i * _shape[1]];
				_ptr = std::shared_ptr<void>(matrix, &pdelete);
			}
			else if constexpr (N == 3)
			{
				tensor_t tensor = static_cast<tensor_t>(palign<matrix_t>(A, _shape[0]));
				matrix_t matrix = static_cast<matrix_t>(palign<vector_t>(A, _shape[0] * _shape[1]));

				for (size_t i = 0; i < _shape[0]; ++i)
				{
					tensor[i] = &matrix[i * _shape[1]];
					for (size_t j = 0; j < _shape[1]; ++j)
						tensor[i][j] = &_buffer[(i * _shape[1] + j) * _shape[2]];
				}
				_ptr = std::shared_ptr<void>(tensor, [](void* ptr){
					pdelete(static_cast<tensor_t>(ptr)[0]); pdelete(ptr);});
			}
			else if constexpr (N == 4)
			{
				tetrad_t tetrad = static_cast<tetrad_t>(palign<tensor_t>(A, _shape[0]));
				tensor_t tensor = static_cast<tensor_t>(palign<matrix_t>(A, _shape[0] * _shape[1]));
				matrix_t matrix = static_cast<matrix_t>(palign<vector_t>(A, _shape[0] * _shape[1] * _shape[2]));

				for (size_t t = 0; t < _shape[0]; ++t)
				{
					tetrad[t] = &tensor[t * _shape[1]];
					for (size_t i = 0; i < _shape[1]; ++i)
					{
						tetrad[t][i] = &matrix[(t * _shape[1] + i) * _shape[2]];
						for (size_t j = 0; j < _shape[2]; ++j)
							tetrad[t][i][j] = &_buffer[((t * _shape[1] + i) * _shape[2] + j) * _shape[3]];
					}
				}

				_ptr = std::shared_ptr<void>(tetrad, [](void* ptr){
					auto p = static_cast<tetrad_t>(ptr);
					pdelete(p[0][0]);
					pdelete(p[0]);
					pdelete(ptr);
				});
			}
		}

	public:
	/**
	 * @brief A rank-3 or rank-4 array that stores only the rows of a box (i, j) in [i0, i1] x [j0, j1]
	 *        dimensions N-3 and N-2, every k
	 *
	 * Indexing covers the whole shape as for a full array; the rows outside the box all point to one
	 * shared row of zeros, so reads there give 0. Writing there would change that row for every
	 * outside cell; padding_zero() tells whether it happened. Rank 4: the same box in every slice t.
	 */
	static carray boxed(std::array<size_t, N> shape, size_t i0, size_t i1, size_t j0, size_t j1)
	{
		static_assert(N == 3 || N == 4, "boxed arrays are rank 3 or 4");
		carray a;
		a._shape = shape_t(new size_t[N]);
		for (size_t d = 0; d < N; ++d)
			a._shape[d] = shape[d];
		const size_t nt = N == 4 ? shape[0] : 1;
		const size_t ni = shape[N - 3];
		const size_t nj = shape[N - 2];
		const size_t nk = shape[N - 1];
		const size_t bi = i1 - i0 + 1;
		const size_t bj = j1 - j0 + 1;
		a._count = nt * bi * bj * nk + nk; // the box rows, then the zero row
		a._buffer = make_shared_aarray<T>(A, a._count);
		fdtd::core::memory::zero_fill(a._buffer.get(), a._count);
		a._boxed = true;
		T *zero_row = &a._buffer[a._count - nk];
		// row tables: one per, t, i in the box, then one zero plane shared by every i outside it
		matrix_t rows = static_cast<matrix_t>(palign<vector_t>(A, (nt * bi + 1) * nj));
		matrix_t zero_plane = &rows[nt * bi * nj];
		for (size_t j = 0; j < nj; ++j)
			zero_plane[j] = zero_row;
		tensor_t planes = static_cast<tensor_t>(palign<matrix_t>(A, nt * ni));
		for (size_t t = 0; t < nt; ++t)
			for (size_t i = 0; i < ni; ++i)
			{
				if (i < i0 || i > i1)
				{
					planes[t * ni + i] = zero_plane;
					continue;
				}
				matrix_t r = &rows[(t * bi + (i - i0)) * nj];
				for (size_t j = 0; j < nj; ++j)
					r[j] = (j < j0 || j > j1) ? zero_row : &a._buffer[((t * bi + (i - i0)) * bj + (j - j0)) * nk];
				planes[t * ni + i] = r;
			}
		if constexpr (N == 3)
			a._ptr = std::shared_ptr<void>(planes, [rows](void *ptr) { pdelete(rows); pdelete(ptr); });
		else
		{
			tetrad_t tetrad = static_cast<tetrad_t>(palign<tensor_t>(A, nt));
			for (size_t t = 0; t < nt; ++t)
				tetrad[t] = &planes[t * ni];
			a._ptr = std::shared_ptr<void>(tetrad, [rows, planes](void *ptr) { pdelete(rows); pdelete(planes); pdelete(ptr); });
		}
		return a;
	}

	/// a boxed array's shared zero row is still zero, always true for a full array
	bool padding_zero() const
	{
		if (!_boxed || !_buffer)
			return true;
		const size_t nk = _shape[N - 1];
		for (size_t k = _count - nk; k < _count; ++k)
			if (_buffer[k] != T{})
				return false;
		return true;
	}

	/// the stored elements, stored() of them, state archives
	T *storage() { return _buffer.get(); }

	/// elements actually stored, fewer than size() for a boxed array
	size_t stored() const { return _buffer ? _count : 0; }

	/**	\brief Buffer RO operator. Used for read-only access to memory. */
	template<typename... IJK>
	const T& operator()(IJK... ijk) const
	{
		static_assert(sizeof...(IJK) == N, "number of indices must match rank of array");

		size_t idx[] = { static_cast<size_t>(ijk)... };

		if constexpr (N == 1) return static_cast<vector_t>(_ptr.get())[idx[0]];
		else if constexpr (N == 2) return static_cast<matrix_t>(_ptr.get())[idx[0]][idx[1]];
		else if constexpr (N == 3) return static_cast<tensor_t>(_ptr.get())[idx[0]][idx[1]][idx[2]];
		else if constexpr (N == 4) return static_cast<tetrad_t>(_ptr.get())[idx[0]][idx[1]][idx[2]][idx[3]];
		else static_assert(N <= 4, "rank of array not supported, for rank > 4");
	}

	/**	\brief Buffer RW operator. Used for read-write access to memory. */
	template<typename... IJK>
	decltype(auto) operator[](IJK&&... ijk)
	{
		size_t idx[] = { static_cast<size_t>(ijk)... };

		if constexpr (N == 1) return static_cast<vector_t>(_ptr.get())[idx[0]];
		else if constexpr (N == 2) return static_cast<matrix_t>(_ptr.get())[idx[0]];
		else if constexpr (N == 3) return static_cast<tensor_t>(_ptr.get())[idx[0]];
		else if constexpr (N == 4) return static_cast<tetrad_t>(_ptr.get())[idx[0]];
	}

	/** \brief Hierarchical view type: T*, T**, T***, T**** for rank 1..4. */
	using view_t = std::conditional_t<N == 1, T*, std::conditional_t<N == 2, T**, std::conditional_t<N == 3, T***, T****>>>;

	/** \brief Release the memory, this reference to it; the array becomes empty. */
	void
	reset()
	{
		_ptr.reset();
		_buffer.reset();
		_shape.reset();
		_count = 0;
		_boxed = false;
	}

	/** \brief Raw hierarchical view for C APIs and legacy kernels, nullptr when empty.
	 *  The view does not own the memory: it is valid while a carray holds it. */
	operator view_t() const
	{
		return _ptr ? static_cast<view_t>(_ptr.get()) : nullptr;
	}

	/** \brief Total number of elements. */
	size_t
	size() const
	{
		if (!_shape)
			return 0;
		return std::accumulate(_shape.get(), _shape.get() + N, size_t{1}, std::multiplies<size_t>());
	}

	/** \brief Extent along dimension d. */
	size_t
	extent(size_t d) const
	{
		return _shape[d];
	}

	/** \brief Fill the whole buffer with value-initialized T, zeros for arithmetic types. */
	carray<T, N, A>&
	zeros()
	{
		fdtd::core::memory::zero_fill(begin(), size_t(end() - begin()));
		return *this;
	}

	/**	\brief Range begin operator. Beginning of contiguous memory block. */
	T*
	begin()
	{
		return _buffer.get();
	}

	/**	\brief Range end operator. End of contiguous memory block. */
	T*
	end()
	{
		return _buffer.get() + _count;
	}

	/** \brief View acquisition. Acquire the hierarchical view (_ptr) for the buffer*/
	decltype(auto)
	get()
	{
		if constexpr (N == 1) return static_cast<vector_t>(_ptr.get());
		else if constexpr (N == 2) return static_cast<matrix_t>(_ptr.get());
		else if constexpr (N == 3) return static_cast<tensor_t>(_ptr.get());
		else if constexpr (N == 4) return static_cast<tetrad_t>(_ptr.get());
	}
};

//Convenience alias for construction.
//Default alignment is suitable for cache efficiency of most architectures
//Specification of alignment per application is recommended.
template<typename T> using cvector = carray<T, 1, 64>;///< shortcut for array of Type T and Rank 1.
template<typename T> using cmatrix = carray<T, 2, 64>;///< shortcut for array of Type T and Rank 2.
template<typename T> using ctensor = carray<T, 3, 64>;///< shortcut for array of Type T and Rank 3.
template<typename T> using ctetrad = carray<T, 4, 64>;///< shortcut for array of Type T and Rank 4.

/**
 * \brief Drop-in replacement for the legacy MemAllocat3DArray(p, NX, NY, NZ):
 *        allocates (NX+1) x (NY+1) x (NZ+1), one ghost cell per axis, and zero-fills.
 *        `.get()` yields the T*** the legacy kernels expect; the memory is released
 *        when the last carray referencing it goes away.
 */
template<typename T>
ctensor<T> field3(int nx, int ny, int nz)
{
	ctensor<T> a(nx + 1, ny + 1, nz + 1);
	a.zeros();
	return a;
}

/**
 * \brief n contiguous fields3 in one row-major block: [n][NX+1][NY+1][NZ+1], zero-filled.
 *        `a[m]` is the T*** of field m; `&a[m][0][0][0]` its contiguous data.
 */
template<typename T>
ctetrad<T> field4(int n, int nx, int ny, int nz)
{
	ctetrad<T> a(n > 0 ? n : 1, nx + 1, ny + 1, nz + 1);
	a.zeros();
	return a;
}

/**
 * \brief Zero-filled vector of n elements, at least one, like `new T[0]`.
 */
template<typename T>
cvector<T> vector1(size_t n)
{
	cvector<T> v(n > 0 ? n : 1);
	v.zeros();
	return v;
}

#endif //__C_ARRAY_H__
