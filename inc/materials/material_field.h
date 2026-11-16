/**
 * @file materials/material_field.h
 * @brief Per-material fields stored contiguously for the materials actually present.
 *
 * Dispersive polarizations and currents are indexed by material code (0 .. M-1) in the
 * kernels: P[m][i][j][k]. Only a few materials occur in a run, so only those are stored, one
 * ctensor [NX+1][NY+1][NZ+1] each, or only the rows of the material's cells, carray::boxed, and
 * `view` maps each material code to its field, nullptr for absent materials, as before.
 *
 *     P.require(m, nx, ny, nz);   // for every material present, any order, repeats ok
 *     P.commit();                 // one block per material, or commit(box_of): only the rows of a box
 *     P[m][i][j][k] ...           // kernels; P converts to T**** for functions taking it
 *     P.reset();                  // release
 */
#ifndef FDTD_MATERIAL_FIELD_H
#define FDTD_MATERIAL_FIELD_H

#include <algorithm>
#include <array>
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

#include "core/carray.h"

/// rows (i, j) of a material's cells: [i0, i1] x [j0, j1], inclusive
struct CellBox
{
	int i0;
	int i1;
	int j0;
	int j1;
};

template <typename T, std::size_t M>
class material_field
{
public:
	/// Mark material m as present with a field of (nx+1) x (ny+1) x (nz+1) cells.
	void require(int m, int nx, int ny, int nz)
	{
		if (std::find(_materials.begin(), _materials.end(), m) == _materials.end())
			_materials.push_back(m);
		_dims = {nx, ny, nz};
	}
	/// Allocate, zero-filled, the field of every required material, full size.
	void commit()
	{
		commit([](int) { return std::optional<CellBox>{}; });
	}
	/// As commit(), with the rows of material m limited to box_of(m) when it returns a box, carray::boxed:
	/// the other rows read 0.
	template <class BoxOf>
	void commit(BoxOf &&box_of)
	{
		_view.fill(nullptr);
		_store.clear();
		std::sort(_materials.begin(), _materials.end());
		for (int m : _materials)
		{
			const std::optional<CellBox> b = box_of(m);
			ctensor<T> a = b ? ctensor<T>::boxed({std::size_t(_dims[0] + 1), std::size_t(_dims[1] + 1), std::size_t(_dims[2] + 1)}, b->i0, b->i1, b->j0, b->j1)
							 : field3<T>(_dims[0], _dims[1], _dims[2]);
			_view[m] = a;
			_store.push_back(std::move(a));
		}
	}
	void reset()
	{
		_store.clear();
		_view.fill(nullptr);
		_materials.clear();
	}
	/// the rows outside every box still read 0, nothing wrote there
	bool padding_zero() const
	{
		return std::all_of(_store.begin(), _store.end(), [](const ctensor<T> &a) { return a.padding_zero(); });
	}
	/// f(array) for the field of every material present, in material order
	template <class F>
	void for_each(F &&f)
	{
		for (auto &a : _store)
			f(a);
	}
	/// elements allocated over all materials
	std::size_t stored() const
	{
		std::size_t n = 0;
		for (const auto &a : _store)
			n += a.stored();
		return n;
	}
	/// Field of material m, nullptr when m is absent.
	T ***&operator[](std::size_t m) { return _view[m]; }
	/// Material -> field table for functions taking T****.
	operator T ****() { return _view.data(); }

private:
	std::vector<ctensor<T>> _store;
	std::array<T ***, M> _view{};
	std::vector<int> _materials;
	std::array<int, 3> _dims{};
};
#endif // FDTD_MATERIAL_FIELD_H
