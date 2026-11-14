// carray sanity checks, run under ASan: ghost-cell factory, zero fill, indexing, sharing.
#include "core/carray.h"
#include <complex>
#include <cstdio>
static int failures = 0;
#define CHECK(c) do { if (!(c)) { std::printf("FAIL %s:%d %s\n", __FILE__, __LINE__, #c); ++failures; } } while (0)
int main()
{
	auto f = field3<double>(3, 4, 5);              // legacy MemAllocat3DArray(p,3,4,5) shape
	CHECK(f.size() == 4u * 5u * 6u && f.extent(2) == 6);
	double*** p = f.get();
	bool allzero = true;
	for (int i = 0; i <= 3; ++i) for (int j = 0; j <= 4; ++j) for (int k = 0; k <= 5; ++k) allzero &= (p[i][j][k] == 0.0);
	CHECK(allzero);
	p[3][4][5] = 7.0; CHECK(f[3][4][5] == 7.0 && f(3, 4, 5) == 7.0 && *(f.end() - 1) == 7.0);
	auto alias = f;                                // shared ownership, same buffer
	alias[0][0][0] = 1.0; CHECK(f(0, 0, 0) == 1.0);
	cmatrix<std::complex<double>> m(3, 7); m.zeros(); m[2][6] = {1, 2}; CHECK(m(2, 6).imag() == 2);
	ctetrad<float> t(2, 3, 4, 5); t.zeros(); t[1][2][3][4] = 3.f; CHECK(t(1, 2, 3, 4) == 3.f && t.size() == 120);
	cvector<int> v(1); v.zeros(); CHECK(v[0] == 0);

	// boxed: rows (i, j) in [1, 2] x [2, 3] stored, every other row reads 0 from one shared zero row
	auto b = ctensor<double>::boxed({4, 5, 6}, 1, 2, 2, 3);
	CHECK(b.size() == 4u * 5u * 6u && b.stored() == 2u * 2u * 6u + 6u && b.padding_zero());
	bool bz = true;
	for (int i = 0; i < 4; ++i) for (int j = 0; j < 5; ++j) for (int k = 0; k < 6; ++k) bz &= (b[i][j][k] == 0.0);
	CHECK(bz);
	b[2][3][5] = 4.0; b[1][2][0] = 2.0;
	CHECK(b(2, 3, 5) == 4.0 && b(1, 2, 0) == 2.0 && b(0, 0, 5) == 0.0 && b(2, 4, 5) == 0.0 && b.padding_zero());
	CHECK(b[0][0] == b[3][4] && b[2][0] == b[0][0]); // the outside rows share one zero row
	b[0][0][1] = 1.0; CHECK(!b.padding_zero() && b(3, 4, 1) == 1.0); // a write outside shows everywhere outside
	auto bt = ctetrad<float>::boxed({2, 4, 5, 6}, 0, 0, 4, 4);
	bt[1][0][4][2] = 5.f; CHECK(bt(1, 0, 4, 2) == 5.f && bt(0, 0, 4, 2) == 0.f && bt(1, 1, 4, 2) == 0.f && bt.stored() == 2u * 6u + 6u && bt.padding_zero());
	std::printf(failures ? "%d FAILURE(S)\n" : "all carray tests passed\n", failures);
	return failures != 0;
}
