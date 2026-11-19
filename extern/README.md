# External modules

Copies of upstream modules, kept unmodified. The solver extends them in its own modules by
derivation; what a derived module has to work around is listed here with the upstream change
that would remove it. Edit the copies only by replacing them with a new upstream version.

| module | source | files | derived in |
|---|---|---|---|
| carray | github.com/cpapakonstantinou/carray, commit 678c610 | `carray.h`, `aligned_memory.h`, `LICENSE` | `inc/core/array.h`: `fdtd::core::array` |
| config | as received, 2026-11-03 | `config.h`, `config.cc` | `inc/io/ini.h`: `fdtd::io::Ini` |
| config | reconstructed | `exceptions.h` | |

`exceptions.h` is included by `config.h` but was not among the files received. The copy here
declares what `config.cc` throws, `cfg::file_error` and `cfg::syntax_error`; replace it with the
upstream file.

## carray

`fdtd::core::array` derives from `carray`. carray allocates full arrays and their row tables;
the derived class adds an empty state, boxed storage and the memory charge.

| derived array needs | workaround now | upstream change |
|---|---|---|
| an empty array, declared now and assigned later | the base holds a one-element placeholder | `carray()`, patch 0001 |
| the shape | a copy of the shape in the derived class | `size()`, `extent(d)`, patch 0001 |
| copies of const arrays | `const_cast` to reach the copy constructor | `carray(const carray&)`, patch 0001 |
| integer indices of any type | an integral-only `operator[]` | the sized constructor takes integral extents, patch 0001 |
| boxed storage, its own memory charge | own view, data and count; every accessor hidden | protected storage, patch 0002 |

Bugs, fixed by patch 0001:

- The element count is accumulated in `int`: an array of more than 2^31 elements is allocated
  too small. The solver's largest arrays stay below this, about 2·10^8 elements per rank at 1 nm.
- `aligned_alloc` is called with sizes that are not a multiple of the alignment. glibc accepts
  them; C11 leaves it undefined, macOS returns NULL, and valgrind and ASan report every call.

Patches: `0001-Empty-carray-size-extent-const-copy-size_t-element-c.patch`,
`0002-Protected-storage-so-derived-classes-can-lay-out-the.patch`, against 678c610.
With both applied, `array` keeps only the boxed layout, the charge and the field factories.

## config

`fdtd::io::Ini` derives from `cfg::Config`, which stores and queries the sections.

| Ini needs | workaround now | upstream change |
|---|---|---|
| parse text held in memory: the root rank reads the file, MPI broadcasts it | Ini parses the lines itself, with Config's `trim`, `split_key_value`, `add_key_val` | `parse(istream, name)`, `load_string`, patch 0001 |
| negative integers, such as `aperture_material = -1` | `Ini::to_integer`, `Ini::get_int` | `to_integer` without `to_unsigned`, patch 0001 |
| `true`/`false` booleans | `Ini::get_bool` | `get_bool` words, patch 0001 |

The solver's own dialect, kept in Ini also after the patch: case-insensitive names, `#` as well
as `;` comments, an error for a key before any section, and the sections in file order for the
geometry blocks, which `std::map` does not keep.

Patch: `0001-Parse-from-a-stream-or-a-string-signed-integers-bool.patch`, against the files as received.
