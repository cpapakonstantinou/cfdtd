# Style Guide


### Code Editor

You can use any editor you like, however it is recommended to use Microsoft Visual Code which is available for Linux, macOS and Windows [here](https://code.visualstudio.com/).

If you chose to use a different editor, **make sure you have the following settings**:

- Indentation by tabs, not spaces
- Allman style braces

Allman style means
- The opening `{` should be on the next line after a statement, not on the same line.
- It should match indentation level of the opening statement
- `if` and `while` statements should be followed by a space before the `()`
 
Example:
```
if (a == 1)
{
    something();
}
```


### Indentation

- All indentation should be using 'tabs', not spaces
- One tab per indent
- Single-condition if statements should still have newlines, indented onto the next


### Line Length

- Limit lines to 79 characters


### Variable and function names

- Names should always be as descriptive as possible. Avoid abbreviating things. Long names are not a bad things since most editors will autocomplete anyway and it make the code much more readable. The only exception is for index variables like `i`, `j`, `k`. (Possibly symbols in equations too?)
- Variables (global, local, struct members) are snake_case: `num_cells_this_process_x`.
- Private members start with an underscore: `_rank`, not `rank_`.
- One space between the type and the name, no padding to align declarations: `int x;`, not `int    x;`.
- One declaration per line: `double x = 1;` and `double y = 2;` on two lines, not `double x = 1, y = 2;`.
- Compile-time constants and macros stay UPPERCASE: `C0`, `PI2`, `NRMAXANFREQUENCIES`.
- Field components keep their conventional capitals (`Ex`, `Hy`); everything else, including the grid and time steps (`dx`, `dt`), is snake_case.
- Functions and methods are snake_case: `load_sim_config`, `Driver::prepare_run`. Types are PascalCase: `SimConfig`, `ChiTensor`. Before naming a function, check that no variable, parameter or C/C++ library function in scope has that name.
- Module functions are named by what they do (`write_exact_multipoles`, `sample_nonlinear_polarisation`), not by the feature. A module variable gets a name a local cannot shadow (`analysis_omega`, not `omega`).
- Configuration keys are snake_case and named for the thing they set, in the section of the option they belong to (`[round_tripole] top_gap`, not a generic `length_2`).


### Files, modules and namespaces

- Each module is a header in `inc/<module>/` that documents the interface and an implementation in `src/<module>/` (e.g. `inc/grid/subdomain.h`, `src/grid/subdomain.cc`). A module owns its state.
- Everything is in the project namespace `fdtd`, each module in its own sub-namespace (`fdtd::grid`, `fdtd::materials`, `fdtd::config`). The modules and their order are on the [Architecture](Architecture.md) page.
- Qualify calls by their origin, so that names cannot collide:

  ```cpp
  std::signal(SIGTERM, on_signal); // C++ standard library, including the C headers it wraps (<csignal>, <cmath>, <cstdio>): std::
  MPI_Allreduce(...);              // third-party C libraries (MPI, HDF5): their own prefix, unqualified
  ::some_c_function();             // a plain C function with no std:: form and no library prefix: ::
  fdtd::some_fdtd_function();      // ours (or the module's own namespace)
  ```
  Include the `<cxxx>` form of C headers and call the `std::` names; `::` is then rarely needed.
- Never `using namespace std;`, in headers or implementation files: write `std::string`, `std::vector`, `std::abs`.
- No type aliases for standard or built-in types: `double`, `std::complex<double>`, `MPI_DOUBLE` are written as they are. Field precision is `double`; single precision or complex fields will come back as a template parameter of the field and update types, not as a global typedef.
- Physical constants are in `inc/core/constants.h` (`fdtd::constants`), the material slots in `inc/materials/slots.h` (`fdtd::materials`).
- No `using namespace` in headers. An implementation file names what it reads from another module with using-declarations (`using fdtd::fields::Ex;`).

### Compiler hints

Use the attribute syntax directly; no macros around it.
- `__attribute__((always_inline)) inline` on per-cell kernels that must fold into their loop; check with a profile (a kernel that shows up with one call per cell is not inlined).
- `[[assume(expr)]]` (C++23, GCC >= 13) for invariants the optimiser may use. Undefined behaviour when false, so only for what the code guarantees; precede it with an `assert(expr)` so debug builds check it. Measure each one: the field loops are memory-bound.
- `[[likely]]` / `[[unlikely]]` on branches whose frequency is known (CPML and TF/SF cells inside a row are rare).

### Math

- Enforce spaces around binary operators (+, -, *, /, %)
- Break long equations onto multiple lines using parentheses


### Function Length

- Not yet defined, aim for minimal (<40 lines, fits on 1 normal 1080p screen.)


### Commenting

- Use Doxygen style of commenting for Docstrings (See Doxygen template Below)
- Minimize inline/block comments within the code by using clear variable/function naming
- No banners: no comment or output line built from runs of `----` or `====`. A section gets a plain one-line comment.
- No asides in parentheses in comments: write `// cell size in m`, not `// cell size (m)`, and `// the halo, one cell per face`, not `// the halo (one cell per face)`. Parentheses stay where they are notation: formulas, index tuples `(i, j, k)`, half-open ranges `[k0, k1)`, calls `apply()`.

**Doxygen Example**
```
/**
* ... Text ...
*/
```


### Conditional Statements

- Be sure to use some form of 'else' or general error catching case.
- Use switch cases when applicable for optimal performance

## Callbacks over per-case copies

Code that differs only in which arrays, component or model it touches is written once and takes the varying part as a callback or a small bundle of arrays, not copied per case. Examples: `stretch::Drude2CP` serves Ex, Ey and Ez (the curl is a callback, the arrays an `AdeRow`), and `plane_dft` serves twelve former plane kernels. When porting legacy code, look for X/Y/Z or linear/nonlinear copies and fold them.

Per-component code is a template on the axis `A` (0 x, 1 y, 2 z). The arrays come from a component bundle (`materials::Component<A>`, built with `pick<A>(x, y, z)`), and whatever involves the other components is written with the cyclic axes `B1 = (A + 1) % 3`, `B2 = (A + 2) % 3` (e.g. the magnetisation cross term `M_{B1} dE_{B2} - M_{B2} dE_{B1}`). Material models are policy types with one `template <int A, class Boundary> static void update(const Cell &, Boundary &&)` (`inc/materials/models/cell.h`); a model is chosen once per stretch through `materials::with_model`, never per cell.
