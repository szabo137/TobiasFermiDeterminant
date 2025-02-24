# Tobias Fermi-Determinant

Calculation of the fermi determinant as given in `UWE.cpp`. The respective Julia code is
given in `Tobias.jl`.

## Usage

### Install Julia

You need a Julia compiler (see [here](https://github.com/JuliaLang/juliaup) for details),
or just run (on Linux)

```bash
curl -fsSL https://install.julialang.org | sh
```

This installs the latest Julia version and adds the path to the compiler to your `PATH`
variable.

### Install Dependencies

Run in the present directory:

```bash
julia --project -e "import Pkg; Pkg.instantiate()"
```

### Test Run

Run in the present directory:

```bash
julia --project -e """include(\"Tobias.jl\"); main()"""
```

### Benchmark Run

Run in the present directory:

```bash
julia --project -e """include(\"Tobias.jl\"); bench()"""
```
