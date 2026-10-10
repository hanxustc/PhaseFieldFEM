# PhaseFieldFEM

Finite element implementation of the nonlocal orientation-field phase-field (NOPF) model for misorientation- and inclination-dependent grain boundaries, built as an application of the [MOOSE](https://mooseframework.inl.gov) framework.

The NOPF model is described in

> X. Han and A. van de Walle, *Nonlocal orientation field phase-field model for misorientation- and inclination-dependent grain boundaries*, Phys. Rev. B **113**, 184114 (2026).

## Contents

| Object | Type | Purpose |
|---|---|---|
| `initialTheta` | UserObject | stores the grain orientations of the initial condition |
| `qpNeighbors` | UserObject | builds the directional quadrature-point neighbor lists used by the nonlocal extrapolation |
| `ghostNeighbors` | RelationshipManager | ghosts the elements within the nonlocal search distance on each processor |
| `nonlocalTheta` | AuxKernel | computes the nonlocal orientation fields `thetaP` and `thetaM` |
| `RegGradAux` | AuxKernel | gradient components or magnitude of `theta` |
| `ADFunctions`, `Functions` | Material | free energy density and its derivatives (AD and non-AD) |
| `ADGrad`, `Grad` | Kernel | gradient term of the residual (AD and non-AD) |
| `ADBulk`, `Bulk` | Kernel | bulk term of the residual (AD and non-AD) |

The time derivative uses the standard MOOSE kernels `ADTimeDerivative` and `TimeDerivative`.

## Requirements

- MOOSE with the `phase_field` module
- PETSc with SuperLU_DIST for parallel direct solves (included in the MOOSE conda environment)

## Build

Place this folder next to a MOOSE checkout (`../moose`), or point to an existing one:

```bash
export MOOSE_DIR=/path/to/moose
make -j 8
```

This produces the executable `quad_test-opt`.

## Example

`examples/polycrystal_local_refinement` contains a polycrystalline grain growth simulation on a 64 x 64 base mesh that is locally refined at the grain boundaries.

```bash
cd examples/polycrystal_local_refinement
mpiexec -n 48 ../../quad_test-opt -i polycrystal_local_refinement.i
```

The initial orientation field is read from `poly.data`.

## Automatic differentiation

The AD and non-AD implementations are selected in the input file by the object types:

| | AD | non-AD |
|---|---|---|
| `[Materials]` | `ADFunctions` | `Functions` |
| `[Kernels]` | `ADTimeDerivative`, `ADGrad`, `ADBulk` | `TimeDerivative`, `Grad`, `Bulk` |
| `[Postprocessors]` | `ADElementIntegralMaterialProperty` | `ElementIntegralMaterialProperty` |

In the non-AD version the Jacobian is computed from the analytical derivatives of the residual.

## License

LGPL-2.1, see `LICENSE`.
