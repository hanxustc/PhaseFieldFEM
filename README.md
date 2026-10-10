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

## Build

```bash
conda activate moose
make -j 28  # -> quad_test-opt
```

## Example

The initial orientation field is read from `poly.data`.

```bash
cd examples
mpiexec -n 48 ../quad_test-opt -i test.i
```

## Automatic differentiation

The AD and non-AD implementations are selected in the input file by the object types:

| | AD | non-AD |
|---|---|---|
| `[Materials]` | `ADFunctions` | `Functions` |
| `[Kernels]` | `ADTimeDerivative`, `ADGrad`, `ADBulk` | `TimeDerivative`, `Grad`, `Bulk` |
| `[Postprocessors]` | `ADElementIntegralMaterialProperty` | `ElementIntegralMaterialProperty` |

## License

LGPL-2.1, see `LICENSE`.
