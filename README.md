# ElementaryFluxModes

[docs-img-stable]: https://img.shields.io/badge/docs-stable-blue
[docs-url-stable]: https://hettiec.github.io/ElementaryFluxModes.jl/stable

[docs-img-dev]: https://img.shields.io/badge/docs-dev-cyan
[docs-url-dev]: https://hettiec.github.io/ElementaryFluxModes.jl/dev

[ci-img]: https://github.com/HettieC/ElementaryFluxModes.jl/actions/workflows/ci.yml/badge.svg
[ci-url]: https://github.com/HettieC/ElementaryFluxModes.jl/actions/workflows/ci.yml

[cov-img]: https://codecov.io/github/HettieC/ElementaryFluxModes.jl/branch/main/graph/badge.svg?token=R9FSE2HZPU
[cov-url]: https://codecov.io/github/HettieC/ElementaryFluxModes.jl

| **Documentation** | **Tests** | **Coverage** |
|:---:|:---:|:---:|
| [![docs-img-stable]][docs-url-stable] [![docs-img-dev]][docs-url-dev] | [![CI status][ci-img]][ci-url] | [![codecov][cov-img]][cov-url] |

This package provides a Julia implementation of the double description algorithm for enumerating elementary flux modes, and can algebraically differentiate the optimal use of these modes to calculate sensitivities.

To use this package, [download and install Julia](https://julialang.org/downloads/), and add the package using the built in package manager
```julia
] add ElementaryFluxModes
```

We also provide the ability to calculate Optimal Flux Modes, analogues to EFMs for models with inhomogeneous flux constraints, such as an ATP maintenance reaction or a fixed glucose intake rate.

The documentation examples were all performed using [COBREXA.jl](https://github.com/LCSB-BioCore/COBREXA.jl), but you only need to input standard features of genome scale models (stoichiometric matrices, reaction indices of fixed fluxes, fixed flux rates), so any tool for solving genome scale models may be used. Any optimization solver that is compatible with [JuMP](https://jump.dev/) can be used to calculate EFM or OFM sensitivities. We have successfully used [Tulip.jl](https://github.com/ds4dm/Tulip.jl) and [HiGHS](https://github.com/jump-dev/HiGHS.jl), which are both open source.

You can test the installation through
```julia
] test ElementaryFluxModes
```

For more information, please see the documentation.

This package is maintained and open for extensions. Please feel free to discuss and suggest changes or ideas via pull requests.

#### Acknowledgements

`ElementaryFluxModes.jl` was developed at the Institute for Quantitative and Theoretical Biology at Heinrich Heine University Düsseldorf (qtb.hhu.de), and at the Luxembourg Centre for Systems Biomedicine of the University of Luxembourg (uni.lu/lcsb).

<img src="docs/src/assets/hhu.svg" alt="HHU logo" height="64px" style="height:64px; width:auto">   <img src="docs/src/assets/qtb.svg" alt="QTB logo" height="64px" style="height:64px; width:auto">   <img src="docs/src/assets/unilu.svg" alt="Uni.lu logo" height="64px">   <img src="docs/src/assets/lcsb.svg" alt="LCSB logo" height="64px">
