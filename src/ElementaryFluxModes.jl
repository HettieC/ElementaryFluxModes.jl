"""
Package `ElementaryFluxModes` provides a Julia implementation of the Double Description
method to calculate extreme rays of convex polyhedral cones. We follow the
method described in Terzer 2009 thesis for the polyhedral cone Ρ = {x ∈ ℜ^d | Ax = 0, x >= 0}.

The package can calculate elementary flux modes (EFMs) of optimal solutions to homogeneous enzyme-constrained genome scale metabolic models (ecGSMMs), and optimal flux modes (OFMs) of optimal solutions of inhomogeneous ecGSMMs.

It is also possible to use differentiation to calculate the sensitivity of the optimal composition of these EFMs or OFMs with respect to model parameters.
"""
module ElementaryFluxModes

using LinearAlgebra
using SparseArrays
using DocStringExtensions
using JuMP
using RowEchelon
using FastDifferentiation

include("utils.jl")
include("DDBinary.jl")
export DDBinary
export get_efms
export get_ofms
include("differentiate_efm.jl")
export differentiate_efm
export differentiate_ofm

end
