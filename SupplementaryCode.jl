

## Install packages, required on first start only

using Pkg
Pkg.add(["InformationGeometry", "Plots", "StaticArrays", "HypergeometricFunctions", "ModelingToolkit"])




# Documentation: https://rafaelarutjunjan.github.io/InformationGeometry.jl/dev/
using InformationGeometry
using Plots

## Simple toy model with non-linear parametrisation
@named ToyData = DataSet([1,2,3], [4,5,6.5], [0.5,0.45,0.6])
@named ToyModel = DataModel(ToyData, (x,p)->(p[1]+p[2])*x + exp(p[1]-p[2]))


# Show MLE
MLE(ToyModel)


# Compute 1σ and 2σ confidence regions
sols = ConfidenceRegions(ToyModel, 1:2)


# Plot confidence regions
VisualizeSols(ToyModel, sols)



# Show dataset with best fit
plot(ToyModel)
# Compute pointwise 1σ confidence bands using computed confidence boundary
ConfidenceBands(ToyModel, sols[1]; plot=true)
# Compute 2σ confidence bands
ConfidenceBands(ToyModel, sols[2]; plot=true)



## SCP data excerpted from http://supernova.lbl.gov/Union/figures/SCPUnion2.1_mu_vs_z.txt
include("SCPData.jl")

@named SCPData = DataSet(SCPredshifts, SCPdistances, SCPuncertainties;
                        xnames=["Redshift z"], ynames=["Distance Modulus μ(z)"])

const d_hubble = 299792458.0 / 70.0


using HypergeometricFunctions
CompactAntiDerivative(u::Real, Ω::Real, w::Real) = -2. * _₂F₁(0.5, -1/(6w), 1. -1/(6w), (1. - 1/Ω) * u^(3w)) / sqrt(Ω * u)
CompactAntiDerivativeSimp(Ω::Real, w::Real) = -2. * _₂F₁(0.5, -1/(6w), 1. -1/(6w), 1. - 1/Ω) / sqrt(Ω)
Integral(z::Real, θ::AbstractVector{<:Real}) = CompactAntiDerivative(z + 1., θ[1], θ[2]) - CompactAntiDerivativeSimp(θ[1], θ[2])

DistanceModulusFunc(z::Real, θ::AbstractVector{<:Real}) = 10. + 5. * log10(d_hubble * (1.0 + z) * Integral(z,θ))
DistanceModulus = ModelMap(DistanceModulusFunc, HyperCube([[0,1],[-10,0]]); pnames=["Matter Density Parameter Ωₘ₀", "Dark Energy E.o.S. ω₀"])


@named SCPModel = DataModel(SCPData, DistanceModulus, [0.2, -1])



plot(SCPModel; leg=:bottomright)


scpsols = ConfidenceRegions(SCPModel, 1:2)
# Visualize 1σ and 2σ confidence regions
VisualizeSols(SCPModel, scpsols)




## SIR model
using ModelingToolkit, StaticArrays
@parameters t β γ
@variables S(t) I(t) R(t);      Dt = Differential(t)
SIReqs = [Dt(S) ~ -β * I * S,    Dt(I) ~ +β * I * S - γ * I,   Dt(R) ~ +γ * I]

SIRstates = [S,I,R];    SIRparams = [β, γ];     SIRobservables = [2]
@named SIRsys = ODESystem(SIReqs, t, SIRstates, SIRparams)

### Influenza dataset from table 2 of https://www.researchgate.net/publication/336701551_On_parameter_estimation_approaches_for_predicting_disease_transmission_through_optimization_deep_learning_and_statistical_inference_methods
days = collect(1:14)
infected = [3, 8, 28, 75, 221, 291, 255, 235, 190, 126, 70, 28, 12, 5]
@named SIRData = DataSet(days, infected, 15ones(14); xnames=["Days"], ynames=["Infected"])

SIRinitial = @MVector [762, 1, 0.]

# SIR model with 2 parameters
SIRDM2 = DataModel(SIRData, SIRsys, SIRinitial, SIRobservables, [0.002, 0.5]; tol=1e-8)


SIR2sols = ConfidenceRegions(SIRDM2, 1:2; tol=1e-9)

VisualizeSols(SIRDM2, SIR2sols)


# 1σ confidence band for SIR predictions
plot(SIRDM2);   ConfidenceBands(SIRDM2, SIR2sols[1]; plot=true)



# SIR model with 3 parameters
SIRDM3 = DataModel(SIRData, SIRsys, p->(@MVector([763-p[1], p[1], 0.0]), p[2:3]), SIRobservables, [0.5, 0.002, 0.5];
                        pnames=["I₀", "β", "γ"], tol=1e-8)


SIR3sols = ConfidenceRegion(SIRDM3, 1; N=50, tol=1e-3)


VisualizeSols(SIRDM3, SIR3sols)


plot(SIRDM3);   ConfidenceBands(SIRDM3, SIR3sols; plot=true)
