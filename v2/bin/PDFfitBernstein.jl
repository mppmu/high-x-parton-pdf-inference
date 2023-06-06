#!/usr/bin/julia
using BAT, DensityInterface
using PartonDensity
using QCDNUM
using Plots, Random, Distributions, ValueShapes, ParallelProcessingTools
using StatsBase, LinearAlgebra
using DelimitedFiles

using ArgParse
import HDF5

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--priorshift"
            help = "Shift (variation) of priors. Only for Dirichlet"
            arg_type = Int
            default = 0
        "--seed", "-s"
            help = "Seed"
            arg_type = Int
            default = 42
        "--nsteps", "-n"
            help = "Number of steps"
            arg_type = Int
            default = 10^3
        "--parametrisation", "-p"
            help = "Parametrisation -- Dirichlet or Valence"
            arg_type = String
            default = "Dirichlet"
        "--pseudodata", "-d"
            help = "Input pseudodata -- file in the pseudodata directory w/o the extension"
            arg_type = String
            default = ""
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end


gr(fmt=:png);


counts_obs_ep = get_data_events(0)
counts_obs_em = get_data_events(1)

nbins = size(counts_obs_ep)[1]

#ENV["JULIA_DEBUG"] = "BAT"

# ### Specify the input PDFs
#
# See the *Input PDF parametrisation and priors* example for more information on the
# definition of the input PDFs.

seed=parsed_args["seed"]
println(seed)
seedtxt=string(seed)
Random.seed!(seed)

pdf_params = BernsteinDirichletPDFParams(initial_U = [-8.], initial_D = [0.5], λ_g1=1.5, λ_g2=-0.4, K_g=5.0,
                                λ_q=-0.25, θ = [0.36, .17, .27, .17, .073, 0.014, 0.002, 0.000001, .003],
                                bspoly_params = [[1,4],[0,4],[0,5]], K_q = 5.) #doesn't have relevance in this scenario

# ### Go from PDFs to counts in ZEUS detector bins
#
# Given the input PDFs, we can then evolve, calculate the cross sections, and fold through
# the ZEUS transfer matrix to get counts in bins. Here, we make use of some simple helper
# functions to do so. For more details, see the *Forward model* example.

# first specify QCDNUM inputs
qcdnum_grid = QCDNUM.GridParams(x_min=[1.0e-3, 1.0e-1, 5.0e-1], x_weights=[1, 2, 2], nx=100, qq_bounds=[1.0e2, 3.0e4], qq_weights=[1.0, 1.0], nq=50, spline_interp=3)
qcdnum_params = QCDNUM.EvolutionParams(order=2, α_S=0.118, q0=100.0, grid_params=qcdnum_grid, n_fixed_flav=5, iqc=1, iqb=1, iqt=1, weight_type=1);


# now SPLINT and quark coefficients
splint_params = QCDNUM.SPLINTParams();
quark_coeffs = QuarkCoefficients();

# initialise QCDNUM
#forward_model_init(qcdnum_grid, qcdnum_params, splint_params)
forward_model_init(qcdnum_params, splint_params)

# Plot the data
scatter(1:nbins, counts_obs_ep, label="Detected counts (eP)", color="blue")
scatter!(1:nbins, counts_obs_em, label="Detected counts (eM)", color="red")
plot!(xlabel="Bin number")

somepdf_params, sim_data = pd_read_sim(string("pseudodata/",parsed_args["pseudodata"],".h5"));


prior = NamedTupleDist(
    θ = Dirichlet([34.0, 17.0, 22.5, 17.0, 7.3, 1.4, 0.2, 10.e-5, 0.3]),
    initial_U = [Uniform(-10., 1.)],
    initial_D = [Uniform(10., 30.)],
    λ_g1 = Uniform(3., 4.5),
    λ_g2 = Uniform(-1, -0.5),
    K_g =  Uniform(5.,9.),
    λ_q = Uniform(-1, -0.5),
    K_q = Uniform(3., 7.),
    bspoly_params = [1,4,0,4,0,5],
)

posterior = PosteriorDensity(get_likelihood(pdf_params, sim_data, qcdnum_params, splint_params, quark_coeffs, true), 
                             prior);
convergence = BrooksGelmanConvergence(threshold=1.3);
burnin = MCMCMultiCycleBurnin(max_ncycles=100);

samples = bat_sample(posterior, MCMCSampling(mcalg=MetropolisHastings(), nsteps=parsed_args["nsteps"], nchains=2)).result

bat_write(string("fitresults/fit-",parsed_args["parametrisation"],"-",parsed_args["priorshift"],"-",seedtxt,"-",parsed_args["pseudodata"],".h5"), samples)
end 

main()
