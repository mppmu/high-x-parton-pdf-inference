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

# # A fit with BAT.jl
#
# In this example we show how to bring the PDF parametrisation and
# forward model together with `BAT.jl` to perform a fit of simulated data.
# This fit is a work in progress and just a starting point for verification
# of the method.


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
rng = MersenneTwister(seed)

if parsed_args["parametrisation"] == "Dirichlet"
weights = [20, 10, 20, 20, 5, 2.5, 1.5, 1.5, 0.5]
θ = rand(rng, Dirichlet(weights))
pdf_params = DirichletPDFParams(K_u=3.5, K_d=4.0, λ_g1=1.5, λ_g2=-0.4, K_g=6.0, λ_q=-0.25, K_q=5.0, θ=θ);
# replace weights with θ to specify the momentum fractions directly
@info "Valence λ:" pdf_params.λ_u pdf_params.λ_d
@info "Momenta:" pdf_params.θ[1],pdf_params.θ[2],pdf_params.θ[3],pdf_params.θ[4]
plot_input_pdfs(pdf_params)
end
if parsed_args["parametrisation"] == "Valence"
weights = [5.0, 5.0, 1.0, 1.0, 1.0, 0.5, 0.5]
λ_u = 0.64;
K_u = 3.38;
λ_d = 0.67;
K_d = 4.73;
θ = get_θ_val(rng, λ_u, K_u, λ_d, K_d, weights)
pdf_params = ValencePDFParams(λ_u=λ_u, K_u=K_u, λ_d=λ_d, K_d=K_d, λ_g1=0.50, λ_g2=-0.63, K_g=4.23, λ_q=-0.23, K_q=5.0, θ=θ);
# replace weights with θ to specify the momentum fractions directly
@info "Valence λ:" pdf_params.λ_u pdf_params.λ_d
@info "Momenta:" pdf_params.θ[1],pdf_params.θ[2],pdf_params.θ[3],pdf_params.θ[4]
plot_input_pdfs(pdf_params)
end


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

# ## Fit the data

#
# loosen the prior on K_g and negative lambdas
#


if parsed_args["parametrisation"] == "Dirichlet"
prior = NamedTupleDist(
    θ = Dirichlet([20, 10, 20, 20, 5, 2.5, 1.5, 1.5, 0.5]),
    K_u = Truncated(Normal(3.5, 0.5), 2., 5.),
    K_d = Truncated(Normal(3.5, 0.5), 2., 5.),
    λ_g1 = Uniform(0., 1.),
    λ_g2 = Uniform(-1.0, -0.1),
    K_g =  Truncated(Normal(4., 1.5), 2., 7.),
    λ_q = Uniform(-1.0, -0.1),
    K_q = Truncated(Normal(4., 1.5), 3., 10.),
    Beta1 =  Truncated(Normal(0, 1), -5, 5),
    Beta2 =  Truncated(Normal(0, 1), -5, 5),
    beta0_1=  Truncated(Normal(0, 1), -5, 5), 
    beta0_2=   Truncated(Normal(0, 1), -5, 5),    
    beta0_3= Truncated(Normal(0, 1), -5, 5), 
    beta0_4=  Truncated(Normal(0, 1), -5, 5), 
    beta0_5=  Truncated(Normal(0, 1), -5, 5), 
    beta0_6=  Truncated(Normal(0, 1), -5, 5), 
    beta0_7=  Truncated(Normal(0, 1), -5, 5), 
    beta0_8=   Truncated(Normal(0, 1), -5, 5)
);

if (parsed_args["priorshift"]==1)
    println("seting prior from Shifted Prior set ",seedtxt)

prior = NamedTupleDist(
    θ = Dirichlet([40, 10, 10, 10, 5, 2.5, 1.5, 1.5, 0.5]),
    K_u = Truncated(Normal(4.5, 0.5), 2, 5),
    K_d = Truncated(Normal(3.5, 0.5), 2., 5.),
    λ_g1 = Uniform(0., 1.),
    λ_g2 = Uniform(-1.0, -0.1),
    K_g =  Truncated(Normal(4., 1.5), 2., 7.),
    λ_q = Uniform(-1.0, -0.1),
    K_q = Truncated(Normal(4., 1.5), 3., 10.),
    Beta1 =  Truncated(Normal(0, 1), -5, 5),
    Beta2 =  Truncated(Normal(0, 1), -5, 5),
    beta0_1=  Truncated(Normal(0, 1), -5, 5), 
    beta0_2=   Truncated(Normal(0, 1), -5, 5),    
    beta0_3= Truncated(Normal(0, 1), -5, 5), 
    beta0_4=  Truncated(Normal(0, 1), -5, 5), 
    beta0_5=  Truncated(Normal(0, 1), -5, 5), 
    beta0_6=  Truncated(Normal(0, 1), -5, 5), 
    beta0_7=  Truncated(Normal(0, 1), -5, 5), 
    beta0_8=   Truncated(Normal(0, 1), -5, 5)
);
elseif (parsed_args["priorshift"]==2)
    println("seting prior from Shifted Prior set ",seedtxt)
prior = NamedTupleDist(
        θ = Dirichlet([20, 10, 30, 30, 5, 2.5, 1.5, 1.5, 0.5]),
        K_u = Truncated(Normal(2.5, 0.5), 2, 5),
    K_d = Truncated(Normal(3.5, 0.5), 2., 5.),
    λ_g1 = Uniform(0., 1.),
    λ_g2 = Uniform(-1.0, -0.1),
    K_g =  Truncated(Normal(4., 1.5), 2., 7.),
    λ_q = Uniform(-1.0, -0.1),
    K_q = Truncated(Normal(4., 1.5), 3., 10.),
    Beta1 =  Truncated(Normal(0, 1), -5, 5),
    Beta2 =  Truncated(Normal(0, 1), -5, 5),
    beta0_1=  Truncated(Normal(0, 1), -5, 5), 
    beta0_2=   Truncated(Normal(0, 1), -5, 5),    
    beta0_3= Truncated(Normal(0, 1), -5, 5), 
    beta0_4=  Truncated(Normal(0, 1), -5, 5), 
    beta0_5=  Truncated(Normal(0, 1), -5, 5), 
    beta0_6=  Truncated(Normal(0, 1), -5, 5), 
    beta0_7=  Truncated(Normal(0, 1), -5, 5), 
    beta0_8=   Truncated(Normal(0, 1), -5, 5)
);
end





end
if parsed_args["parametrisation"] == "Valence"
##FIXME!!!
weights = [5.0, 5.0, 1.0, 1.0, 1.0, 0.5, 0.5]
prior = NamedTupleDist(
    θ_tmp=Dirichlet(weights),
    λ_u=Truncated(Normal(pdf_params.λ_u, 1), 0, 1),
    K_u=Truncated(Normal(pdf_params.K_u, 1), 2, 10),
    λ_d=Truncated(Normal(pdf_params.λ_d, 1), 0, 1),
    K_d=Truncated(Normal(pdf_params.K_d, 1), 2, 10),
    λ_g1=Truncated(Normal(pdf_params.λ_g1, 1), -1, 0),
    λ_g2=Truncated(Normal(pdf_params.λ_g2, 1), -1, 0),
    K_g=Truncated(Normal(pdf_params.K_g, 1), 2, 10),
    λ_q=Truncated(Normal(pdf_params.λ_q, 0.1), -1, 0),
    K_q=Truncated(Normal(pdf_params.K_q, 0.5), 3, 7),
    Beta1 =  Truncated(Normal(0, 1), -5, 5),
    Beta2 =  Truncated(Normal(0, 1), -5, 5),
    beta0_1=  Truncated(Normal(0, 1), -5, 5), 
    beta0_2=   Truncated(Normal(0, 1), -5, 5),    
    beta0_3= Truncated(Normal(0, 1), -5, 5), 
    beta0_4=  Truncated(Normal(0, 1), -5, 5), 
    beta0_5=  Truncated(Normal(0, 1), -5, 5), 
    beta0_6=  Truncated(Normal(0, 1), -5, 5), 
    beta0_7=  Truncated(Normal(0, 1), -5, 5), 
    beta0_8=   Truncated(Normal(0, 1), -5, 5)    
);
end

if parsed_args["parametrisation"] == "Bernstein"

#prior = NamedTupleDist(
#    θ = Dirichlet([28.0, 12.5, 20.0, 20.0, 10.0, 1.4, 0.2, 10.e-5, 0.3]),
#    initial_U = [Truncated(Normal(30., 15.), 0, 80)],
#    initial_D = [Uniform(0., 20.)],
#    U_weights = [Truncated(Normal(30., 15.), 0, 80)],
#    D_weights = [Uniform(0., 20.)],

#    λ_g1 = Uniform(1., 2.0),
#    λ_g2 = Uniform(-0.5, -0.3),
#    K_g =  Uniform(5.,9.),
#    λ_q = Uniform(-0.5, -0.),
#    K_q = Uniform(3., 7.),
#    bspoly_params = [[0,3],[0,4],[1,4]],
#    )

prior = NamedTupleDist(
    θ = Dirichlet([34.0, 17.0, 22.5, 17.0, 7.3, 1.4, 0.2, 10.e-5, 0.3]),
    initial_U = Uniform(-10., 1.),
    initial_D = Uniform(10., 30.),
    λ_g1 = Uniform(3., 4.5),
    λ_g2 = Uniform(-1, -0.5),
    K_g =  Uniform(5.,9.),
    λ_q = Uniform(-1, -0.5),
    K_q = Uniform(3., 7.),
    #bspoly_params = [[0, 3], [0, 4], [1, 4], [0, 5]],
        bspoly_params = [1,4,0,4,0,5],
        Beta1 =  Truncated(Normal(0, 1), -5, 5),
    Beta2 =  Truncated(Normal(0, 1), -5, 5),
    beta0_1=  Truncated(Normal(0, 1), -5, 5), 
    beta0_2=   Truncated(Normal(0, 1), -5, 5),    
    beta0_3= Truncated(Normal(0, 1), -5, 5), 
    beta0_4=  Truncated(Normal(0, 1), -5, 5), 
    beta0_5=  Truncated(Normal(0, 1), -5, 5), 
    beta0_6=  Truncated(Normal(0, 1), -5, 5), 
    beta0_7=  Truncated(Normal(0, 1), -5, 5), 
    beta0_8=   Truncated(Normal(0, 1), -5, 5)
    )
end




# The `@critical` macro is used because `forward_model()` is currently not thread safe, so
# this protects it from being run in parallel.

likelihood = let d = sim_data
    counts_obs_ep = d["counts_obs_ep"]
    counts_obs_em = d["counts_obs_em"]
    nbins = d["nbins"]
    logfuncdensity(function (params)

       if parsed_args["parametrisation"] == "Bernstein"
       pdf_params = BernsteinDirichletPDFParams(

            initial_U=params.initial_U, 
            initial_D=params.initial_D, 
           # U_weights=params.U_weights,
           # D_weights=params.D_weights, 
            λ_g1=params.λ_g1, λ_g2=params.λ_g2, K_g=params.K_g, λ_q=params.λ_q, 
            K_q=params.K_q,
             bspoly_params = params.bspoly_params,
            
             θ=params.θ)
       end


       if parsed_args["parametrisation"] == "Dirichlet"
         pdf_params = DirichletPDFParams(K_u=params.K_u, K_d=params.K_d, λ_g1=params.λ_g1, λ_g2=params.λ_g2, K_g=params.K_g, λ_q=params.λ_q, K_q=params.K_q, θ=params.θ)
       end

       if parsed_args["parametrisation"] == "Valence"
         θ = get_scaled_θ(params.λ_u, params.K_u, params.λ_d,params.K_d, Vector(params.θ_tmp))
         pdf_params = ValencePDFParams(λ_u=params.λ_u, K_u=params.K_u, λ_d=params.λ_d, K_d=params.K_d, λ_g1=params.λ_g1, λ_g2=params.λ_g2, K_g=params.K_g, λ_q=params.λ_q, K_q=params.K_q, 
         θ=θ)
       end
          
       ParErrs = [params.beta0_1,params.beta0_2,params.beta0_3,params.beta0_4, params.beta0_5,params.beta0_6,params.beta0_7,params.beta0_8]
            
            counts_pred_ep, counts_pred_em = @critical  forward_model(pdf_params, qcdnum_params, splint_params, quark_coeffs
            ,ParErrs 
            );

            ll_value = 0.0
            for i in 1:nbins

                if counts_pred_ep[i] < 0
                   @debug "counts_pred_ep[i] < 0, setting to 0" i counts_pred_ep[i]
                   counts_pred_ep[i] = 0
                end

                if counts_pred_em[i] < 0
                   @debug "counts_pred_em[i] < 0, setting to 0" i counts_pred_em[i]
                   counts_pred_em[i] = 0
                end
                counts_pred_ep[i] =counts_pred_ep[i]*(1+0.018*params.Beta1)
                counts_pred_em[i] =counts_pred_em[i]*(1+0.018*params.Beta2)                

                ll_value += logpdf(Poisson(counts_pred_ep[i]), counts_obs_ep[i])
                ll_value += logpdf(Poisson(counts_pred_em[i]), counts_obs_em[i])
            end

            return ll_value
    end)
end

# We can now run the MCMC sampler. We will start by using the
# Metropolis-Hastings algorithm as implemented in `BAT.jl`.
# To demonstrate, we run a short chain of 10^3 iterations, which
# should take around 10 minutes.

posterior = PosteriorDensity(likelihood, prior);
mcalg = MetropolisHastings(proposal=BAT.MvTDistProposal(10.0))
convergence = BrooksGelmanConvergence(threshold=1.3);
burnin = MCMCMultiCycleBurnin(max_ncycles=50);

samples = bat_sample(posterior, MCMCSampling(mcalg=mcalg, nsteps=parsed_args["nsteps"], nchains=2)).result;

# Let's save the result for further analysis


bat_write(string("fitresults/fit-",parsed_args["parametrisation"],"-",parsed_args["priorshift"],"-",seedtxt,"-",parsed_args["pseudodata"],".h5"), samples)
end 

main()



