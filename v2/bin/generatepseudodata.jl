#!/usr/bin/julia
using BAT, DensityInterface
using PartonDensity
using QCDNUM
using Plots, Random, Distributions, ValueShapes, ParallelProcessingTools
using StatsBase, LinearAlgebra

using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--seed", "-s"
            help = "Seed"
            arg_type = Int
            default = 42
        "--parametrisation", "-p"
            help = "Parametrisation -- Dirichlet or Valence"
            arg_type = String
            default = "Dirichlet"
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
seed=parsed_args["seed"]
println(seed)
seedtxt=string(seed)
rng = MersenneTwister(seed)
Random.MersenneTwister(seed)

if parsed_args["parametrisation"] == "Dirichlet"

#weights = [30.0, 15.0, 12.0, 6.0, 3.6, 0.85, 0.85, 0.85, 0.85]
#θ = rand(rng, Dirichlet(weights))
θ = [ 0.228, 0.104, 0.249, 0.249, 0.104, 0.052, 0.010, 0.005, 0.0005]
pdf_params = DirichletPDFParams(K_u=3.7, K_d=3.7, λ_g1=0.5, λ_g2=-0.5, K_g=5.0,λ_q=-0.5, K_q=6.0, θ=θ);

end
if parsed_args["parametrisation"] == "Valence"
weights = [5.0, 5.0, 1.0, 1.0, 1.0, 0.5, 0.5]
λ_u = 0.64;
K_u = 3.38;
λ_d = 0.67;
K_d = 4.73;
θ = get_θ_val(rng, λ_u, K_u, λ_d, K_d, weights)
pdf_params = ValencePDFParams(λ_u=λ_u, K_u=K_u, λ_d=λ_d, K_d=K_d, λ_g1=0.50, λ_g2=-0.63, K_g=4.23, λ_q=-0.23, K_q=5.0, θ=θ);
end


if parsed_args["parametrisation"] == "Bernstein"

θ=[.33, .13, .27, .17, .073, 0.014, 0.002, 0.000001, .003]
bspoly_params = [[0,3],[0,4],[1,4]]
λ_g1=1.5;
λ_g2=-0.4;
K_g=6.0;
λ_q=-0.25;
K_q=5.0;
pdf_params = BernsteinDirichletPDFParams( λ_g1=λ_g2,
            λ_g2=λ_g2, K_g=K_g, λ_q=λ_q, K_q=K_q, 
            bspoly_params = bspoly_params, θ=θ)
#initial_U = [1.], initial_D = [1.], λ_g1=1.5, λ_g2=-0.4, K_g=6.0,
#                                λ_q=-0.25, θ = [.33, .13, .27, .17, .073, 0.014, 0.002, 0.000001, .003],
#                                bspoly_params = [[0,3],[0,4],[1,4]]
end

plot_input_pdfs(pdf_params)


qcdnum_grid = QCDNUM.GridParams(x_min=[1.0e-3, 1.0e-1, 5.0e-1], x_weights=[1, 2, 2], nx=100, qq_bounds=[1.0e2, 3.0e4], qq_weights=[1.0, 1.0], nq=50, spline_interp=3)
qcdnum_params = QCDNUM.EvolutionParams(order=2, α_S=0.118, q0=100.0, grid_params=qcdnum_grid,n_fixed_flav=5, iqc=1, iqb=1, iqt=1, weight_type=1);
    
splint_params = QCDNUM.SPLINTParams();
quark_coeffs = QuarkCoefficients();
forward_model_init(qcdnum_params, splint_params)



counts_pred_ep, counts_pred_em = forward_model(pdf_params, qcdnum_params,splint_params, quark_coeffs);
    
nbins = size(counts_pred_ep)[1]
counts_obs_ep = zeros(UInt64, nbins)
counts_obs_em = zeros(UInt64, nbins)

for i in 1:nbins
    counts_obs_ep[i] = rand(Poisson(counts_pred_ep[i]))
    counts_obs_em[i] = rand(Poisson(counts_pred_em[i]))
end


plot(1:nbins, counts_pred_ep, label="Expected counts (eP)", color="blue")
plot!(1:nbins, counts_pred_em, label="Expected counts (eM)", color="red")
scatter!(1:nbins, counts_obs_ep, label="Detected counts (eP)", color="blue")
scatter!(1:nbins, counts_obs_em, label="Detected counts (eM)", color="red")
plot!(xlabel="Bin number")

sim_data = Dict{String,Any}()
sim_data["nbins"] = nbins;
sim_data["counts_obs_ep"] = counts_obs_ep;
sim_data["counts_obs_em"] = counts_obs_em;

pd_write_sim(string("pseudodata/simulation-",parsed_args["parametrisation"],"-",seedtxt,".h5"), pdf_params, sim_data)
#QCDNUM.save_params("output/qcdnum_params_dir.h5", qcdnum_params)
#QCDNUM.save_params("output/splint_params_dir.h5", splint_params)
end

main()