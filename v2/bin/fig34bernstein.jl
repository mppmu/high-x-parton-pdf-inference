#!/usr/bin/julia
using BAT, DensityInterface
using PartonDensity
using QCDNUM
using Plots, Random, Distributions, ValueShapes, ParallelProcessingTools
using StatsBase, LinearAlgebra
using SpecialFunctions, Printf
const sf = SpecialFunctions;
using DelimitedFiles
using LaTeXStrings
using HypothesisTests
using Statistics
using Measures
using ArgParse
import HDF5
include("priors.jl")
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
        "--parametrisation", "-p"
            help = "Parametrisation -- Dirichlet or Valence"
            arg_type = String
            default = "Dirichlet"
        "--pseudodata", "-d"
            help = "Input pseudodata -- file in the pseudodata directory w/o the extension"
            arg_type = String
            default = ""
        "--fitresults", "-f"
            help = "Input fitresults -- file in the pseudodata directory w/o the extension"
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
c1 = :teal
c2 = :royalblue4
c3 = :midnightblue
c4 = :grey
color_scheme = :viridis
font_family = "Computer Modern"
default(fontfamily = "Computer Modern")
Plots.scalefontsizes()
Plots.scalefontsizes(1.2);
# Results
seed=parsed_args["seed"]
println(seed)
seedtxt=string(seed)

#Sim data!!!
pdf_params, sim_data=pd_read_sim(string("pseudodata/", parsed_args["pseudodata"], ".h5"))

#Fit results!!!
samples_data = bat_read(string("fitresults/", parsed_args["fitresults"], ".h5")).result;


counts_obs_ep_data=sim_data["counts_obs_ep"]
counts_obs_em_data=sim_data["counts_obs_em"]

nbins = size(counts_obs_ep_data)[1]


prob_ep_gen = zeros(nbins)
prob_em_gen = zeros(nbins)

prob_ep_sim = zeros(nbins)
prob_em_sim = zeros(nbins)

prob_ep_data = zeros(nbins)
prob_em_data = zeros(nbins)

mode_pars_data = mode(samples_data)

qcdnum_grid = QCDNUM.GridParams(x_min=[1.0e-3, 1.0e-1, 5.0e-1], x_weights=[1, 2, 2], nx=100,qq_bounds=[1.0e2, 3.0e4], qq_weights=[1.0, 1.0], nq=50, spline_interp=3)
qcdnum_params = QCDNUM.EvolutionParams(order=2, α_S=0.118, q0=100.0, grid_params=qcdnum_grid,n_fixed_flav=5, iqc=1, iqb=1, iqt=1, weight_type=1);
splint_params = QCDNUM.SPLINTParams();
quark_coeffs = QuarkCoefficients()


Ns = 300000 # Number of samples from posterior
rn = MersenneTwister(seed);
sub_samples = BAT.bat_sample(rn, samples_data, BAT.OrderedResampling(nsamples=Ns)).result;

forward_model_init(qcdnum_params, splint_params)



# Use +2 to avoid lightest colors (not easy to see)
color_scheme = :viridis

c1 = :teal
c2 = :royalblue4
c3 = :midnightblue
c4 = :grey

Plots.scalefontsizes()
Plots.scalefontsizes(1.2);
alpha = 0.6
prior_alpha = 0.2;

# Get some prior samples for plotting
prior=get_priors(parsed_args)
prior_samples=bat_sample(prior).result;

xlims_initial_U = (-15, 5) # initial_U
#xlims_D_u = (-10.0, 20.0) # (0.29, 0.37)
xlims_D_u =  (0.15, 0.45)

intervals = [0.68, 0.95]
labels = [L"~~\mathrm{Posterior}~68~\%", L"~~\mathrm{Posterior}~95~\%"]
prior_labels = [L"~~\mathrm{Prior}~68~\%", L"~~\mathrm{Prior}~95~\%"]
colors = [c3, c1]
prior_colors = [:grey40, :grey50]


if parsed_args["parametrisation"] == "Dirichlet"
#weights = [30.0, 15.0, 12.0, 6.0, 3.6, 0.85, 0.85, 0.85, 0.85]
#θ = rand(rng, Dirichlet(weights))
θ_true = [ 0.228, 0.104, 0.249, 0.249, 0.104, 0.052, 0.010, 0.005, 0.0005]
θ_sum=sum(θ_true[1:9])
θ_true=θ_true/θ_sum
K_u_true=3.7
K_d_true=3.7
λ_g1_true=0.5
λ_g2_true=-0.5
K_g_true=5.0
λ_q_true=-0.5 
K_q_true=6.0
end
if parsed_args["parametrisation"] == "Valence"
#weights = [5.0, 5.0, 1.0, 1.0, 1.0, 0.5, 0.5]
λ_u_true = 0.64;
K_u_true = 3.38;
λ_d_true = 0.67;
K_d_true = 4.73;
#FIXME!!!
θ_true=[0.22, 0.10, 0.24, 0.24, 0.10,0.05, 0.01, 0.005, 0.0005]
θ_sum=sum(θ_true[1:9])
θ_true=θ_true/θ_sum
end
if parsed_args["parametrisation"] == "Bernstein"
θ_true=[.33, .13, .27, .17, .073, 0.014, 0.002, 0.000001, .003]
θ_sum=sum(θ_true[1:9])
θ_true=θ_true/θ_sum
bspoly_params_true = [1,4,0,4,0,5]
λ_g1_true=1.5;
λ_g2_true=-0.4;
K_g_true=6.0;
λ_q_true=-0.25;
K_q_true=5.0;
initial_U_true = [-8.];
initial_D_true = [15.0]
end

plot(framestyle=:axes, size=(500, 400), fontfamily=font_family, 
    layout=@layout([a b; c{0.55w, 0.6h} d]), grid=false
    , right_margin=1mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=0mm
)

# Joint posterior
plot!(prior_samples, (:(initial_U), :(θ[1])),
    subplot=3, 
    xlabel=L"A_3", ylabel=L"\Delta_u",
    seriestype=:smallest_intervals_contourf, smoothing=4, 
    marginalmode=false, intervals=intervals, fillcolors=reverse(prior_colors), linewidth=0, 
    alpha=prior_alpha
    , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , right_margin=-2mm
    , left_margin=6mm
    , top_margin=0mm
    , bottom_margin=5mm
        , ylims=(0.18,0.45)
)
plot!(samples_data, (:(initial_U), :(θ[1])), 
    subplot=3,     
    xlabel=L"A_3", ylabel=L"\Delta_u",
    seriestype=:smallest_intervals_contourf, smoothing=2, 
    marginalmode=false, intervals=intervals, fillcolors=reverse(colors), linewidth=0, alpha=alpha,
    xlims=xlims_initial_U  #, ylims=xlims_D_u
    , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , right_margin=-2mm
    , left_margin=6mm
    , top_margin=0mm
    , bottom_margin=5mm
    , yticks=(0.2:0.1:0.4,["0.2","0.3","0.4"])
    , ylims=(0.18,0.45)
)
plot!([initial_U_true[1]],[θ_true[1]],
    seriestype = :scatter, subplot = 3, color = "red", label = " Truth", legend = :topleft,lw=0, 
    foreground_color_legend=:transparent, background_color_legend=:transparent,  lc=:red, markerstrokecolor=:red, legendfontsize=14)



# initial_U marginal
plot!(prior_samples, :initial_U,
    subplot=1, 
    legend=false, marginalmode=false, 
    seriestype=:smallest_intervals, intervals=intervals, 
    colors=prior_colors,alpha=prior_alpha
    , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , right_margin=-2mm
    , left_margin=6mm
    , top_margin=0mm
    , bottom_margin=5mm
)
plot!(samples_data, :initial_U, 
    legend=false, xlabel="", ylabel=L"P(A_3)", subplot=1, 
    xlims=xlims_initial_U, ylims=(0, 0.4), seriestype=:smallest_intervals, 
    marginalmode=false, intervals=intervals, colors=colors, alpha=alpha
    , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , right_margin=-2mm
    , left_margin=6mm
    , top_margin=0mm
    , bottom_margin=5mm
    , yticks=(0.0:0.2:0.4,["0","0.2","0.4"])
)
vline!([initial_U_true[1]], color="red", label=" Truth", lw=0.5)

# Delta_u marginal
plot!(prior_samples, :(θ[1]),
    subplot=4,
    legend=false, marginalmode=false, 
    seriestype=:smallest_intervals, intervals=intervals,
    colors=prior_colors, alpha=prior_alpha, 
    orientation=:horizontal
    , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , right_margin=1mm
    , left_margin=5mm
    , top_margin=0mm
    , bottom_margin=5mm
        , ylims=(0.18,0.45)
)
plot!(samples_data, :(θ[1]), 
    subplot=4,
    legend=false, ylabel="", xlabel=L"P(\Delta_u)", 
    #ylims=xlims_D_u, 
    xlims=(0, 2*45), 
    seriestype=:smallest_intervals, intervals=intervals, marginalmode=false, 
    colors=colors, alpha=alpha, orientation=:horizontal
    , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , xticks=(0:20:80,["0","20","40","60","80"])
    , right_margin=1mm
    , left_margin=5mm
    , top_margin=0mm
    , bottom_margin=5mm
    , yticks=(0.2:0.1:0.4,["0.2","0.3","0.4"])
        , ylims=(0.18,0.45)
)
hline!([θ_true[1]], color="red", label=" Truth", subplot=4, lw=0.5)

# Legend
plot!(prior_samples, (:(initial_U), :(θ[1])), 
    subplot=2,
    xlabel=L"A_3", ylabel=L"\Delta_u",
    seriestype=:smallest_intervals,
    marginalmode=false, intervals=intervals, interval_labels=prior_labels, 
    colors=reverse(prior_colors), linewidth=0, 
    alpha=prior_alpha+0.2
    , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , right_margin=1mm
    , left_margin=5mm
    , top_margin=0mm
    , bottom_margin=5mm
    , yticks=(0.2:0.1:0.4,["0.2","0.3","0.4"])
)
p = plot!(samples_data, (:(initial_U), :(θ[1])),
    subplot=2,
    seriestype=:smallest_intervals,
    marginalmode=false, intervals=intervals, colors=reverse(colors), 
    interval_labels=labels,
    linewidth=0, alpha=alpha+0.2, legend=:bottomleft, foreground_color_legend=:transparent, background_color_legend=:transparent,
    framestyle=:none, xlims=(0, 1.), ylims=(0, 0.1)
    , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=16
    , right_margin=1mm
    , left_margin=5mm
    , top_margin=0mm
    , bottom_margin=5mm
    , yticks=(0.2:0.1:0.4,["0.2","0.3","0.4"])
)

filename = string("figures/fig4-",parsed_args["fitresults"], "_v2.pdf")
savefig(p, filename)

end

main()
