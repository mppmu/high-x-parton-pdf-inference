#!/usr/bin/julia
using BAT, DensityInterface
using PartonDensity
using QCDNUM
using Plots, Colors , Random, Distributions, ValueShapes, ParallelProcessingTools
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
include(string(dirname(pathof(PartonDensity)),"/../utils/priors.jl"))
include(string(dirname(pathof(PartonDensity)),"/../data/ZEUS_I1787035/ZEUS_I1787035.jl"))

#using bla
PWIDTH=1000
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
# Results
seed=parsed_args["seed"]
println(seed)
seedtxt=string(seed)

#Sim data!!!
pdf_params, sim_data, MD_TEMP=pd_read_sim(string("pseudodata/", parsed_args["pseudodata"], ".h5"),MD_G)

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

q2_edges_all = Any[]
x_edges_all = Any[]
for i in 1:nbins
    (q2_edges, x_edges) = ([MD_TEMP.m_q2bins_M_begin[i], MD_TEMP.m_q2bins_M_end[i]], [MD_TEMP.m_xbins_M_begin[i], MD_TEMP.m_xbins_M_end[i]])
    push!(q2_edges_all, q2_edges)
    push!(x_edges_all, x_edges)
end

q2_edges_unique = copy(q2_edges_all)
n_q2_bins = length(unique!(q2_edges_unique))

# get x, counts for each q2 range
counts_em_qsel = Any[]
counts_ep_qsel = Any[]
x_values = Any[]
for q2r in 1:n_q2_bins
    bin_sel = findall(==(q2_edges_unique[q2r]), q2_edges_all)
    push!(counts_em_qsel, counts_obs_em_data[bin_sel])
    push!(counts_ep_qsel, counts_obs_ep_data[bin_sel])
    push!(x_values, [mean(_) for _ in x_edges_all[bin_sel]])
end
# Use +2 to avoid lightest colors (not easy to see)



Ns = 1000 # Number of samples from posterior
rn = MersenneTwister(seed);
sub_samples = BAT.bat_sample(samples_data, BAT.OrderedResampling(nsamples=Ns)).result;

forward_model_init(qcdnum_params, splint_params)

counts_em_sampled = zeros(UInt64, (length(sub_samples), nbins))
counts_ep_sampled = zeros(UInt64, (length(sub_samples), nbins))
chisqep = zeros( length(sub_samples))
chisqem = zeros( length(sub_samples))


rng = MersenneTwister(seed);
nsyst=8
sys_err_params = rand(MvNormal(zeros(nsyst), zeros(nsyst)))

for s in eachindex(sub_samples)

    pdf_params_s = DirichletPDFParams(K_u=sub_samples.v.K_u[s], K_d=sub_samples.v.K_d[s], K_q=sub_samples.v.K_q[s],
                                      λ_g1=sub_samples.v.λ_g1[s], 
                                      λ_g2=sub_samples.v.λ_g2[s],
                                      K_g=sub_samples.v.K_g[s], 
                                      λ_q=sub_samples.v.λ_q[s], 
                                      θ=Vector(sub_samples.v.θ[s]))
        
    counts_ep_pred_s, counts_em_pred_s = forward_model(pdf_params_s,    qcdnum_params, splint_params,quark_coeffs, MD_TEMP, sys_err_params)
    
    for j in 1:nbins
        
        counts_ep_pred_s[j] *= 1 + 0.018 * sub_samples.v.Beta1[s]
        counts_em_pred_s[j] *= 1 + 0.018 * sub_samples.v.Beta2[s]
        
        counts_em_sampled[s, j] = rand(Poisson(counts_em_pred_s[j]))
        counts_ep_sampled[s, j] = rand(Poisson(counts_ep_pred_s[j]))
    end
    
end

color_scheme = :viridis
cmap = palette(color_scheme, n_q2_bins+2)

c1 = :teal
c2 = :royalblue4
c3 = :midnightblue
c4 = :grey
prior_alpha = 0.2;
alpha_xf=0.3;
colors = [c3, c1]
prior_colors = [:grey40, :grey50]


# Get some prior samples for plotting
prior=get_priors(parsed_args)
prior_samples=bat_sample(prior).result;

xlims_K_u = (2.0, 7.0) # (3.2, 4.4)
xlims_D_u = (0., 0.5) # (0.29, 0.37)
xlims_K_d = (2.0, 7.0) # (3.2, 4.4)
xlims_D_d = (0., 0.5) # (0.29, 0.37)
intervals = [0.68, 0.95]
labels = [L"~~\mathrm{Posterior}~68~\%", L"~~\mathrm{Posterior}~95~\%"]
prior_labels = [L"~~\mathrm{Prior}~68~\%", L"~~\mathrm{Prior}~95~\%"]

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




plot(framestyle=:axes, size=(PWIDTH,PWIDTH/2), fontfamily=font_family, 
    leftmargin=3Plots.mm, bottommargin=5Plots.mm, rightmargin=9mm,
    layout=@layout([a b c{0.20w}]),
    xtickfontsize=14,ytickfontsize=14,yguidefontsize=18,xguidefontsize=18
   , grid=false
)

plot!(inset=(1, bbox(0.23, 0.75, 0.55, 0.25, :bottom)))
plot!(inset=(2, bbox(0.23, 0.75, 0.55, 0.25, :bottom)))

plot!(xlim=(5e-2, 1), ylim=(0, 900), xlabel="\$x\$", ylabel="Counts", xscale=:log, 
    grid=false, legend=:false, foreground_color_legend=nothing,
    legendfontsize=10, thickness_scaling=1, 
    xticks=([0.1, 10^-0.5, 1.0],[L"$10^{-1}$",L"$10^{-0.5}$",L"$1$"]),
    subplot=1
    ,left_margin=13mm
    ,right_margin=1mm
    ,bottom_margin=8mm
    )
plot!(xlim=(5e-2, 1), ylim=(0, 900), xlabel="\$x\$", ylabel="", xscale=:log, 
    grid=false, legend=false, 
    xticks=([0.1, 10^-0.5, 1.0],[L"$10^{-1}$",L"$10^{-0.5}$",L"$1$"]),
    subplot=2
    ,left_margin=6mm
    ,right_margin=8mm
    ,bottom_margin=8mm
    )

for sp in [4, 5]
    plot!(xlim=(0.3, 1.0), ylim=(0, 50), grid=false, subplot=sp)
end

# Data
for i in 1:n_q2_bins
    label = @sprintf "  \$%g - %g\$" q2_edges_unique[i][1] q2_edges_unique[i][2]
    # Main plots and inset plots
    scatter!(x_values[i], counts_em_qsel[i], label="", color=cmap[i], markerstrokewidth=0, subplot=1)
    scatter!(x_values[i], counts_ep_qsel[i], label="", color=cmap[i], markerstrokewidth=0, subplot=2)
    scatter!(x_values[i], counts_em_qsel[i], label="", color=cmap[i], markerstrokewidth=0, subplot=4)
    scatter!(x_values[i], counts_ep_qsel[i], label="", color=cmap[i], markerstrokewidth=0, subplot=5)
    # For legend (invisible, trick to create space)
    scatter!(x_values[i], counts_em_qsel[i], label=label, color=cmap[i], markerstrokewidth=0, subplot=3)
end

# Samples
for s in eachindex(sub_samples)
    
    counts_em_s = Any[]
    counts_ep_s = Any[]
    for q2r in 1:n_q2_bins
        bin_sel = findall(==(q2_edges_unique[q2r]), q2_edges_all)
        push!(counts_em_s, counts_em_sampled[s, bin_sel])
        push!(counts_ep_s, counts_ep_sampled[s, bin_sel])
    end

    for i in 1:n_q2_bins
        scatter!(x_values[i], counts_em_s[i], color=cmap[i], markerstrokewidth=0, alpha=0.01, label="", subplot=1)
        scatter!(x_values[i], counts_ep_s[i], color=cmap[i], markerstrokewidth=0, alpha=0.01, label="", subplot=2)
        scatter!(x_values[i], counts_em_s[i], color=cmap[i], markerstrokewidth=0, alpha=0.01, label="", subplot=4)
        scatter!(x_values[i], counts_ep_s[i], color=cmap[i], markerstrokewidth=0, alpha=0.01, label="", subplot=5)
    end

end

# Draw box around zoomed region and connecting line
for sp in [1, 2]
    plot!([0.5, 0.5], [0, 50], color="red", linewidth=2, linestyle=:solid, subplot=sp)
    plot!([1.0, 1.0], [0, 50], color="red", linewidth=2, linestyle=:solid, subplot=sp)
    plot!([0.5, 1.0], [50, 50], color="red", linewidth=2, linestyle=:solid, subplot=sp)
    plot!([0.5, 1.0], [0, 0], color="red", linewidth=2, linestyle=:solid, subplot=sp)
end

plot!(legend=:left, foreground_color_legend=nothing, framestyle=:none,
    subplot=3, xlim=(1,2), ylim=(0, 900), legendfontsize=14, thickness_scaling=1,
    left_margin=-12mm, right_margin=3mm)
annotate!(1.4, 980*0.9, text("\$Q^2\$ [GeV\$^2\$]", 14, font_family), subplot=3)
annotate!(0.14, 600*0.9, text(L"$e^{-}$", 22, font_family), subplot=1)
p = annotate!(0.14, 600*0.9, text(L"$e^{+}$", 22, font_family), subplot=2)
p

savefig( string("figures/fig7-",parsed_args["fitresults"],"_v2.pdf"))
end

main()

