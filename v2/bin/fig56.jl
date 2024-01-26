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



Ns = 10000 # Number of samples from posterior
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

        chisqep[s]+=(counts_ep_pred_s[j]-counts_ep_sampled[s, j])^2/counts_ep_pred_s[j]
        chisqem[s]+=(counts_em_pred_s[j]-counts_em_sampled[s, j])^2/counts_em_pred_s[j]

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

bins = 60

# Combine results for Δ_sea and Δ_g
# NB: for sea, our θ for quark/antiquark are already x2 - no need to do so here!
# in our x_q_x() in parametrisations_base.jl we take the input θ and halve it. 
comb_prior_samples = bat_transform(v -> (Δ_g = v.θ[3] + v.θ[4], Δ_sea = sum(v.θ[5:9])), prior_samples).result
comb_samples = bat_transform(v -> (Δ_g = v.θ[3] + v.θ[4], Δ_sea = sum(v.θ[5:9])), samples_data).result;

plot(framestyle=:axes, size=(PWIDTH, PWIDTH), fontfamily=font_family,  layout=@layout([a b; c d]), grid=false)
plot!(samples_data, (:(θ[1]), :(θ[2])), subplot=1, xlabel=L"\Delta_{u}", ylabel=L"\Delta_{d}",
    seriestype=:smallest_intervals_contourf, smoothing=4, 
    marginalmode=false, intervals=intervals, fillcolors=reverse(colors), linewidth=0, 
    alpha=prior_alpha
    , xtickfontsize=14,ytickfontsize=14,yguidefontsize=18,xguidefontsize=18,legendfontsize=14
    ,ylims=(0,0.32),yticks=(0:0.1:0.3,["0","0.1","0.2","0.3"])
    ,xlims=(0.1,0.42),xticks=(0.1:0.1:0.4,["0.1","0.2","0.3","0.4"])
    , right_margin=2mm
    , left_margin=8mm
    , top_margin=0mm
    , bottom_margin=5mm
)
plot!(prior_samples, (:(θ[1]), :(θ[2])), subplot=1, xlabel=L"\Delta_{u}", ylabel=L"\Delta_{d}",
   seriestype=:smallest_intervals_contourf, smoothing=4, 
    marginalmode=false, intervals=intervals, fillcolors=reverse(prior_colors), linewidth=0, 
    alpha=prior_alpha
, xtickfontsize=14,ytickfontsize=14,yguidefontsize=18,xguidefontsize=18,legendfontsize=14   
    
)
plot!([θ_true[1]],[θ_true[2]], subplot=1, color="red",seriestype=:scatter, label=" Truth", lw=0, 
foreground_color_legend=false,   lc=:red, markerstrokecolor=:red, legendfontsize=14)

##############################################3
comb_prior_samples = bat_transform(v -> (Δ_g = v.θ[3] + v.θ[4], Δ_u = v.θ[1]), prior_samples).result
comb_samples = bat_transform(v -> (Δ_g = v.θ[3] + v.θ[4], Δ_u = v.θ[1]), samples_data).result;

plot!(comb_samples, (:(Δ_u), :(Δ_g)),subplot=2,xlabel=L"\Delta_{u}", ylabel=L"\Delta_{g}",
   seriestype=:smallest_intervals_contourf, smoothing=4, 
    marginalmode=false, intervals=intervals, fillcolors=reverse(colors), linewidth=0, 
    alpha=prior_alpha
, xtickfontsize=14,ytickfontsize=14,yguidefontsize=18,xguidefontsize=18,legendfontsize=14
    ,ylims=(0.25,0.72),yticks=(0.3:0.1:0.72,["0.3","0.4","0.5","0.6","0.7"])
        ,xlims=(0.1,0.42),xticks=(0.1:0.1:0.4,["0.1","0.2","0.3","0.4"])
        , right_margin=2mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=5mm
    
)
plot!(comb_prior_samples, (:(Δ_u), :(Δ_g)), subplot=2,xlabel=L"\Delta_{u}", ylabel=L"\Delta_{g}",
   seriestype=:smallest_intervals_contourf, smoothing=4, 
    marginalmode=false, intervals=intervals, fillcolors=reverse(prior_colors), linewidth=0, 
    alpha=prior_alpha
, xtickfontsize=14,ytickfontsize=14,yguidefontsize=18,xguidefontsize=18,legendfontsize=14
)
θ_g_true=θ_true[3]+θ_true[4]
plot!([θ_true[1]],[θ_g_true], subplot=2,color="red",seriestype=:scatter, label=:none, lw=1, markerstrokecolor=:red)


##############################################3

comb_prior_samples = bat_transform(v -> (Δ_sea = sum(v.θ[5:9]), Δ_u = v.θ[1]), prior_samples).result
comb_samples = bat_transform(v -> (Δ_sea = sum(v.θ[5:9]), Δ_u = v.θ[1]), samples_data).result;

plot!(comb_samples, (:(Δ_u), :(Δ_sea)),subplot=3,xlabel=L"\Delta_{u}", ylabel=L"\Delta_\mathrm{sea}",
   seriestype=:smallest_intervals_contourf, smoothing=4, 
    marginalmode=false, intervals=intervals, fillcolors=reverse(colors), linewidth=0, 
    alpha=prior_alpha
, xtickfontsize=14,ytickfontsize=14,yguidefontsize=18,xguidefontsize=18,legendfontsize=14
    ,ylims=(0,0.32),yticks=(0:0.1:0.3,["0","0.1","0.2","0.3"])
        ,xlims=(0.1,0.42),xticks=(0.1:0.1:0.4,["0.1","0.2","0.3","0.4"]) 
        , right_margin=-2mm
    , left_margin=8mm
    , top_margin=0mm
    , bottom_margin=5mm
)
plot!(comb_prior_samples, (:(Δ_u), :(Δ_sea)), subplot=3,xlabel=L"\Delta_{u}", ylabel=L"\Delta_\mathrm{sea}",
   seriestype=:smallest_intervals_contourf, smoothing=4, 
    marginalmode=false, intervals=intervals, fillcolors=reverse(prior_colors), linewidth=0, 
    alpha=prior_alpha
, xtickfontsize=14,ytickfontsize=14,yguidefontsize=18,xguidefontsize=18,legendfontsize=14
)
θ_sea_true=sum(θ_true[5:9])
plot!([θ_true[1]],[θ_sea_true], subplot=3,color="red",seriestype=:scatter, label=:none, lw=2, markerstrokecolor=:red)


##############################################3

comb_prior_samples = bat_transform(v -> (Δ_sea = sum(v.θ[5:9]), Δ_g = v.θ[3] + v.θ[4]), prior_samples).result
comb_samples = bat_transform(v -> (Δ_sea = sum(v.θ[5:9]), Δ_g = v.θ[3] + v.θ[4]), samples_data).result;

plot!(comb_samples, (:(Δ_g), :(Δ_sea)),subplot=4, xlabel=L"\Delta_{g}", ylabel=L"\Delta_\mathrm{sea}",
   seriestype=:smallest_intervals_contourf, smoothing=4, 
    marginalmode=false, intervals=intervals, fillcolors=reverse(colors), linewidth=0, 
    alpha=prior_alpha
, xtickfontsize=14,ytickfontsize=14,yguidefontsize=18,xguidefontsize=18,legendfontsize=14
     ,ylims=(0,0.32),yticks=(0:0.1:0.3,["0","0.1","0.2","0.3"])
    ,xlims=(0.25,0.72),xticks=(0.3:0.1:0.72,["0.3","0.4","0.5","0.6","0.7"])
    , right_margin=2mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=5mm

)
plot!(comb_prior_samples, (:(Δ_g), :(Δ_sea)),subplot=4,xlabel=L"\Delta_{g}", ylabel=L"\Delta_\mathrm{sea}",
   seriestype=:smallest_intervals_contourf, smoothing=4, 
    marginalmode=false, intervals=intervals, fillcolors=reverse(prior_colors), linewidth=0, 
    alpha=prior_alpha
, xtickfontsize=14,ytickfontsize=14,yguidefontsize=18,xguidefontsize=18,legendfontsize=14
)
Δ_sea_true = sum(θ_true[5:9])
Δ_g_true = sum(θ_true[3:4])
plot!([Δ_g_true],[Δ_sea_true], color="red",subplot=4,seriestype=:scatter, label=:none, lw=2, markerstrokecolor=:red)


filename = string("figures/fig5-momentum-corr-", parsed_args["fitresults"], "_v2.pdf")
savefig(filename)
################################################################################################
function x_uv_x(x::Real, λ_u::Real, K_u::Real)
    A_u = 2 / sf.beta(λ_u, K_u + 1)
    return A_u * x^λ_u * (1 - x)^K_u
end

function x_dv_x(x::Real, λ_d::Real, K_d::Real)
    A_d = 1 / sf.beta(λ_d, K_d + 1)
    return A_d * x^λ_d * (1 - x)^K_d
end

function x_g_x(x::Real, λ_g1::Real, λ_g2::Real, K_g::Real,
    K_q::Real, w1::Real, w2::Real)
    A_g1 = w1 / sf.beta(λ_g1 + 1, K_g + 1)
    x_g1_x = A_g1 * x^λ_g1 * (1 - x)^K_g
    A_g2 = w2 / sf.beta(λ_g2 + 1, K_q + 1)
    x_g2_x = A_g2 * x^λ_g2 * (1 - x)^K_q
    return x_g1_x + x_g2_x
end

function x_q_x(x::Real, λ_q::Real, K_q::Real, w::Real)
    A_q = (w / 2) / sf.beta(λ_q + 1, K_q + 1)
    return A_q * x^λ_q * (1 - x)^K_q
end


function wrap_xuval(p::NamedTuple, x::Real)
    pdf_p = DirichletPDFParams(K_u=p.K_u, K_d=p.K_d, λ_g1=p.λ_g1, K_q=p.K_q,λ_g2=p.λ_g2, K_g=p.K_g, λ_q=p.λ_q, θ=p.θ)
    return PartonDensity.x_uv_x(x,pdf_p.λ_u, pdf_p.K_u)
end

function wrap_xdval(p::NamedTuple, x::Real)
    pdf_p = DirichletPDFParams(K_u=p.K_u, K_d=p.K_d, λ_g1=p.λ_g1, K_q=p.K_q,λ_g2=p.λ_g2, K_g=p.K_g, λ_q=p.λ_q, θ=p.θ)
    return PartonDensity.x_dv_x(x,pdf_p.λ_d, pdf_p.K_d)
end

function wrap_xg(p::NamedTuple, x::Real)
    scale = 1.0
    pdf_p = DirichletPDFParams(K_u=p.K_u, K_d=p.K_d, λ_g1=p.λ_g1, K_q=p.K_q,λ_g2=p.λ_g2, K_g=p.K_g, λ_q=p.λ_q, θ=p.θ)
    return PartonDensity.x_g_x(x, pdf_p.λ_g1, pdf_p.λ_g2, pdf_p.K_g, pdf_p.K_q, pdf_p.θ[3], pdf_p.θ[4]) * scale
end

function wrap_xsea(p::NamedTuple, x::Real)
    scale = 0.1
    pdf_p = DirichletPDFParams(K_u=p.K_u, K_d=p.K_d, λ_g1=p.λ_g1, K_q=p.K_q,λ_g2=p.λ_g2, K_g=p.K_g, λ_q=p.λ_q, θ=p.θ)
    x_ubar_x = PartonDensity.x_q_x(x, pdf_p.λ_q, pdf_p.K_q, pdf_p.θ[5])
    x_dbar_x = PartonDensity.x_q_x(x, pdf_p.λ_q, pdf_p.K_q, pdf_p.θ[6])
    x_sbar_x = PartonDensity.x_q_x(x, pdf_p.λ_q, pdf_p.K_q, pdf_p.θ[7])
    x_cbar_x = PartonDensity.x_q_x(x, pdf_p.λ_q, pdf_p.K_q, pdf_p.θ[8])
    x_bbar_x = PartonDensity.x_q_x(x, pdf_p.λ_q, pdf_p.K_q, pdf_p.θ[9])
    return 2 * (x_ubar_x + x_dbar_x + x_sbar_x + x_cbar_x + x_bbar_x) * scale
end
     


    x_grid = range(0.05, 1.0, length=95)

plot(framestyle=:axes, size=(PWIDTH,PWIDTH), fontfamily=font_family, 
    layout=@layout([a b; c d]), grid=false

    , right_margin=0mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=0mm

)


λ_u_true = θ_true[1]*(1+K_u_true)/(2-θ_true[1])
p = plot!(x_grid, [x_uv_x(x, λ_u_true, K_u_true) for x in x_grid],  lw=3, c=:red, 
    subplot=1
    , label=false
     , ylims=(0.00001, 65.0),yscale=:log10, xlims=(0.00,1.0)
    , foreground_color_legend=false
    , right_margin=0mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm
    , xtickfontsize=14,ytickfontsize=14,yguidefontsize=18,xguidefontsize=18,legendfontsize=14
        , xticks=(0.0:0.2:1.0,["0","0.2","0.4","0.6","0.8","1"]) 
        , yticks=([0.0001,0.001,0.01,0.1,1,10],[L"10^{-4}",L"10^{-3}",L"10^{-2}","0.1","1","10"]) 
    
)

Ns = 100 # Number of samples from posterior
sub_samples = BAT.bat_sample(samples_data, BAT.OrderedResampling(nsamples=Ns)).result;
for s in eachindex(sub_samples)

    pdf_params_s = DirichletPDFParams(K_u=sub_samples.v.K_u[s], K_d=sub_samples.v.K_d[s], K_q=sub_samples.v.K_q[s],
                                      λ_g1=sub_samples.v.λ_g1[s], 
                                      λ_g2=sub_samples.v.λ_g2[s],
                                      K_g=sub_samples.v.K_g[s], 
                                      λ_q=sub_samples.v.λ_q[s], 
                                      θ=Vector(sub_samples.v.θ[s]))
    λ_u= pdf_params_s.θ[1]*(1+pdf_params_s.K_u)/(2-pdf_params_s.θ[1])

    plot!(x_grid, [x_uv_x(x, λ_u, pdf_params_s.K_u) for x in x_grid],lw=0.5
    ,alpha=alpha_xf
    ,label=:none, subplot=1, ylims=(0.00001, 65.0),yscale=:log10)

end

    p = plot!(xlabel=L"x", ylabel=L"xu_v", subplot=1)


p = plot!(x_grid, [x_uv_x(x, λ_u_true, K_u_true) for x in x_grid], label=L"~xu_v \; \mathrm{true}", lw=3, c=:red, 
    subplot=1
     , ylims=(0.00001, 65.0),yscale=:log10, xlims=(0.00,1.0)
    , foreground_color_legend=false
    , right_margin=0mm
    , left_margin=5mm
    , top_margin=0mm
    , bottom_margin=-1mm
    , xtickfontsize=14,ytickfontsize=14,yguidefontsize=18,xguidefontsize=18,legendfontsize=14
        , xticks=(0.0:0.2:1.0,["0","0.2","0.4","0.6","0.8","1"]) 
        , yticks=([0.0001,0.001,0.01,0.1,1,10],[L"10^{-4}",L"10^{-3}",L"10^{-2}","0.1","1","10"])
)

#
λ_d_true = θ_true[2]*(1+K_d_true)/(1-θ_true[2])

p = plot!(x_grid, [x_dv_x(x, λ_d_true, K_d_true) for x in x_grid],   lw=3, c=:red,
    subplot=2
        , label=false
    , ylims=(0.00001, 65.0),yscale=:log10, xlims=(0.00,1.0)
    , foreground_color_legend=false
    , right_margin=1mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=5mm
    , xtickfontsize=14,ytickfontsize=14,yguidefontsize=18,xguidefontsize=18,legendfontsize=14
        , xticks=(0.0:0.2:1.0,["0","0.2","0.4","0.6","0.8","1"]) 
        , yticks=([0.0001,0.001,0.01,0.1,1,10],[L"10^{-4}",L"10^{-3}",L"10^{-2}","0.1","1","10"]) 
)


for s in eachindex(sub_samples)
    

    pdf_params_s = DirichletPDFParams(K_u=sub_samples.v.K_u[s], K_d=sub_samples.v.K_d[s], K_q=sub_samples.v.K_q[s],
                                      λ_g1=sub_samples.v.λ_g1[s], 
                                      λ_g2=sub_samples.v.λ_g2[s],
                                      K_g=sub_samples.v.K_g[s], 
                                      λ_q=sub_samples.v.λ_q[s], 
                                      θ=Vector(sub_samples.v.θ[s]))
    
        λ_d= pdf_params_s.θ[2]*(1+pdf_params_s.K_d)/(1-pdf_params_s.θ[2])

    plot!(x_grid, [x_dv_x(x, λ_d, pdf_params_s.K_d) for x in x_grid],ylims=(0.00001, 65.0),yscale=:log10,lw=0.5
    ,alpha=alpha_xf
        ,color=:gray
    ,label=:none, subplot=2)

end

    p = plot!(xlabel=L"x", ylabel=L"xd_v", subplot=2)

p = plot!(x_grid, [x_dv_x(x, λ_d_true, K_d_true) for x in x_grid], label=L"~xd_v \; \mathrm{true}",  lw=3, c=:red,
    subplot=2
    , ylims=(0.00001, 65.0),yscale=:log10, xlims=(0.00,1.0)
    , foreground_color_legend=false
    , right_margin=1mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=5mm
    , xtickfontsize=14,ytickfontsize=14,yguidefontsize=18,xguidefontsize=18,legendfontsize=14
        , xticks=(0.0:0.2:1.0,["0","0.2","0.4","0.6","0.8","1"]) 
        , yticks=([0.0001,0.001,0.01,0.1,1,10],[L"10^{-4}",L"10^{-3}",L"10^{-2}","0.1","1","10"]) 
)

#


    p = plot!(x_grid, [x_g_x(x, pdf_params.λ_g1, pdf_params.λ_g2, pdf_params.K_g, pdf_params.K_q, pdf_params.θ[3], pdf_params.θ[4])
                       for x in x_grid],   lw=3, c=:red,
    subplot=3
        , label=false
     , ylims=(0.00001, 65.0),yscale=:log10, xlims=(0.0,1.0)
    , foreground_color_legend=false
    , right_margin=-2mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=5mm
    , xtickfontsize=14,ytickfontsize=14,yguidefontsize=18,xguidefontsize=18,legendfontsize=14
        , xticks=(0.0:0.2:1.0,["0","0.2","0.4","0.6","0.8","1"]) 
        , yticks=([0.0001,0.001,0.01,0.1,1,10],[L"10^{-4}",L"10^{-3}",L"10^{-2}","0.1","1","10"]) 
)


for s in eachindex(sub_samples)
    
   pdf_params_s = DirichletPDFParams(K_u=sub_samples.v.K_u[s], K_d=sub_samples.v.K_d[s], K_q=sub_samples.v.K_q[s],
                                      λ_g1=sub_samples.v.λ_g1[s], 
                                      λ_g2=sub_samples.v.λ_g2[s],
                                      K_g=sub_samples.v.K_g[s], 
                                      λ_q=sub_samples.v.λ_q[s], 
                                      θ=Vector(sub_samples.v.θ[s]))

    p = plot!(x_grid, [x_g_x(x, pdf_params_s.λ_g1, pdf_params_s.λ_g2, pdf_params_s.K_g, pdf_params_s.K_q, pdf_params_s.θ[3], pdf_params_s.θ[4])
                       for x in x_grid], ylims=(0.00001, 65.0),yscale=:log10,lw=0.5
                       ,alpha=alpha_xf
                       ,label=:none,subplot=3)
end

    p = plot!(xlabel=L"x", ylabel=L"xg", subplot=3)

    p = plot!(x_grid, [x_g_x(x, pdf_params.λ_g1, pdf_params.λ_g2, pdf_params.K_g, pdf_params.K_q, pdf_params.θ[3], pdf_params.θ[4])
                       for x in x_grid], label=L"~xg \; \mathrm{true}",  lw=3, c=:red,
    subplot=3
     , ylims=(0.00001, 65.0),yscale=:log10, xlims=(0.0,1.0)
    , foreground_color_legend=false
    , right_margin=-2mm
    , left_margin=5mm
    , top_margin=0mm
    , bottom_margin=5mm
    , xtickfontsize=14,ytickfontsize=14,yguidefontsize=18,xguidefontsize=18,legendfontsize=14
        , xticks=(0.0:0.2:1.0,["0","0.2","0.4","0.6","0.8","1"]) 
        , yticks=([0.0001,0.001,0.01,0.1,1,10],[L"10^{-4}",L"10^{-3}",L"10^{-2}","0.1","1","10"]) 
)




#
    p = plot!(x_grid, [x_q_x(x, pdf_params.λ_q, pdf_params.K_q, pdf_params.θ[5]) for x in x_grid],
    ylims=(0.00001, 65.0),yscale=:log10, xlims=(0.0,1.0),   lw=3, c=:red,
    subplot=4
    , label=false
    , foreground_color_legend=false
    , right_margin=1mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=5mm
    , xtickfontsize=14,ytickfontsize=14,yguidefontsize=18,xguidefontsize=18,legendfontsize=14
    , xticks=(0.0:0.2:1.0,["0","0.2","0.4","0.6","0.8","1"]) 
        , yticks=([0.0001,0.001,0.01,0.1,1,10],[L"10^{-4}",L"10^{-3}",L"10^{-2}","0.1","1","10"]) 
)




for s in eachindex(sub_samples)
    
   pdf_params_s = DirichletPDFParams(K_u=sub_samples.v.K_u[s], K_d=sub_samples.v.K_d[s], K_q=sub_samples.v.K_q[s],
                                      λ_g1=sub_samples.v.λ_g1[s], 
                                      λ_g2=sub_samples.v.λ_g2[s],
                                      K_g=sub_samples.v.K_g[s], 
                                      λ_q=sub_samples.v.λ_q[s], 
                                      θ=Vector(sub_samples.v.θ[s]))


        p = plot!(x_grid, [x_q_x(x, pdf_params_s.λ_q, pdf_params_s.K_q, pdf_params_s.θ[5]) for x in x_grid],ylims=(0.00001, 65.0),yscale=:log10, lw=0.5
        ,alpha=alpha_xf
        ,label=:none,subplot=4)

end

    p = plot!(xlabel=L"x", ylabel=L"x\bar{u}", subplot=4)


    p = plot!(x_grid, [x_q_x(x, pdf_params.λ_q, pdf_params.K_q, pdf_params.θ[5]) for x in x_grid],
    ylims=(0.00001, 65.0),yscale=:log10, xlims=(0.0,1.0), label=L"~x\bar{u} \; \mathrm{true}",  lw=3, c=:red,
    subplot=4

    , foreground_color_legend=false
    , right_margin=1mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm
    , xtickfontsize=14,ytickfontsize=14,yguidefontsize=18,xguidefontsize=18,legendfontsize=14
    , xticks=(0.0:0.2:1.0,["0","0.2","0.4","0.6","0.8","1"]) 
        , yticks=([0.0001,0.001,0.01,0.1,1,10],[L"10^{-4}",L"10^{-3}",L"10^{-2}","0.1","1","10"]) 
)

filename = string("figures/fig6-parton-xf(x)-", parsed_args["fitresults"], "_v2.pdf")
savefig(filename)


end

main()

