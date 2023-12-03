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

#c2=colorant"#CCE5E5"
#c1=colorant"#93A0AB"


c1 = :midnightblue
c2 = :teal
c3 = :grey
c4 = :grey
c5 = :grey

color_scheme = :viridis
font_family = "Computer Modern"
default(fontfamily = "Computer Modern")

alpha_posterior = 0.4;

# Results
seed=parsed_args["seed"]
println(seed)
seedtxt=string(seed)

#Sim data!!!
pdf_params, sim_data, MD_TEMP=pd_read_sim(string("pseudodata/", parsed_args["pseudodata"], ".h5"),MD_G)

#Fit results!!!
samples_data = bat_read(string("fitresults/", parsed_args["fitresults"], ".h5")).result;



samples_data1 = bat_transform(v -> (Δ_g = v.θ[3] + v.θ[4], Δ_sea = sum(v.θ[5:9]), K_u = v.K_u, K_q=v.K_q, K_g = v.K_g, K_d=v.K_d, 

#θ=v.θ 
beta0_1=v.beta0_1,
beta0_2=v.beta0_2,
beta0_3=v.beta0_3,
beta0_4=v.beta0_4,
beta0_5=v.beta0_5,
beta0_6=v.beta0_6,
beta0_7=v.beta0_7,
beta0_8=v.beta0_8


), samples_data).result



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


# As in PDF-fit-dirichlet.ipynb
qcdnum_grid = QCDNUM.GridParams(x_min=[1.0e-3, 1.0e-1, 5.0e-1], x_weights=[1, 2, 2], nx=100,qq_bounds=[1.0e2, 3.0e4], qq_weights=[1.0, 1.0], nq=50, spline_interp=3)
qcdnum_params = QCDNUM.EvolutionParams(order=2, α_S=0.118, q0=100.0, grid_params=qcdnum_grid,n_fixed_flav=5, iqc=1, iqb=1, iqt=1, weight_type=1);
splint_params = QCDNUM.SPLINTParams();
quark_coeffs = QuarkCoefficients()

q2_edges_all = Any[]
x_edges_all = Any[]
for i in 1:nbins
    #(q2_edges, x_edges) = get_bin_info(i, quiet=true);
    
    #(q2_edges, x_edges) = ([BinQ2low[n], BinQ2high[n]], [Binxlow[n], Binxhigh[n]])
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

Ns = 30000 # Number of samples from posterior
rn = MersenneTwister(seed);
sub_samples = BAT.bat_sample(samples_data, BAT.OrderedResampling(nsamples=Ns)).result;

forward_model_init(qcdnum_params, splint_params)

counts_em_sampled = zeros(UInt64, (length(sub_samples), nbins))
counts_ep_sampled = zeros(UInt64, (length(sub_samples), nbins))
chisqep = zeros( length(sub_samples))
chisqem = zeros( length(sub_samples))


rng = MersenneTwister(seed);
nsyst=8
sys_err_params = rand(rng, MvNormal(zeros(nsyst), zeros(nsyst)))



# Use +2 to avoid lightest colors (not easy to see)



cmap = palette(color_scheme, n_q2_bins+2)
alpha = 0.6
prior_alpha = 0.2;

# Get some prior samples for plotting
prior=get_priors(parsed_args)

prior_samples =bat_sample(prior).result;
prior_samples1 = bat_transform(v -> (Δ_g = v.θ[3] + v.θ[4], Δ_sea = sum(v.θ[5:9]), K_u = v.K_u, K_q=v.K_q, K_g = v.K_g, K_d=v.K_d ,  

#θ=v.θ
beta0_1=v.beta0_1,
beta0_2=v.beta0_2,
beta0_3=v.beta0_3,
beta0_4=v.beta0_4,
beta0_5=v.beta0_5,
beta0_6=v.beta0_6,
beta0_7=v.beta0_7,
beta0_8=v.beta0_8,
beta0_9=v.Beta1,
beta0_10=v.Beta2

), prior_samples).result



xlims_K_u = (2.0, 7.0) # (3.2, 4.4)
xlims_D_u = (0., 0.5) # (0.29, 0.37)
xlims_K_d = (2.0, 7.0) # (3.2, 4.4)
xlims_D_d = (0., 0.5) # (0.29, 0.37)
intervals = [0.68, 0.95]
labels = [L"~~\mathrm{Posterior}~68~\%", L"~~\mathrm{Posterior}~95~\%"]
prior_labels = [L"~~\mathrm{Prior}~68~\%", L"~~\mathrm{Prior}~95~\%"]
colors = [c3, c1]


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
weights = [5.0, 5.0, 1.0, 1.0, 1.0, 0.5, 0.5]
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
bspoly_params = [[0,3],[0,4],[1,4]]
λ_g1_true=1.5;
λ_g2_true=-0.4;
K_g_true=6.0;
λ_q_true=-0.25;
K_q_true=5.0;
U_weights_true = [1.]
D_weights_true = [1.]
end

NNN = 8


l = @layout [
        [a1{0.14w}  b1{0.12w} c1{0.12w} d1{0.12w} e1{0.12w} f1{0.12w} g1{0.12w} h1{0.12w}]
        [a2{0.14w} b2{0.12w} c2{0.12w} d2{0.12w} e2{0.12w} f2{0.12w} g2{0.12w} h2{0.12w}]
        [a3{0.14w} b3{0.12w} c3{0.12w} d3{0.12w} e3{0.12w} f3{0.12w} g3{0.12w} h3{0.12w} ]
        [a4{0.14w} b4{0.12w} c4{0.12w} d4{0.12w} e4{0.12w} f4{0.12w}  g4{0.12w} h4{0.12w}]
        [a5{0.14w} b5{0.12w} c5{0.12w} d5{0.12w} e5{0.12w} f5{0.12w} g5{0.12w} h5{0.12w}]
        [a6{0.14w} b6{0.12w} c6{0.12w} d6{0.12w} e6{0.12w} f6{0.12w} g6{0.12w} h6{0.12w}]
        [a7{0.14w} b7{0.12w} c7{0.12w} d7{0.12w} e7{0.12w} f7{0.12w} g7{0.12w} h7{0.12w}]
        [a8{0.14w} b8{0.12w} c8{0.12w} d8{0.12w} e8{0.12w} f8{0.12w} g8{0.12w} h8{0.12w}]

]





p=plot(size=(PWIDTH,PWIDTH),
layout=l,
    colors=[c1, c2, c3],
    frame=:box,
    legend = :none, 
    framestyle = :box,
    ylims=(1.9, 10.1),
    xlims=(1.9, 10.1), 
   # plot_titlevspan=0.001,
  #  xlabel=["" "" "" "" "" "" "" "" "" "" "" "" L"p(K_{u})" L"K_{d}" L"K_{q}" L"K_{g}"],
  #  ylabel=[L"p(K_{u})" "" "" "" L"K_{d}" "" "" "" L"K_{q}" "" "" "" L"K_{g}" "" "" ""],
    # xticks=(-3.0:1:3.0,["-3","-2","-1","0","1","2","3"]),    
  #  yticks=(0:2:12,["","","","","","",""])  ,
    xticks=:none,
    yticks=:none,
    top_margin=-3mm,
    bottom_margin=-1mm,
    right_margin=-1mm,
    left_margin=-1mm,
    xtickfontsize=10,
    ytickfontsize=10,
    yguidefontsize=12,
    xguidefontsize=12,
    fontfamily=font_family,
    grid=false
)
#plot!(samples_data, :(θ[1]),subplot=1, xlabel="",ylabel=L"θ[1]",colors=[c1, c2, c3],alpha=alpha_posterior, bins=100)
#for i in 2:NNN
#plot!(legend=false,grid=false,foreground_color_subplot=:white,subplot=i)
#end
j=1
plot!(samples_data, :(beta0_1) ,subplot=(j-1)*NNN+1, xlabel="",ylabel="",colors=[c1, c2, c3],alpha=alpha_posterior, bins=100)


j=2
plot!(samples_data, (:(beta0_1), :(beta0_2)),subplot=(j-1)*NNN+1, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data,           :(beta0_2) ,subplot=(j-1)*NNN+2, xlabel="",ylabel="",colors=[c1, c2, c3],alpha=alpha_posterior, bins=100)

j=3
plot!(samples_data, (:(beta0_1), :(beta0_3)),subplot=(j-1)*NNN+1, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data, (:(beta0_2), :(beta0_3)),subplot=(j-1)*NNN+2, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data,           :(beta0_3) ,subplot=(j-1)*NNN+3, xlabel="",ylabel="",colors=[c1, c2, c3],alpha=alpha_posterior, bins=100)


j=4
plot!(samples_data, (:(beta0_1), :(beta0_4)),subplot=(j-1)*NNN+1, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data, (:(beta0_2), :(beta0_4)),subplot=(j-1)*NNN+2, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data, (:(beta0_3), :(beta0_4)),subplot=(j-1)*NNN+3, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data,           :(beta0_4) ,subplot=(j-1)*NNN+4, xlabel="",ylabel="",colors=[c1, c2, c3],alpha=alpha_posterior, bins=100)


j=5
plot!(samples_data, (:(beta0_1), :(beta0_5)),subplot=(j-1)*NNN+1, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data, (:(beta0_2), :(beta0_5)),subplot=(j-1)*NNN+2, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data, (:(beta0_3), :(beta0_5)),subplot=(j-1)*NNN+3, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data, (:(beta0_4), :(beta0_5)),subplot=(j-1)*NNN+4, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data,           :(beta0_5) ,subplot=(j-1)*NNN+5, xlabel="",ylabel="",colors=[c1, c2, c3],alpha=alpha_posterior, bins=100)


j=6
plot!(samples_data, (:(beta0_1), :(beta0_6)),subplot=(j-1)*NNN+1, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data, (:(beta0_2), :(beta0_6)),subplot=(j-1)*NNN+2, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data, (:(beta0_3), :(beta0_6)),subplot=(j-1)*NNN+3, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data, (:(beta0_4), :(beta0_6)),subplot=(j-1)*NNN+4, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data, (:(beta0_5), :(beta0_6)),subplot=(j-1)*NNN+5, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data,           :(beta0_6) ,subplot=(j-1)*NNN+6, xlabel="",ylabel="",colors=[c1, c2, c3],alpha=alpha_posterior, bins=100)

j=7
plot!(samples_data, (:(beta0_1), :(beta0_7)),subplot=(j-1)*NNN+1, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data, (:(beta0_2), :(beta0_7)),subplot=(j-1)*NNN+2, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data, (:(beta0_3), :(beta0_7)),subplot=(j-1)*NNN+3, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data, (:(beta0_4), :(beta0_7)),subplot=(j-1)*NNN+4, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data, (:(beta0_5), :(beta0_7)),subplot=(j-1)*NNN+5, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data, (:(beta0_6), :(beta0_7)),subplot=(j-1)*NNN+6, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data,           :(beta0_7) ,subplot=(j-1)*NNN+7, xlabel="",ylabel="",colors=[c1, c2, c3],alpha=alpha_posterior, bins=100)

j=8
plot!(samples_data, (:(beta0_1), :(beta0_8)),subplot=(j-1)*NNN+1, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data, (:(beta0_2), :(beta0_8)),subplot=(j-1)*NNN+2, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data, (:(beta0_3), :(beta0_8)),subplot=(j-1)*NNN+3, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data, (:(beta0_4), :(beta0_8)),subplot=(j-1)*NNN+4, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data, (:(beta0_5), :(beta0_8)),subplot=(j-1)*NNN+5, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data, (:(beta0_6), :(beta0_8)),subplot=(j-1)*NNN+6, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data, (:(beta0_7), :(beta0_8)),subplot=(j-1)*NNN+7, xlabel="",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
plot!(samples_data,           :(beta0_8) ,subplot=(j-1)*NNN+8, xlabel="",ylabel="",colors=[c1, c2, c3],alpha=alpha_posterior, bins=100)




for j in 1:NNN-1
#for i in 1:j
#plot!(samples_data, (:(θ[i]), :(θ[j])),subplot=j*(NNN-1)+i, xlabel="",ylabel=L"K_d",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior,markerstrokewidth = 0, marginalmode=false)
#end
#xx=string("(θ[",j,"])")
#plot!(samples_data, :xx,subplot=j*(NNN-1)+j, xlabel="",ylabel="",colors=[c1, c2, c3],alpha=alpha_posterior, bins=100)
for i in j+1:NNN
plot!(legend=false,grid=false,foreground_color_subplot=:white,subplot=(j-1)*NNN+i)
end
end




#plot!(samples_data1, (:(θ[1]), :(θ[7])),subplot=7*NNN+1, xlabel=L"K_u",ylabel=L"\Delta_{\mathrm{sea}}",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior, marginalmode=false)
#plot!(samples_data1, (:(K_d), :(Δ_sea)),subplot=7*NNN+2, xlabel=L"K_d",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior, marginalmode=false)
#plot!(samples_data1, (:(K_q), :(Δ_sea)),subplot=7*NNN+3, xlabel=L"K_q",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior, marginalmode=false)
#plot!(samples_data1, (:(K_g), :(Δ_sea)),subplot=7*NNN+4, xlabel=L"K_g",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior, marginalmode=false)
#plot!(samples_data1, (:(Δ_g),:(Δ_sea)), subplot=7*NNN+5, xlabel=L"\Delta_g",ylabel="",fillcolors=reverse([c1, c2, c3]),alpha=alpha_posterior, marginalmode=false)
#plot!(samples_data1, :(Δ_sea),subplot=5*NNN+6, xlabel=L"\Delta_{\mathrm{sea}}",ylabel="",colors=[c1, c2, c3],alpha=alpha_posterior, bins=100)

for j in 2:NNN
plot!(p[(j-1)*NNN+1],left_margin=3mm,yticks=(-3.0:3:3.0,["-3","0","3"]))
end

for i in 1:NNN
plot!(p[(NNN-1)*NNN+i],xticks=(-3.0:3:3.0,["-3","0","3"]))
end

plot!(p[1],ylabel=L"\theta_1")
plot!(p[1+NNN],ylabel=L"\theta_2")
plot!(p[1+2*NNN],ylabel=L"\theta_3")
plot!(p[1+3*NNN],ylabel=L"\theta_4")
plot!(p[1+4*NNN],ylabel=L"\theta_5")
plot!(p[1+5*NNN],ylabel=L"\theta_6")
plot!(p[1+6*NNN],ylabel=L"\theta_7")
plot!(p[1+7*NNN],ylabel=L"\theta_8")


plot!(p[7*NNN+1],xlabel=L"\delta_1")
plot!(p[7*NNN+2],xlabel=L"\delta_2")
plot!(p[7*NNN+3],xlabel=L"\delta_3")
plot!(p[7*NNN+4],xlabel=L"\delta_4")
plot!(p[7*NNN+5],xlabel=L"\delta_5")
plot!(p[7*NNN+6],xlabel=L"\delta_6")
plot!(p[7*NNN+7],xlabel=L"\delta_7")
plot!(p[7*NNN+8],xlabel=L"\delta_8")





for i in 1:NNN
for j in 1:NNN
plot!(p[NNN*(i-1)+j],xlims=(-4.0,4.0))
end
end

for i in 1:NNN
for j in 1:NNN
plot!(p[(i-1)*NNN+j],ylims=(-4.0,4.0))
end
end

for i in 1:NNN
plot!(p[(i-1)*NNN+i],ylims=(0.0,1.2))
annotate!(p[(i-1)*NNN+i],1,1.0,text("Counts",10))
end
for i in 1:NNN*NNN
plot!(p[0+i],right_margin=-3mm)
end



labels = [L"~~\mathrm{Posterior}~68~\%", L"~~\mathrm{Posterior}~95~\%",L"~~\mathrm{Posterior}~99~\%"]
prior_labels = [L"~~\mathrm{Prior}~68~\%", L"~~\mathrm{Prior}~95~\%",L"~~\mathrm{Posterior}~99~\%"]

plot!(legend=false,label="xx",
    marginalmode=false,  interval_labels=prior_labels, 
    colors=reverse([c1, c2, c3]), linewidth=0,
    grid=false,foreground_color_subplot=:white,subplot=NNN-2
    , legendfontsize=14   
    , xlims=(0.0,1.0)
    , ylims=(0.0,1.0)
)
plot!(legend=false,label="xx",
    marginalmode=false,  interval_labels=prior_labels, 
    colors=reverse([c1, c2, c3]), linewidth=0,
    grid=false,foreground_color_subplot=:white,subplot=2*NNN-2
    , legendfontsize=14   
    , xlims=(0.0,1.0)
    , ylims=(0.0,1.0)
)
plot!(legend=false,label="xx",
    marginalmode=false,  interval_labels=prior_labels, 
    colors=reverse([c1, c2, c3]), linewidth=0,
    grid=false,foreground_color_subplot=:white,subplot=3*NNN-2
    , legendfontsize=14   
    , xlims=(0.0,1.0)
    , ylims=(0.0,1.0)
)
plot!(legend=false,label="xx",
    marginalmode=false,  interval_labels=prior_labels, 
    colors=reverse([c1, c2, c3]), linewidth=0,
    grid=false,foreground_color_subplot=:white,subplot=NNN-1
    , legendfontsize=14   
    , xlims=(0.0,1.0)
    , ylims=(0.0,1.0)
)
plot!(legend=false,label="xx",
    marginalmode=false,  interval_labels=prior_labels, 
    colors=reverse([c1, c2, c3]), linewidth=0,
    grid=false,foreground_color_subplot=:white,subplot=2*NNN-1
    , legendfontsize=14   
    , xlims=(0.0,1.0)
    , ylims=(0.0,1.0)
)
plot!(legend=false,label="xx",
    marginalmode=false,  interval_labels=prior_labels, 
    colors=reverse([c1, c2, c3]), linewidth=0,
    grid=false,foreground_color_subplot=:white,subplot=3*NNN-1
    , legendfontsize=14   
    , xlims=(0.0,1.0)
    , ylims=(0.0,1.0)
)
#rectangle(w, h, x, y) = Shape(x .+ [0.0,w,w,0.0], y .+ [0.0,0.0,h,h])
#plot!(rectangle(0.05,0.05,0.5,0.5),subplot=4)
plot!(Shape([0.00,0.00,0.2,0.2],[0.0,0.2,0.2,0.0]),subplot=NNN-2,fillcolor=c3,alpha=alpha_posterior)
plot!(Shape([0.00,0.00,0.2,0.2],[0.4,0.6,0.6,0.4]),subplot=2*NNN-2,fillcolor=c2,alpha=alpha_posterior)
plot!(Shape([0.00,0.00,0.2,0.2],[1.0,0.8,0.8,1.0]),subplot=3*NNN-2,fillcolor=c1,alpha=alpha_posterior)
annotate!(p[NNN-1],0.0,0.10,text(L"~~~~\mathrm{Posterior}~99~\%",18))
annotate!(p[2*NNN-1],0.0,0.50,text(L"~~~~\mathrm{Posterior}~95~\%",18))
annotate!(p[3*NNN-1],0.0,0.90,text(L"~~~~\mathrm{Posterior}~68~\%",18))


plot(p)

savefig(string("figures/fig10-corner-",parsed_args["fitresults"],"_v3.pdf"))

end

main()
