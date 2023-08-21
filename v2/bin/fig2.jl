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
include("priors.jl")
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

c1 = :teal
c2 = :royalblue4
c3 = :midnightblue
c4 = :grey

c5 = :grey

#c2=colorant"#CCE5E5"
#c1=colorant"#93A0AB"

c2 = :teal
c1 = :midnightblue
c3 = :grey

color_scheme = :viridis
font_family = "Computer Modern"
default(fontfamily = "Computer Modern")

# Results
seed=parsed_args["seed"]
println(seed)
seedtxt=string(seed)

#Sim data!!!
pdf_params, sim_data=pd_read_sim(string("pseudodata/", parsed_args["pseudodata"], ".h5"))

#Fit results!!!
samples_data = bat_read(string("fitresults/", parsed_args["fitresults"], ".h5")).result;



samples_data1 = bat_transform(v -> (Δ_g = v.θ[3] + v.θ[4], Δ_sea = sum(v.θ[5:9]), K_u = v.K_u, K_q=v.K_q, K_g = v.K_g, K_d=v.K_d ), samples_data).result



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
    (q2_edges, x_edges) = get_bin_info(i, quiet=true);
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

Ns = 100000 # Number of samples from posterior
rn = MersenneTwister(seed);
sub_samples = BAT.bat_sample(rn, samples_data, BAT.OrderedResampling(nsamples=Ns)).result;

forward_model_init(qcdnum_params, splint_params)

counts_em_sampled = zeros(UInt64, (length(sub_samples), nbins))
counts_ep_sampled = zeros(UInt64, (length(sub_samples), nbins))
chisqep = zeros( length(sub_samples))
chisqem = zeros( length(sub_samples))


rng = MersenneTwister(seed);
sys_err_params = rand(rng, MvNormal(zeros(PartonDensity.nsyst), zeros(PartonDensity.nsyst)))



# Use +2 to avoid lightest colors (not easy to see)



cmap = palette(color_scheme, n_q2_bins+2)
alpha = 0.6
prior_alpha = 0.2;

# Get some prior samples for plotting
prior=get_priors(parsed_args)

prior_samples =bat_sample(prior).result;
prior_samples1 = bat_transform(v -> (Δ_g = v.θ[3] + v.θ[4], Δ_sea = sum(v.θ[5:9]), K_u = v.K_u, K_q=v.K_q, K_g = v.K_g, K_d=v.K_d ), prior_samples).result



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

NNN = 6

l = @layout [[grid(5,5)]]
#pgfplotsx()



l = @layout [
        [a1{0.20w}  b1{0.18w} c1{0.18w} d1{0.18w} e1{0.18w}]
        [a2{0.20w} b2{0.18w} c2{0.18w} d2{0.18w} e2{0.18w}]
        [a3{0.20w} b3{0.18w} c3{0.18w} d3{0.18w} e3{0.18w}]
        [a4{0.20w} b4{0.18w} c4{0.18w} d4{0.18w} e4{0.18w} ]
        [a5{0.20w} b5{0.18w} c5{0.18w} d5{0.18w} e5{0.18w} ]

]

l = @layout [
        [a1{0.18w}  b1{0.16w} c1{0.16w} d1{0.16w} e1{0.16w} f1{0.16w}]
        [a2{0.18w} b2{0.16w} c2{0.16w} d2{0.16w} e2{0.16w} f2{0.16w}]
        [a3{0.18w} b3{0.16w} c3{0.16w} d3{0.16w} e3{0.16w} f3{0.16w} ]
        [a4{0.18w} b4{0.16w} c4{0.16w} d4{0.16w} e4{0.16w} f4{0.16w} ]
        [a5{0.18w} b5{0.16w} c5{0.16w} d5{0.16w} e5{0.16w} f5{0.16w}]
        [a6{0.18w} b6{0.16w} c6{0.16w} d6{0.16w} e6{0.16w} f6{0.16w}]

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
  #  xticks=(0:2:12,["","","","","","",""]),    
  #  yticks=(0:2:12,["","","","","","",""])  ,
    xticks=:none,
    yticks=:none,
    top_margin=-3mm,
    bottom_margin=-1mm,
    right_margin=-4mm,
    left_margin=-4mm,
    xtickfontsize=14,
    ytickfontsize=14,
    yguidefontsize=16,
    xguidefontsize=16,
    fontfamily=font_family,
    grid=false
)
plot!(samples_data, :(K_u),subplot=1, xlabel="",ylabel=L"K_u",colors=[c1, c2, c3])
plot!(legend=false,grid=false,foreground_color_subplot=:white,subplot=2)
plot!(legend=false,grid=false,foreground_color_subplot=:white,subplot=3)
plot!(legend=false,grid=false,foreground_color_subplot=:white,subplot=4)
plot!(legend=false,grid=false,foreground_color_subplot=:white,subplot=5)
#plot!(legend=false,grid=false,foreground_color_subplot=:white,subplot=4)



plot!(samples_data, (:(K_u), :(K_d)),subplot=NNN+1, xlabel="",ylabel=L"K_d",colors=[c1, c2, c3],markerstrokewidth = 0)
plot!(samples_data, :(K_d),subplot=NNN+2, xlabel="",ylabel="",colors=[c1, c2, c3])
plot!(legend=false,grid=false,foreground_color_subplot=:white,subplot=NNN+3)
plot!(legend=false,grid=false,foreground_color_subplot=:white,subplot=NNN+4)
plot!(legend=false,grid=false,foreground_color_subplot=:white,subplot=NNN+5)
plot!(legend=false,grid=false,foreground_color_subplot=:white,subplot=NNN+6)

plot!(samples_data, (:(K_u), :(K_q)),subplot=2*NNN+1, xlabel="",ylabel=L"K_q",colors=[c1, c2, c3])
plot!(samples_data, (:(K_d), :(K_q)),subplot=2*NNN+2, xlabel="",ylabel="",colors=[c1, c2, c3])
plot!(samples_data, :(K_q),subplot=2*NNN+3, xlabel="",ylabel="",colors=[c1, c2, c3])
plot!(legend=false,grid=false,foreground_color_subplot=:white,subplot=2*NNN+4)
plot!(legend=false,grid=false,foreground_color_subplot=:white,subplot=2*NNN+5)
plot!(legend=false,grid=false,foreground_color_subplot=:white,subplot=2*NNN+6)

plot!(samples_data, (:(K_u), :(K_g)),subplot=3*NNN+1, xlabel="",ylabel=L"K_g",colors=[c1, c2, c3])
plot!(samples_data, (:(K_d), :(K_g)),subplot=3*NNN+2, xlabel="",ylabel="",colors=[c1, c2, c3])
plot!(samples_data, (:(K_q), :(K_g)),subplot=3*NNN+3, xlabel="",ylabel="",colors=[c1, c2, c3])
plot!(samples_data, :(K_g),subplot=3*NNN+4, xlabel="",ylabel="",colors=[c1, c2, c3])
plot!(legend=false,grid=false,foreground_color_subplot=:white,subplot=3*NNN+5)
plot!(legend=false,grid=false,foreground_color_subplot=:white,subplot=3*NNN+6)


plot!(samples_data1, (:(K_u), :(Δ_g)),subplot=4*NNN+1, xlabel="",ylabel=L"\Delta_g",colors=[c1, c2, c3])
plot!(samples_data1, (:(K_d), :(Δ_g)),subplot=4*NNN+2, xlabel="",ylabel="",colors=[c1, c2, c3])
plot!(samples_data1, (:(K_q), :(Δ_g)),subplot=4*NNN+3, xlabel="",ylabel="",colors=[c1, c2, c3])
plot!(samples_data1, (:(K_g), :(Δ_g)),subplot=4*NNN+4, xlabel="",ylabel="",colors=[c1, c2, c3])
plot!(samples_data1, :(Δ_g),subplot=4*NNN+5, xlabel="",ylabel="",colors=[c1, c2, c3])
plot!(legend=false,grid=false,foreground_color_subplot=:white,subplot=4*NNN+6)


plot!(samples_data1, (:(K_u), :(Δ_sea)),subplot=5*NNN+1, xlabel=L"K_u",ylabel=L"\Delta_{\mathrm{sea}}",colors=[c1, c2, c3])
plot!(samples_data1, (:(K_d), :(Δ_sea)),subplot=5*NNN+2, xlabel=L"K_d",ylabel="",colors=[c1, c2, c3])
plot!(samples_data1, (:(K_q), :(Δ_sea)),subplot=5*NNN+3, xlabel=L"K_q",ylabel="",colors=[c1, c2, c3])
plot!(samples_data1, (:(K_g), :(Δ_sea)),subplot=5*NNN+4, xlabel=L"K_g",ylabel="",colors=[c1, c2, c3])
plot!(samples_data1, (:(Δ_g),:(Δ_sea)), subplot=5*NNN+5, xlabel=L"\Delta_g",ylabel="",colors=[c1, c2, c3])
plot!(samples_data1, :(Δ_sea),subplot=5*NNN+6, xlabel=L"\Delta_{\mathrm{sea}}",ylabel="",colors=[c1, c2, c3])





for i in 0:NNN-1
plot!(p[NNN*i+1],xlims=(2.75,4.5))#1.5,10.5
plot!(p[NNN*i+2],xlims=(2.0,5.5))#1.5,10.5
plot!(p[NNN*i+3],xlims=(2.0,9.5))
plot!(p[NNN*i+4],xlims=(1.0,7.0))
plot!(p[NNN*i+5],xlims=(0.35,0.65))
plot!(p[NNN*i+6],xlims=(0.05,0.35))
end
for i in 1:NNN
plot!(p[0+i],ylims=(0.0,10.5))
plot!(p[1*NNN+i],ylims=(1.5,6.5))
plot!(p[2*NNN+i],ylims=(1.0,9.0))
plot!(p[3*NNN+i],ylims=(1.0,9.0))
plot!(p[4*NNN+i],ylims=(0.33,0.65))
plot!(p[5*NNN+i],ylims=(0.05,0.35))
end
for i in 1:NNN*NNN

plot!(p[0+i],right_margin=-3mm)
end

plot!(p[NNN*(NNN-1)+ 1 ],xticks=(3.0:0.5:4.0,["3","3.5","4"]),bottom_margin=5.0mm)
plot!(p[NNN*(NNN-1)+ 2 ],xticks=(2:1:5,["","3", "4","5"]),bottom_margin=5.0mm)
plot!(p[NNN*(NNN-1)+ 3 ],xticks=(4:2:9,["4","6","8"]),bottom_margin=5.0mm)
plot!(p[NNN*(NNN-1)+ 4 ],xticks=(2:2:6,["2","4","6"]),bottom_margin=5.0mm)
plot!(p[NNN*(NNN-1)+ 5 ],xticks=(0.4:0.1:0.6,["0.4","0.5","0.6"]),bottom_margin=5.0mm)
plot!(p[NNN*(NNN-1)+ 6 ],xticks=(0.1:0.1:0.3,["0.1","0.2","0.3"]),bottom_margin=5.0mm)
#plot!(p[      1],yticks=(2:2:10,["","4","6","8","10"]),left_margin=7.0mm)
plot!(p[  NNN+1],yticks=(2:2:6,["2","4","6"]),left_margin=7.0mm)
plot!(p[2*NNN+1],yticks=(2:2:9,["2","4","6","8"]),left_margin=7.0mm)
plot!(p[3*NNN+1],yticks=(2:2:9,["2","4","6","8"]),left_margin=7.0mm)
plot!(p[4*NNN+1],yticks=(0.4:0.1:0.6,["0.4","0.5","0.6"]),left_margin=7.0mm)
plot!(p[5*NNN+1],yticks=(0.1:0.1:0.3,["0.1","0.2","0.3"]),left_margin=7.0mm)

ky=0.8
kx=0.5
plot!(p[1],ylims=(0.0,4.0),right_margin=0.1mm)
annotate!(p[1],2.75+1.75*kx,ky*4.0,text("Counts",14))
plot!(p[NNN+2],ylims=(0.0,1.3),right_margin=0.1mm)
annotate!(p[NNN+2],2.0+3.5*kx,ky*1.3,text("Counts",14))
plot!(p[2*NNN+3],ylims=(0.0,0.5),right_margin=0.1mm)
annotate!(p[2*NNN+3],3.0+6.5*kx,ky*0.5,text("Counts",14))
plot!(p[3*NNN+4],ylims=(0.0,0.5),right_margin=0.1mm)
annotate!(p[3*NNN+4],2.0+5.0*kx,ky*0.5,text("Counts",14))
plot!(p[4*NNN+5],ylims=(0.0,20),right_margin=0.1mm)
annotate!(p[4*NNN+5],0.35+0.3*kx,ky*20,text("Counts",14))

plot!(p[5*NNN+6],ylims=(0.0,30),right_margin=0.1mm)
annotate!(p[5*NNN+6],0.05+0.3*kx,ky*30,text("Counts",14))



labels = [L"~~\mathrm{Posterior}~68~\%", L"~~\mathrm{Posterior}~95~\%",L"~~\mathrm{Posterior}~99~\%"]
prior_labels = [L"~~\mathrm{Prior}~68~\%", L"~~\mathrm{Prior}~95~\%",L"~~\mathrm{Posterior}~99~\%"]

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
plot!(legend=false,label="xx",
    marginalmode=false,  interval_labels=prior_labels, 
    colors=reverse([c1, c2, c3]), linewidth=0,
    grid=false,foreground_color_subplot=:white,subplot=NNN
    , legendfontsize=14   
    , xlims=(0.0,1.0)
    , ylims=(0.0,1.0)
)
plot!(legend=false,label="xx",
    marginalmode=false,  interval_labels=prior_labels, 
    colors=reverse([c1, c2, c3]), linewidth=0,
    grid=false,foreground_color_subplot=:white,subplot=2*NNN
    , legendfontsize=14   
    , xlims=(0.0,1.0)
    , ylims=(0.0,1.0)
)
plot!(legend=false,label="xx",
    marginalmode=false,  interval_labels=prior_labels, 
    colors=reverse([c1, c2, c3]), linewidth=0,
    grid=false,foreground_color_subplot=:white,subplot=3*NNN
    , legendfontsize=14   
    , xlims=(0.0,1.0)
    , ylims=(0.0,1.0)
)
#rectangle(w, h, x, y) = Shape(x .+ [0.0,w,w,0.0], y .+ [0.0,0.0,h,h])
#plot!(rectangle(0.05,0.05,0.5,0.5),subplot=4)
plot!(Shape([0.00,0.00,0.2,0.2],[0.0,0.2,0.2,0.0]),subplot=NNN-1,fillcolor=c3)
plot!(Shape([0.00,0.00,0.2,0.2],[0.4,0.6,0.6,0.4]),subplot=2*NNN-1,fillcolor=c2)
plot!(Shape([0.00,0.00,0.2,0.2],[1.0,0.8,0.8,1.0]),subplot=3*NNN-1,fillcolor=c1)
annotate!(p[NNN],0.0,0.10,text(L"~~\mathrm{Posterior}~99~\%",18))
annotate!(p[2*NNN],0.0,0.50,text(L"~~\mathrm{Posterior}~95~\%",18))
annotate!(p[3*NNN],0.0,0.90,text(L"~~\mathrm{Posterior}~68~\%",18))


plot(p)

savefig(string("figures/fig2-corner-",parsed_args["fitresults"],"_v3.pdf"))

end

main()
