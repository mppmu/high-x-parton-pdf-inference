#!/usr/bin/julia
using BAT, DensityInterface
#Pkg.add("QCDNUM")
#Pkg.add(url="https://github.com/cescalara/PartonDensity.jl.git")
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

Ns = 10000 # Number of samples from posterior
rn = MersenneTwister(seed);
sub_samples = BAT.bat_sample(rn, samples_data, BAT.OrderedResampling(nsamples=Ns)).result;

forward_model_init(qcdnum_params, splint_params)

counts_em_sampled = zeros(UInt64, (length(sub_samples), nbins))
counts_ep_sampled = zeros(UInt64, (length(sub_samples), nbins))
chisqep = zeros( length(sub_samples))
chisqem = zeros( length(sub_samples))


rng = MersenneTwister(seed);
sys_err_params = rand(rng, MvNormal(zeros(PartonDensity.nsyst), zeros(PartonDensity.nsyst)))

for s in eachindex(sub_samples)

    pdf_params_s = DirichletPDFParams(K_u=sub_samples.v.K_u[s], K_d=sub_samples.v.K_d[s], K_q=sub_samples.v.K_q[s],
                                      λ_g1=sub_samples.v.λ_g1[s], 
                                      λ_g2=sub_samples.v.λ_g2[s],
                                      K_g=sub_samples.v.K_g[s], 
                                      λ_q=sub_samples.v.λ_q[s], 
                                      θ=Vector(sub_samples.v.θ[s]))
        
    counts_ep_pred_s, counts_em_pred_s = forward_model(pdf_params_s,    
                                                              qcdnum_params, 
                                                              splint_params,
                                                              quark_coeffs, 
                                                              sys_err_params)
    
    for j in 1:nbins
        
        counts_ep_pred_s[j] *= 1 + 0.018 * sub_samples.v.Beta1[s]
        counts_em_pred_s[j] *= 1 + 0.018 * sub_samples.v.Beta2[s]
        
        counts_em_sampled[s, j] = rand(Poisson(counts_em_pred_s[j]))
        counts_ep_sampled[s, j] = rand(Poisson(counts_ep_pred_s[j]))

        chisqep[s]+=(counts_ep_pred_s[j]-counts_ep_sampled[s, j])^2/counts_ep_pred_s[j]
        chisqem[s]+=(counts_em_pred_s[j]-counts_em_sampled[s, j])^2/counts_em_pred_s[j]

    end
    
end

#
# Chi squared for the actual data
#
pdfpars(params)=   DirichletPDFParams(
    K_u=params.K_u, K_d=params.K_d, λ_g1=params.λ_g1, λ_g2=params.λ_g2, K_q=params.K_q,
    K_g=params.K_g, λ_q=params.λ_q, θ=Vector(params.θ))
mode_pars_data = mode(samples_data)
println(mode_pars_data)

pdf_params = pdfpars(mode_pars_data)
println(pdf_params)
sys_err_params =[
    mode_pars_data.beta0_1,mode_pars_data.beta0_2,mode_pars_data.beta0_3,mode_pars_data.beta0_4,
    mode_pars_data.beta0_5,mode_pars_data.beta0_6,mode_pars_data.beta0_7,mode_pars_data.beta0_8]
    counts_pred_ep_data, counts_pred_em_data = forward_model(pdf_params, qcdnum_params, splint_params, quark_coeffs, sys_err_params)

for i in 1:nbins     
    counts_pred_ep_data[i] =counts_pred_ep_data[i]*(1+0.018*mode_pars_data.Beta1)
    counts_pred_em_data[i] =counts_pred_em_data[i]*(1+0.018*mode_pars_data.Beta2)
end


#
# Calculate the Poisson probabilities for the different data results
#
chisqep_data=0.
chisqem_data=0.
        for j in 1:nbins
            pred=counts_pred_ep_data[j]   
            best=floor(counts_pred_ep_data[j])
            prob_ep_data[j] = pdf(Poisson(pred), counts_obs_ep_data[j])/pdf(Poisson(pred), best)
            chisqep_data+=(counts_obs_ep_data[j]-pred)^2/pred
            if ( (counts_obs_ep_data[j]-pred)^2/pred>4) 
                get_bin_info(j) 
                println(j," positron ",pred," ",counts_obs_ep_data[j]," ",(counts_obs_ep_data[j]-pred)^2/pred)
            end
            pred=counts_pred_em_data[j]
            best=floor(counts_pred_em_data[j])
            prob_em_data[j] = pdf(Poisson(pred), counts_obs_em_data[j])/pdf(Poisson(pred), best)
            chisqem_data+=(counts_obs_em_data[j]-pred)^2/pred
            if ( (counts_obs_em_data[j]-pred)^2/pred>4) 
                get_bin_info(j) 
                println(j," electron ",pred," ",counts_obs_em_data[j]," ",(counts_obs_em_data[j]-pred)^2/pred)
            end 
       end
println(chisqep_data," ",chisqem_data)

countep=0
countem=0
count_tot=0
ncounts=0

for s in eachindex(sub_samples)
    ncounts+=1
    if (chisqep[s]>chisqep_data) countep+=1 end
    if (chisqem[s]>chisqem_data) countem+=1 end
    if (chisqem[s]+chisqep[s]>chisqem_data+chisqep_data) count_tot+=1 end

end
pvep=string(float(countep)/ncounts)
pvem=string(float(countem)/ncounts)
println(" p-value for ep fit: ",pvep," p-value for em fit: ",pvem," p-value for total fit: ",float(count_tot)/ncounts) 


p1=histogram(chisqep,bins=100,xlabel=L"\chi^2_P", ylabel="Entries", fontfamily=font_family,color=c1 , grid=false)
p1=plot!([chisqep_data],seriestype = :vline,lw=5,legend=:none, fontfamily=font_family 
, xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
        ,ylims=(0, 500), xlims=(50,260)
        ,color=c3, grid=false
)
p2=histogram(chisqem,bins=50:2:250,legend=:false,xlabel=L"\chi^2_P", ylabel="Entries", fontfamily=font_family,color=c1 , grid=false)

p2=plot!([chisqem_data],seriestype = :vline,lw=5
, xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    ,ylims=(0, 500), xlims=(50,260)
,color=c3, grid=false
)
annotate!(p1,220.0,350,text(L"$e^{+}p$",26))
annotate!(p2,220.0,350,text(L"$e^{-}p$",26))
plot(p1,p2,layout=(2,1))

filename = string("figures/fig8-chisq-pvalue-",parsed_args["fitresults"],"_v2.pdf")
savefig(filename)

# Use +2 to avoid lightest colors (not easy to see)

c1 = :teal
c2 = :royalblue4
c3 = :midnightblue
c4 = :grey
cmap = palette(color_scheme, n_q2_bins+2)
Plots.scalefontsizes()
Plots.scalefontsizes(1.2);
alpha = 0.6
prior_alpha = 0.2;

# Get some prior samples for plotting

if parsed_args["parametrisation"] == "Dirichlet"
if (parsed_args["priorshift"]==0)
    println("seting prior from Shifted Prior set ",seedtxt)

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
end

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

prior = NamedTupleDist(
    θ = Dirichlet([34.0, 17.0, 22.5, 17.0, 7.3, 1.4, 0.2, 10.e-5, 0.3]),
    initial_U = Uniform(-10., 1.),
    initial_D = Uniform(10., 30.),
    λ_g1 = Uniform(3., 4.5),
    λ_g2 = Uniform(-1, -0.5),
    K_g =  Uniform(5.,9.),
    λ_q = Uniform(-1, -0.5),
    K_q = Uniform(3., 7.),
    bspoly_params = [[0, 3], [0, 4], [1, 4], [0, 5]],
    #    bspoly_params = [1,4,0,4,0,5],
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

prior_samples=bat_sample(prior).result;

xlims_K_u = (2.0, 7.0) # (3.2, 4.4)
xlims_D_u = (0., 0.5) # (0.29, 0.37)
xlims_K_d = (2.0, 7.0) # (3.2, 4.4)
xlims_D_d = (0., 0.5) # (0.29, 0.37)
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



l = @layout [[grid(4,4)]]
#pgfplotsx()
l0 = @layout [
    [a1{0.28w}  _ _ _ ]
        [a2{0.28w} b2{0.24w} _ _ ]
        [a3{0.28w} b3{0.24w} c3{0.24w} _ ]
        [a4{0.28w} b4{0.24w} c4{0.24w} d4{0.24w} ]
]
l = @layout [
    #[a1{0.28w}  _ _ _ ]
      #  [a2{0.28w} b2{0.24w} _ _ ]
      #  [a3{0.28w} b3{0.24w} c3{0.24w} _ ]
      #  [a4{0.28w} b4{0.24w} c4{0.24w} d4{0.24w} ]

        [a1{0.28w}  b1{0.24w} c1{0.24w} d1{0.24w} ]
        [a2{0.28w} b2{0.24w} c2{0.24w} d2{0.24w} ]
        [a3{0.28w} b3{0.24w} c3{0.24w} d3{0.24w} ]
        [a4{0.28w} b4{0.24w} c4{0.24w} d4{0.24w} ]

]

NNN=4
p=plot(size=(1000,800),
layout=l,
#samples_data,vsel=[:K_u,:K_d, :K_q, :K_g], 
    #cmap=:viridis,
    colors=[c1, c2, :red],
  #  colors=[:chartreuse2, :yellow, :red]
    frame=:box,
    legend = :none, framestyle = :box,
    ylims=(1.9, 10.1),xlims=(1.9, 10.1), 
    plot_titlevspan=0.001,
     xlabel=["" "" "" "" "" "" "" "" "" "" "" "" L"p(K_{u})" L"K_{d}" L"K_{q}" L"K_{g}"],
    ylabel=[L"p(K_{u})" "" "" "" L"K_{d}" "" "" "" L"K_{q}" "" "" "" L"K_{g}" "" "" ""],
    xticks=(0:2:12,["","","","","","",""]),    
    yticks=(0:2:12,["","","","","","",""])  ,
    top_margin=-3mm,
    bottom_margin=-3.5mm,
    right_margin=-4mm,
    left_margin=-4mm,
    xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16
    , fontfamily=font_family
    ,grid=false

)
plot!(samples_data, :(K_u),subplot=1, xlabel="",ylabel=L"K_u",colors=[c1, c2, :red])
plot!(legend=false,grid=false,foreground_color_subplot=:white,subplot=2)
plot!(legend=false,grid=false,foreground_color_subplot=:white,subplot=3)
plot!(legend=false,grid=false,foreground_color_subplot=:white,subplot=4)
plot!(samples_data, (:(K_u), :(K_d)),subplot=5, xlabel="",ylabel=L"K_d",colors=[c1, c2, :red])
plot!(samples_data, :(K_d),subplot=6, xlabel="",ylabel="",colors=[c1, c2, :red])
plot!(legend=false,grid=false,foreground_color_subplot=:white,subplot=7)
plot!(legend=false,grid=false,foreground_color_subplot=:white,subplot=8)
plot!(samples_data, (:(K_u), :(K_q)),subplot=9, xlabel="",ylabel=L"K_q",colors=[c1, c2, :red])
plot!(samples_data, (:(K_d), :(K_q)),subplot=10, xlabel="",ylabel="",colors=[c1, c2, :red])
plot!(samples_data, :(K_q),subplot=11, xlabel="",ylabel="",colors=[c1, c2, :red])
plot!(legend=false,grid=false,foreground_color_subplot=:white,subplot=12)
plot!(samples_data, (:(K_u), :(K_g)),subplot=13, xlabel=L"K_u",ylabel=L"K_g",colors=[c1, c2, :red])
plot!(samples_data, (:(K_d), :(K_g)),subplot=14, xlabel=L"K_d",ylabel="",colors=[c1, c2, :red])
plot!(samples_data, (:(K_q), :(K_g)),subplot=15, xlabel=L"K_q",ylabel="",colors=[c1, c2, :red])
plot!(samples_data, :(K_g),subplot=16, xlabel=L"K_g",ylabel="",colors=[c1, c2, :red])

for i in 0:NNN-1
plot!(p[NNN*i+1],xlims=(2.75,4.5))#1.5,10.5
plot!(p[NNN*i+2],xlims=(2.0,5.5))#1.5,10.5
plot!(p[NNN*i+3],xlims=(3.0,9.5))
plot!(p[NNN*i+4],xlims=(2.0,7.0))
end
for i in 1:NNN
plot!(p[0+i],ylims=(1.5,10.5))
plot!(p[1*NNN+i],ylims=(1.5,6.5))
plot!(p[2*NNN+i],ylims=(3.0,9.0))
plot!(p[3*NNN+i],ylims=(2.0,7.0))
end
for i in 1:NNN*NNN

plot!(p[0+i],right_margin=-3mm)
end

plot!(p[13],xticks=(3.0:0.5:4.0,["3","3.5","4"]),bottom_margin=3mm)
plot!(p[14],xticks=(2:1:5,["","3", "4","5"]),bottom_margin=3mm)
plot!(p[15],xticks=(4:2:9,["4","6","8"]),bottom_margin=3mm)
plot!(p[16],xticks=(3:1:6,["3","4","5","6"]),bottom_margin=3mm)
plot!(p[1],yticks=(2:2:10,["","4","6","8","10"]),left_margin=3mm)
plot!(p[5],yticks=(2:2:6,["2","4","6"]),left_margin=3mm)
plot!(p[9],yticks=(4:2:9,["4","6","8"]),left_margin=3mm)
plot!(p[13],yticks=(2:2:7,["2","4","6"]),left_margin=3mm)
plot!(p[1],ylims=(0.0,3.0))
annotate!(p[1],2.75+1.75*0.7,0.8*3.0,text("Counts",14))
plot!(p[6],ylims=(0.0,1))
annotate!(p[6],2.0+3.5*0.7,0.8*1.0,text("Counts",14))
plot!(p[11],ylims=(0.0,0.5))
annotate!(p[11],3.0+6.5*0.7,0.8*0.5,text("Counts",14))
plot!(p[16],ylims=(0.0,0.5))
annotate!(p[16],2.0+5.0*0.7,0.8*0.5,text("Counts",14))


plot(p)

savefig(string("figures/fig2-corner-",parsed_args["fitresults"],"_v3.pdf"))

end

main()
