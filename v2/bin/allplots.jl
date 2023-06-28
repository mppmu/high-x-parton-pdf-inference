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


plot(framestyle=:axes, size=(500, 400), fontfamily=font_family, 
    layout=@layout([a b; c{0.55w, 0.6h} d]), grid=false
    , right_margin=0mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=0mm
)

# Joint posterior
plot!(prior_samples, (:(K_u), :(θ[1])), xlabel=L"K_u", ylabel=L"\Delta_u",
    seriestype=:smallest_intervals_contourf, smoothing=4, 
    marginalmode=false, intervals=intervals, fillcolors=reverse(prior_colors), linewidth=0, 
    alpha=prior_alpha,
    subplot=3
     , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , right_margin=-2mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm

)
plot!(samples_data, (:(K_u), :(θ[1])), xlabel=L"K_u", ylabel=L"\Delta_u",
    seriestype=:smallest_intervals_contourf, smoothing=2, 
    marginalmode=false, intervals=intervals, fillcolors=reverse(colors), linewidth=0, alpha=alpha,
    subplot=3, xlims=xlims_K_u, ylims=xlims_D_u
 , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , right_margin=-2mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm

)
p = plot!([K_u_true],[θ_true[1]], color="red",subplot=3, seriestype=:scatter, label=" Truth", lw=0, foreground_color_legend=false, markersize=3, thickness_scaling=1.0, lc=:red, markerstrokecolor=:red, legendfontsize=18)



# K_u marginal
plot!(prior_samples, :K_u, legend=false, marginalmode=false, 
    seriestype=:smallest_intervals, intervals=intervals, 
    colors=prior_colors, subplot=1, alpha=prior_alpha
 , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , right_margin=-2mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm

)
plot!(samples_data, :K_u, legend=false, xlabel="", ylabel=L"P(K_u)", subplot=1, 
    xlims=xlims_K_u, ylims=(0, 4), seriestype=:smallest_intervals, 
    marginalmode=false, intervals=intervals, colors=colors, alpha=alpha
 , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , right_margin=-2mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm

)
vline!([K_u_true], color="red", label=" Truth", lw=0.5)

# Delta_u marginal
plot!(prior_samples, :(θ[1]), legend=false, marginalmode=false, 
    seriestype=:smallest_intervals, intervals=intervals,
    colors=prior_colors, subplot=4, alpha=prior_alpha, 
    orientation=:horizontal
 , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , right_margin=-2mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm

)
plot!(samples_data, :(θ[1]), legend=false, ylabel="", xlabel=L"P(\Delta_u)", 
    subplot=4, ylims=xlims_D_u, xlims=(0, 55), 
    seriestype=:smallest_intervals, intervals=intervals, marginalmode=false, 
    colors=colors, alpha=alpha, orientation=:horizontal
 , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , xticks=(0:20:40,["0","20","40"])
    , right_margin=-2mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm

)
hline!([θ_true[1]], color="red", label=" Truth", subplot=4, lw=0.5)

# Legend
plot!(prior_samples, (:(K_u), :(θ[1])), xlabel=L"K_u", ylabel=L"\Delta_u",
    seriestype=:smallest_intervals,
    marginalmode=false, intervals=intervals, interval_labels=prior_labels, 
    colors=reverse(prior_colors), linewidth=0, 
    alpha=prior_alpha+0.2,
    subplot=2
 , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , right_margin=-2mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm

)
p = plot!(samples_data, (:(K_u), :(θ[1])),
    seriestype=:smallest_intervals,
    marginalmode=false, intervals=intervals, colors=reverse(colors), 
    interval_labels=labels,
    linewidth=0, alpha=alpha+0.2, legend=:bottomleft, foreground_color_legend=false,
    framestyle=:none, subplot=2, xlims=(0, 1.), ylims=(0, 0.1)
 , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=16
    , right_margin=7mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm


)

filename = string("figures/fig3-",parsed_args["fitresults"], "_v2.pdf")
savefig(p, filename)




plot(framestyle=:axes, size=(500, 400), fontfamily=font_family, 
    layout=@layout([a b; c{0.55w, 0.6h} d]), grid=false
    , right_margin=0mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=0mm
)

# Joint posterior
plot!(prior_samples, (:(K_d), :(θ[2])), xlabel=L"K_d", ylabel=L"\Delta_d",
    seriestype=:smallest_intervals_contourf, smoothing=4, 
    marginalmode=false, intervals=intervals, fillcolors=reverse(prior_colors), linewidth=0, 
    alpha=prior_alpha,
    subplot=3
 , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , right_margin=-2mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm
)
plot!(samples_data, (:(K_d), :(θ[2])), xlabel=L"K_d", ylabel=L"\Delta_d",
    seriestype=:smallest_intervals_contourf, smoothing=2, 
    marginalmode=false, intervals=intervals, fillcolors=reverse(colors), linewidth=0, alpha=alpha,
    subplot=3, xlims=xlims_K_d, ylims=xlims_D_d
 , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , right_margin=-2mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm
)
p = plot!([K_d_true],[θ_true[2]], color="red",subplot=3, seriestype=:scatter, label=" Truth", lw=0, foreground_color_legend=false, markersize=3, thickness_scaling=1.0, lc=:red, markerstrokecolor=:red, legendfontsize=18)


# K_d marginal
plot!(prior_samples, :K_d, legend=false, marginalmode=false, 
    seriestype=:smallest_intervals, intervals=intervals, 
    colors=prior_colors, subplot=1, alpha=prior_alpha
 , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , right_margin=-2mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm
)
plot!(samples_data, :K_d, legend=false, xlabel="", ylabel=L"P(K_d)", subplot=1, 
    xlims=xlims_K_u, ylims=(0, 2), seriestype=:smallest_intervals, 
    marginalmode=false, intervals=intervals, colors=colors, alpha=alpha
 , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , right_margin=-2mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm
)
vline!([K_d_true], color="red", label=" Truth", lw=0.5)

# Delta_u marginal
plot!(prior_samples, :(θ[2]), legend=false, marginalmode=false, 
    seriestype=:smallest_intervals, intervals=intervals,
    colors=prior_colors, subplot=4, alpha=prior_alpha, 
    orientation=:horizontal
 , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , right_margin=-2mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm
)
plot!(samples_data, :(θ[2]), legend=false, ylabel="", xlabel=L"P(\Delta_d)", 
    subplot=4, ylims=xlims_D_u, xlims=(0, 25), 
    seriestype=:smallest_intervals, intervals=intervals, marginalmode=false, 
    colors=colors, alpha=alpha, orientation=:horizontal
 , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , xticks=(5:10:25,["5","15","25"])
    , right_margin=-2mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm
)
hline!([θ_true[2]], color="red", label=" Truth", subplot=4, lw=0.5)

# Legend
plot!(prior_samples, (:(K_d), :(θ[2])), xlabel=L"K_d", ylabel=L"\Delta_d",
    seriestype=:smallest_intervals,
    marginalmode=false, intervals=intervals, interval_labels=prior_labels, 
    colors=reverse(prior_colors), linewidth=0, 
    alpha=prior_alpha+0.2,
    subplot=2
 , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , right_margin=-2mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm
)
p = plot!(samples_data, (:(K_d), :(θ[2])),
    seriestype=:smallest_intervals,
    marginalmode=false, intervals=intervals, colors=reverse(colors), 
    interval_labels=labels,
    linewidth=0, alpha=alpha+0.2, legend=:bottomleft, foreground_color_legend=false,
    framestyle=:none, subplot=2, xlims=(0, 1.), ylims=(0, 0.1)
 , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , right_margin=7mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm
)



filename = string("figures/fig3-d-",parsed_args["fitresults"],"_v2.pdf")
savefig(filename)



prior_color = :grey40
prior_alpha = 0.2
posterior_color = c2
posterior_alpha = 0.7
bins = 60

# Combine results for Δ_sea and Δ_g
# NB: for sea, our θ for quark/antiquark are already x2 - no need to do so here!
# in our x_q_x() in parametrisations_base.jl we take the input θ and halve it. 
comb_prior_samples = bat_transform(v -> (Δ_g = v.θ[3] + v.θ[4], Δ_sea = sum(v.θ[5:9])), prior_samples).result
comb_samples = bat_transform(v -> (Δ_g = v.θ[3] + v.θ[4], Δ_sea = sum(v.θ[5:9])), samples_data).result;

plot(framestyle=:axes, size=(700, 600), fontfamily=font_family, 
    layout=@layout([a b; c d]), grid=false)
plot!(samples_data, (:(θ[1]), :(θ[2])), subplot=1, xlabel=L"\Delta_{u}", ylabel=L"\Delta_{d}",
   seriestype=:smallest_intervals_contourf, smoothing=4, 
    marginalmode=false, intervals=intervals, fillcolors=reverse(colors), linewidth=0, 
    alpha=prior_alpha
, xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    ,ylims=(0,0.32),yticks=(0:0.1:0.3,["0","0.1","0.2","0.3"])
        ,xlims=(0.1,0.42),xticks=(0.1:0.1:0.4,["0.1","0.2","0.3","0.4"])
        , right_margin=-2mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm
)
plot!(prior_samples, (:(θ[1]), :(θ[2])), subplot=1, xlabel=L"\Delta_{u}", ylabel=L"\Delta_{d}",
   seriestype=:smallest_intervals_contourf, smoothing=4, 
    marginalmode=false, intervals=intervals, fillcolors=reverse(prior_colors), linewidth=0, 
    alpha=prior_alpha
, xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14   
    
)
plot!([θ_true[1]],[θ_true[2]], subplot=1, color="red",seriestype=:scatter, label=" Truth", lw=0, foreground_color_legend=false, markersize=3, thickness_scaling=1.0, lc=:red, markerstrokecolor=:red, legendfontsize=18)

comb_prior_samples = bat_transform(v -> (Δ_g = v.θ[3] + v.θ[4], Δ_u = v.θ[1]), prior_samples).result
comb_samples = bat_transform(v -> (Δ_g = v.θ[3] + v.θ[4], Δ_u = v.θ[1]), samples_data).result;

plot!(comb_samples, (:(Δ_u), :(Δ_g)),subplot=2,xlabel=L"\Delta_{u}", ylabel=L"\Delta_{g}",
   seriestype=:smallest_intervals_contourf, smoothing=4, 
    marginalmode=false, intervals=intervals, fillcolors=reverse(colors), linewidth=0, 
    alpha=prior_alpha
, xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    ,ylims=(0.25,0.72),yticks=(0.3:0.1:0.72,["0.3","0.4","0.5","0.6","0.7"])
        ,xlims=(0.1,0.42),xticks=(0.1:0.1:0.4,["0.1","0.2","0.3","0.4"])
        , right_margin=-2mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm
    
)
plot!(comb_prior_samples, (:(Δ_u), :(Δ_g)), subplot=2,xlabel=L"\Delta_{u}", ylabel=L"\Delta_{g}",
   seriestype=:smallest_intervals_contourf, smoothing=4, 
    marginalmode=false, intervals=intervals, fillcolors=reverse(prior_colors), linewidth=0, 
    alpha=prior_alpha
, xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
)
θ_g_true=θ_true[3]+θ_true[4]
plot!([θ_true[1]],[θ_g_true], subplot=2,color="red",seriestype=:scatter, label=:none, lw=1, markerstrokecolor=:red, markersize=3)


comb_prior_samples = bat_transform(v -> (Δ_sea = sum(v.θ[5:9]), Δ_u = v.θ[1]), prior_samples).result
comb_samples = bat_transform(v -> (Δ_sea = sum(v.θ[5:9]), Δ_u = v.θ[1]), samples_data).result;

plot!(comb_samples, (:(Δ_u), :(Δ_sea)),subplot=3,xlabel=L"\Delta_{u}", ylabel=L"\Delta_\mathrm{sea}",
   seriestype=:smallest_intervals_contourf, smoothing=4, 
    marginalmode=false, intervals=intervals, fillcolors=reverse(colors), linewidth=0, 
    alpha=prior_alpha
, xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    ,ylims=(0,0.32),yticks=(0:0.1:0.3,["0","0.1","0.2","0.3"])
        ,xlims=(0.1,0.42),xticks=(0.1:0.1:0.4,["0.1","0.2","0.3","0.4"]) 
        , right_margin=-2mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm
)
plot!(comb_prior_samples, (:(Δ_u), :(Δ_sea)), subplot=3,xlabel=L"\Delta_{u}", ylabel=L"\Delta_\mathrm{sea}",
   seriestype=:smallest_intervals_contourf, smoothing=4, 
    marginalmode=false, intervals=intervals, fillcolors=reverse(prior_colors), linewidth=0, 
    alpha=prior_alpha
, xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
)
θ_sea_true=sum(θ_true[5:9])
plot!([θ_true[1]],[θ_sea_true], subplot=3,color="red",seriestype=:scatter, label=:none, lw=2, markerstrokecolor=:red, markersize=3)


comb_prior_samples = bat_transform(v -> (Δ_sea = sum(v.θ[5:9]), Δ_g = v.θ[3] + v.θ[4]), prior_samples).result
comb_samples = bat_transform(v -> (Δ_sea = sum(v.θ[5:9]), Δ_g = v.θ[3] + v.θ[4]), samples_data).result;

plot!(comb_samples, (:(Δ_g), :(Δ_sea)),subplot=4, xlabel=L"\Delta_{g}", ylabel=L"\Delta_\mathrm{sea}",
   seriestype=:smallest_intervals_contourf, smoothing=4, 
    marginalmode=false, intervals=intervals, fillcolors=reverse(colors), linewidth=0, 
    alpha=prior_alpha
, xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
     ,ylims=(0,0.32),yticks=(0:0.1:0.3,["0","0.1","0.2","0.3"])
    ,xlims=(0.25,0.72),xticks=(0.3:0.1:0.72,["0.3","0.4","0.5","0.6","0.7"])
    , right_margin=-2mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm

)
plot!(comb_prior_samples, (:(Δ_g), :(Δ_sea)),subplot=4,xlabel=L"\Delta_{g}", ylabel=L"\Delta_\mathrm{sea}",
   seriestype=:smallest_intervals_contourf, smoothing=4, 
    marginalmode=false, intervals=intervals, fillcolors=reverse(prior_colors), linewidth=0, 
    alpha=prior_alpha
, xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
)
Δ_sea_true = sum(θ_true[5:9])
Δ_g_true = sum(θ_true[3:4])
plot!([Δ_g_true],[Δ_sea_true], color="red",subplot=4,seriestype=:scatter, label=:none, lw=2, markerstrokecolor=:red, markersize=3)


filename = string("figures/fig5-momentum-corr-", parsed_args["fitresults"], "_v2.pdf")
savefig(filename)

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
    scale = 0.1
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

plot(framestyle=:axes, size=(600, 500), fontfamily=font_family, 
    layout=@layout([a b; c d]), grid=false

    , right_margin=0mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=0mm

)


λ_u_true = θ_true[1]*(1+K_u_true)/(2-θ_true[1])
p = plot!(x_grid, [x_uv_x(x, λ_u_true, K_u_true) for x in x_grid], label=L"~xu_v \; \mathrm{true}", lw=3, c=:red, 
    subplot=1
     , ylims=(0, 0.65), xlims=(0.00,1.0)
    , foreground_color_legend=false
    , right_margin=0mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm
    , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
        , xticks=(0.0:0.2:1.0,["0","0.2","0.4","0.6","0.8","1"]) 
        , yticks=(0.0:0.2:0.6,["0","0.2","0.4","0.6"]) 
    
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

    plot!(x_grid, [x_uv_x(x, λ_u, pdf_params_s.K_u) for x in x_grid],lw=0.5,alpha=0.3,label=:none, subplot=1)

end

    p = plot!(xlabel=L"x", ylabel=L"xu_v", subplot=1)


p = plot!(x_grid, [x_uv_x(x, λ_u_true, K_u_true) for x in x_grid], label=L"~xu_v \; \mathrm{true}", lw=3, c=:red, 
    subplot=1
     , ylims=(0, 0.65), xlims=(0.00,1.0)
    , foreground_color_legend=false
    , right_margin=0mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm
    , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
        , xticks=(0.0:0.2:1.0,["0","0.2","0.4","0.6","0.8","1"]) 
        , yticks=(0.0:0.2:0.6,["0","0.2","0.4","0.6"]) 
    
)



#
λ_d_true = θ_true[2]*(1+K_d_true)/(1-θ_true[2])

p = plot!(x_grid, [x_dv_x(x, λ_d_true, K_d_true) for x in x_grid], label=L"~xd_v \; \mathrm{true}",  lw=3, c=:red,
    subplot=2
    , ylims=(0, 0.65), xlims=(0.00,1.0)
    , foreground_color_legend=false
    , right_margin=1mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm
    , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
        , xticks=(0.0:0.2:1.0,["0","0.2","0.4","0.6","0.8","1"]) 
        , yticks=(0.0:0.2:0.6,["0","0.2","0.4","0.6"]) 
)


for s in eachindex(sub_samples)
    

    pdf_params_s = DirichletPDFParams(K_u=sub_samples.v.K_u[s], K_d=sub_samples.v.K_d[s], K_q=sub_samples.v.K_q[s],
                                      λ_g1=sub_samples.v.λ_g1[s], 
                                      λ_g2=sub_samples.v.λ_g2[s],
                                      K_g=sub_samples.v.K_g[s], 
                                      λ_q=sub_samples.v.λ_q[s], 
                                      θ=Vector(sub_samples.v.θ[s]))
    
        λ_d= pdf_params_s.θ[2]*(1+pdf_params_s.K_d)/(1-pdf_params_s.θ[2])

    plot!(x_grid, [x_dv_x(x, λ_d, pdf_params_s.K_d) for x in x_grid],ylims=(0, 0.65),lw=0.5,alpha=0.3,label=:none, subplot=2)

end

    p = plot!(xlabel=L"x", ylabel=L"xd_v", subplot=2)

p = plot!(x_grid, [x_dv_x(x, λ_d_true, K_d_true) for x in x_grid], label=L"~xd_v \; \mathrm{true}",  lw=3, c=:red,
    subplot=2
    , ylims=(0, 0.65), xlims=(0.00,1.0)
    , foreground_color_legend=false
    , right_margin=1mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm
    , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
        , xticks=(0.0:0.2:1.0,["0","0.2","0.4","0.6","0.8","1"]) 
        , yticks=(0.0:0.2:0.6,["0","0.2","0.4","0.6"]) 
)

#


    p = plot!(x_grid, [0.1*x_g_x(x, pdf_params.λ_g1, pdf_params.λ_g2, pdf_params.K_g, pdf_params.K_q, pdf_params.θ[3], pdf_params.θ[4])
                       for x in x_grid], label=L"~xg/10 \; \mathrm{true}",  lw=3, c=:red,
    subplot=3
     , ylims=(0, 0.65), xlims=(0.0,1.0)
    , foreground_color_legend=false
    , right_margin=-2mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm
    , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
        , xticks=(0.0:0.2:1.0,["0","0.2","0.4","0.6","0.8","1"]) 
        , yticks=(0.0:0.2:0.6,["0","0.2","0.4","0.6"]) 
)


for s in eachindex(sub_samples)
    
   pdf_params_s = DirichletPDFParams(K_u=sub_samples.v.K_u[s], K_d=sub_samples.v.K_d[s], K_q=sub_samples.v.K_q[s],
                                      λ_g1=sub_samples.v.λ_g1[s], 
                                      λ_g2=sub_samples.v.λ_g2[s],
                                      K_g=sub_samples.v.K_g[s], 
                                      λ_q=sub_samples.v.λ_q[s], 
                                      θ=Vector(sub_samples.v.θ[s]))

    p = plot!(x_grid, [0.1*x_g_x(x, pdf_params_s.λ_g1, pdf_params_s.λ_g2, pdf_params_s.K_g, pdf_params_s.K_q, pdf_params_s.θ[3], pdf_params_s.θ[4])
                       for x in x_grid], ylims=(0, 0.65),lw=0.5,alpha=0.3,label=:none,subplot=3)
end

    p = plot!(xlabel=L"x", ylabel=L"xg/10", subplot=3)

    p = plot!(x_grid, [0.1*x_g_x(x, pdf_params.λ_g1, pdf_params.λ_g2, pdf_params.K_g, pdf_params.K_q, pdf_params.θ[3], pdf_params.θ[4])
                       for x in x_grid], label=L"~xg/10 \; \mathrm{true}",  lw=3, c=:red,
    subplot=3
     , ylims=(0, 0.65), xlims=(0.0,1.0)
    , foreground_color_legend=false
    , right_margin=-2mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm
    , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
        , xticks=(0.0:0.2:1.0,["0","0.2","0.4","0.6","0.8","1"]) 
        , yticks=(0.0:0.2:0.6,["0","0.2","0.4","0.6"]) 
)




#
    p = plot!(x_grid, [x_q_x(x, pdf_params.λ_q, pdf_params.K_q, pdf_params.θ[5]) for x in x_grid],
    ylims=(0, 0.65), xlims=(0.0,1.0), label=L"~x\bar{u} \; \mathrm{true}",  lw=3, c=:red,
    subplot=4

    , foreground_color_legend=false
    , right_margin=1mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm
    , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , xticks=(0.0:0.2:1.0,["0","0.2","0.4","0.6","0.8","1"]) 
        , yticks=(0.0:0.2:0.6,["0","0.2","0.4","0.6"]) 
)




for s in eachindex(sub_samples)
    
   pdf_params_s = DirichletPDFParams(K_u=sub_samples.v.K_u[s], K_d=sub_samples.v.K_d[s], K_q=sub_samples.v.K_q[s],
                                      λ_g1=sub_samples.v.λ_g1[s], 
                                      λ_g2=sub_samples.v.λ_g2[s],
                                      K_g=sub_samples.v.K_g[s], 
                                      λ_q=sub_samples.v.λ_q[s], 
                                      θ=Vector(sub_samples.v.θ[s]))


        p = plot!(x_grid, [x_q_x(x, pdf_params_s.λ_q, pdf_params_s.K_q, pdf_params_s.θ[5]) for x in x_grid],ylims=(0, 0.65), lw=0.5,alpha=0.3,label=:none,subplot=4)

end

    p = plot!(xlabel=L"x", ylabel=L"x\bar{u}", subplot=4)


    p = plot!(x_grid, [x_q_x(x, pdf_params.λ_q, pdf_params.K_q, pdf_params.θ[5]) for x in x_grid],
    ylims=(0, 0.65), xlims=(0.0,1.0), label=L"~x\bar{u} \; \mathrm{true}",  lw=3, c=:red,
    subplot=4

    , foreground_color_legend=false
    , right_margin=1mm
    , left_margin=0mm
    , top_margin=0mm
    , bottom_margin=-1mm
    , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , xticks=(0.0:0.2:1.0,["0","0.2","0.4","0.6","0.8","1"]) 
        , yticks=(0.0:0.2:0.6,["0","0.2","0.4","0.6"]) 
)

filename = string("figures/fig6-parton-xf(x)-", parsed_args["fitresults"], "_v2.pdf")
savefig(filename)







plot(framestyle=:axes, size=(1200, 500), fontfamily=font_family, 
    leftmargin=6Plots.mm, bottommargin=5Plots.mm, rightmargin=5Plots.mm,
    layout=@layout([a b c{0.15w}]),
    xtickfontsize=14,ytickfontsize=14,yguidefontsize=18,xguidefontsize=18
   , grid=false
)

plot!(inset=(1, bbox(0.23, 0.75, 0.55, 0.25, :bottom)))
plot!(inset=(2, bbox(0.23, 0.75, 0.55, 0.25, :bottom)))

plot!(xlim=(5e-2, 1), ylim=(0, 1000), xlabel="\$x\$", ylabel="Counts", xscale=:log, 
    grid=false, legend=:false, foreground_color_legend=nothing,
    legendfontsize=10, thickness_scaling=1, 
    xticks=([0.1, 10^-0.5, 1.0],[L"$10^{-1}$",L"$10^{-0.5}$",L"$1$"]),
    subplot=1)
plot!(xlim=(5e-2, 1), ylim=(0, 1000), xlabel="\$x\$", ylabel="", xscale=:log, 
    grid=false, legend=false, 
    xticks=([0.1, 10^-0.5, 1.0],[L"$10^{-1}$",L"$10^{-0.5}$",L"$1$"]),
    subplot=2)

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
    subplot=3, xlim=(1,2), ylim=(0, 1000), legendfontsize=14, thickness_scaling=1,
    left_margin=-10Plots.mm)
annotate!(1.4, 980, text("\$Q^2\$ [GeV\$^2\$]", 14, font_family), subplot=3)
annotate!(0.13, 600, text(L"$e^{-}$", 22, font_family), subplot=1)
p = annotate!(0.13, 600, text(L"$e^{+}$", 22, font_family), subplot=2)
p

filename = string("figures/fig7-",parsed_args["fitresults"],"_v2.pdf")
savefig(filename)
println("Done")
end

main()
