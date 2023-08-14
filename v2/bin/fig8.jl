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

c2=colorant"#CCE5E5"
c1=colorant"#93A0AB"
color_scheme = :viridis
font_family = "Computer Modern"
default(fontfamily = "Computer Modern")
#Plots.scalefontsizes()
#Plots.scalefontsizes(1.2);
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
    counts_ep_pred_s, counts_em_pred_s = forward_model(pdf_params_s,qcdnum_params, splint_params, quark_coeffs, sys_err_params)
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
pdfpars(params)=   DirichletPDFParams( K_u=params.K_u, K_d=params.K_d, λ_g1=params.λ_g1, λ_g2=params.λ_g2, K_q=params.K_q,
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



p1=histogram(chisqep,bins=100,xlabel=L"\chi^2_P", ylabel="Entries", fontfamily=font_family,color=c1, linecolor=c1, grid=false)
p1=plot!([chisqep_data],seriestype = :vline,lw=5,legend=:none, fontfamily=font_family 
, xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
        ,ylims=(0, 500), xlims=(65,260)
        ,color=:red, grid=false,left_margin=14mm,bottom_margin=5.5mm
)
p2=histogram(chisqem,bins=50:2:250,legend=:false,xlabel=L"\chi^2_P", ylabel="Entries", fontfamily=font_family,color=c1, linecolor=c1, grid=false)

p2=plot!([chisqem_data],seriestype = :vline,lw=5
, xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    ,ylims=(0, 500), xlims=(65,260)
,color=:red, grid=false,left_margin=14mm,bottom_margin=5.5mm
)
annotate!(p1,220.0,350,text(L"$e^{+}p$",26))
annotate!(p2,220.0,350,text(L"$e^{-}p$",26))
#plot(p1,p2,layout=(2,1))

plot(p1,p2,layout=(2,1),#size=(800,400),
    #top_margin=-3mm,
    #bottom_margin=-3.5mm,
    #right_margin=-4mm,
    left_margin=1mm
)

filename = string("figures/fig8-chisq-pvalue-",parsed_args["fitresults"],"_v2.pdf")
savefig(filename)
end

main()

