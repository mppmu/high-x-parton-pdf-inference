#!/usr/bin/julia

using BAT, DensityInterface
using PartonDensity
using QCDNUM
using Plots, Random, Distributions, ValueShapes, ParallelProcessingTools
using StatsBase, LinearAlgebra
using SpecialFunctions, Printf
using Distributions
import HDF5
using DelimitedFiles
using LaTeXStrings
using HypothesisTests
using Measures

using ArgParse
import HDF5

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

Plots.gr(format="png")
color_scheme = :viridis
default(fontfamily = "Computer Modern")
Plots.scalefontsizes()
Plots.scalefontsizes(1.2);

samples_sim = bat_read(string("fitresults/", parsed_args["fitresults"], ".h5")).result
pdf_params_gen, sim_data = pd_read_sim(string("pseudodata/", parsed_args["pseudodata"], ".h5"))


qcdnum_grid = QCDNUM.GridParams(x_min=[1.0e-3, 1.0e-1, 5.0e-1], x_weights=[1, 2, 2], nx=100, qq_bounds=[1.0e2, 3.0e4], qq_weights=[1.0, 1.0], nq=50, spline_interp=3)
qcdnum_params = QCDNUM.EvolutionParams(order=2, α_S=0.118, q0=100.0, grid_params=qcdnum_grid, n_fixed_flav=5, iqc=1, iqb=1, iqt=1, weight_type=1);
# now SPLINT and quark coefficients
splint_params = QCDNUM.SPLINTParams();
quark_coeffs = QuarkCoefficients();
# initialise QCDNUM
forward_model_init( qcdnum_params, splint_params)

nbins=sim_data["nbins"]
counts_obs_ep_sim= sim_data["counts_obs_ep"]
counts_obs_em_sim= sim_data["counts_obs_em"]



prob_ep_gen = zeros(nbins)
prob_em_gen = zeros(nbins)

prob_ep_sim = zeros(nbins)
prob_em_sim = zeros(nbins)



mode_pars_sim= mode(samples_sim)


pdfpars(params)=   DirichletPDFParams(K_u=params.K_u, K_d=params.K_d, λ_g1=params.λ_g1, λ_g2=params.λ_g2,K_g=params.K_g, λ_q=params.λ_q, K_q=params.K_q,  θ=Vector(params.θ))
        
#
# First the generated values
#
pdf_params = pdfpars(pdf_params_gen)
counts_pred_ep_gen, counts_pred_em_gen = forward_model(pdf_params, qcdnum_params, splint_params, quark_coeffs)
#
# Now the results from fitting the simulation 
#
pdf_params = pdfpars(mode_pars_sim)
counts_pred_ep_sim, counts_pred_em_sim = forward_model(pdf_params, qcdnum_params, splint_params, quark_coeffs)
#

#
# Calculate the Poisson probabilities for the different data results
#
        for j in 1:nbins
            pred=counts_pred_ep_gen[j]   
            best=floor(counts_pred_ep_gen[j])
            prob_ep_gen[j] = pdf(Poisson(pred), counts_obs_ep_sim[j])/pdf(Poisson(pred), best)
            pred=counts_pred_em_gen[j]   
            best=floor(counts_pred_em_gen[j])
            prob_em_gen[j] = pdf(Poisson(pred), counts_obs_em_sim[j])/pdf(Poisson(pred), best)

#           
            pred=counts_pred_ep_sim[j]   
            best=floor(counts_pred_ep_sim[j])
            prob_ep_sim[j] = pdf(Poisson(pred), counts_obs_ep_sim[j])/pdf(Poisson(pred), best)
            pred=counts_pred_em_sim[j]   
            best=floor(counts_pred_em_sim[j])
            prob_em_sim[j] = pdf(Poisson(pred), counts_obs_em_sim[j])/pdf(Poisson(pred), best)
#           
       end

p1=scatter(counts_obs_ep_sim,  markersize=3, markershape=:circle,label=L"$e^+$p Simulated",ylab="Counts", grid=false,
   xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
)
p1=scatter!(counts_pred_ep_gen, 
   markershape = :diamond,     markeralpha = 0.4,     markercolor = :green,  markersize = 2,
    label=L"$e^+p$ Generated", grid=false
      , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    )
p2=scatter(counts_obs_em_sim,color="red", markersize=2,markershape=:circle,label=L"$e^-p$ Simulated", grid=false,markerstrokecolor=:red
  , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
)
p2=scatter!(counts_pred_em_gen,color="red",    
    markershape = :diamond,     markeralpha = 0.4,     markercolor = :green,  markersize = 2,
    label=L"$e^-p$ Generated",
    linewidth = 2,markerstrokecolor=:red
      , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    )
p3=scatter(prob_ep_gen,color="blue",yscale=:log10,markersize=2,markershape=:circle,label="",xlab="Bin",ylab="Scaled Probability", grid=false,markerstrokecolor=:blue
  , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
)
p4=scatter(prob_em_gen,color="red",yscale=:log10, markersize=2,markershape=:circle,label="",xlab="Bin", grid=false,markerstrokecolor=:red
  , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
)
plot(p1,p2,p3,p4,layout=(2,2))

p1=scatter(counts_obs_ep_sim, color="blue",  markersize=2, markershape=:circle,label=L" $e^{+}p$ Observed",ylab="Counts", grid=false,markerstrokecolor=:blue
,ylims=(0, 900), xlims=(-1,155)
    , size=(600, 500)
      , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
,left_margin=15.0mm,bottom_margin=0.0mm, yticks=(0:200:800,["0","200","400","600","800"]), xticks=(0:50:150,["","","",""])
)
p1=plot!(counts_pred_ep_sim,color="blue",    linealpha = 0.5,
    linewidth = 2,label=L" $e^{+}p$ Predicted", grid=false,markerstrokecolor=:blue
,ylims=(0, 900), xlims=(-1,155)
  , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , foreground_color_legend=false
,left_margin=15.0mm,bottom_margin=0.0mm, yticks=(0:200:800,["0","200","400","600","800"]),xticks=(0:50:150,["","","",""])
)
p2=scatter(counts_obs_em_sim,color="red", markersize=2,markershape=:circle,label=L" $e^{-}p$ Observed", grid=false,markerstrokecolor=:red
,ylims=(0, 900), xlims=(-1,155)
  , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , foreground_color_legend=false
,left_margin=0.0mm,bottom_margin=0.0mm, yticks=(0:200:800,["","","","",""]), xticks=(0:50:150,["","","",""])
)
p2=plot!(counts_pred_em_sim,color="red",    linealpha = 0.5,label=L" $e^{-}p$ Predicted", grid=false,markerstrokecolor=:red
,linewidth = 2
,ylims=(0, 900), xlims=(-1,155)
  , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , foreground_color_legend=false
,left_margin=0.0mm,bottom_margin=0.0mm, yticks=(0:200:800,["","","","",""]), xticks=(0:50:150,["","","",""])
)
p3=scatter(prob_ep_sim,color="blue",yscale=:log10,markersize=2,markershape=:circle,label="",xlab="Bin", grid=false,markerstrokecolor=:blue
,
    ylab="Scaled probability"
,ylims=(0.003, 1.1), xlims=(-1,155)
  , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , foreground_color_legend=false
,left_margin=15.0mm,bottom_margin=9.0mm, yticks=([0.01,0.1,1.0],["0.01","0.1","1.0"]), xticks=(0:50:150,["0","50","100","150"])
)
p4=scatter(prob_em_sim,color="red",yscale=:log10, markersize=2,markershape=:circle,label="",xlab="Bin", grid=false,markerstrokecolor=:red
,ylims=(0.003, 1.1), xlims=(-1,155)
      , xtickfontsize=14,ytickfontsize=14,yguidefontsize=16,xguidefontsize=16, legendfontsize=14
    , foreground_color_legend=false
    ,left_margin=0.0mm,bottom_margin=9.0mm, yticks=([0.01,0.1,1.0],["","",""]), xticks=(0:50:150,["0","50","100","150"])
)


plot(p1,p2,p3,p4,layout=(2,2))

savefig(string("figures/fig9-Data-GoF-",parsed_args["fitresults"],".pdf"))
end
main()









