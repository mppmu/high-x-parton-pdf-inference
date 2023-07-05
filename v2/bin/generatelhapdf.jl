#!/usr/bin/julia
using BAT, DensityInterface
using PartonDensity
using QCDNUM
using Plots, Random, Distributions, ValueShapes, ParallelProcessingTools
using StatsBase, LinearAlgebra
using Printf
using ArgParse

function func(ipdf, x)::Float64
pdf_params = DirichletPDFParams(K_u=3.7, K_d=3.7, λ_g1=0.5, λ_g2=-0.5, K_g=5.0,λ_q=-0.5, K_q=6.0,θ=[0.22, 0.10, 0.24, 0.24, 0.10,0.05, 0.01, 0.005, 0.0005])
r::Float64 = get_input_pdf_func(pdf_params)(ipdf, x)
return r
end

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
        "--fitresults", "-f"
            help = "Input fitresults -- file in the pseudodata directory w/o the extension"
            arg_type = String
            default = "fit-Dirichlet-0-45-simulation-Dirichlet-42"    
            
            
    end

    return parse_args(s)
end


function main()
parsed_args = parse_commandline()
println("Parsed args:")
for (arg,val) in parsed_args
  println("  $arg  =>  $val")
end
seed=parsed_args["seed"]
println(seed)
seedtxt=string(seed)
rng = MersenneTwister(seed)

if parsed_args["parametrisation"] == "Dirichlet"
  θ = [ 0.228, 0.104, 0.249, 0.249, 0.104, 0.052, 0.010, 0.005, 0.0005]
  θ_sum=sum(θ[1:9])
  θ=θ/θ_sum  
  pdf_params0 = DirichletPDFParams(K_u=3.7, K_d=3.7, λ_g1=0.5, λ_g2=-0.5, K_g=5.0,λ_q=-0.5, K_q=6.0, θ=θ);
end
if parsed_args["parametrisation"] == "Valence"
#  weights = [5.0, 5.0, 1.0, 1.0, 1.0, 0.5, 0.5]
#  θ = get_θ_val(rng, λ_u, K_u, λ_d, K_d, weights)
  λ_u = 0.64;
  K_u = 3.38;
  λ_d = 0.67;
  K_d = 4.73;
  θ=[0.22, 0.10, 0.24, 0.24, 0.10,0.05, 0.01, 0.005, 0.0005]
  θ_sum=sum(θ[1:9])
  θ=θ/θ_sum
  pdf_params1 = ValencePDFParams(λ_u=λ_u, K_u=K_u, λ_d=λ_d, K_d=K_d, λ_g1=0.50, λ_g2=-0.63, K_g=4.23, λ_q=-0.23, K_q=5.0, θ=θ);
end


if parsed_args["parametrisation"] == "Bernstein"
  θ=[.33, .13, .27, .17, .073, 0.014, 0.002, 0.000001, .003]
  θ_sum=sum(θ[1:9])
  θ=θ/θ_sum  
  bspoly_params = [1,4,0,4,0,5];
  λ_g1=1.5;
  λ_g2=-0.4;
  K_g=6.0;
  λ_q=-0.25;
  K_q=5.0;
  initial_U = [-8.];
  initial_D = [15.0];
  vec_bspp = Vector(bspoly_params)
  bspoly_params = [[vec_bspp[Int(2 * i - 1)], vec_bspp[Int(2 * i)]] for i in 1:length(vec_bspp)/2]
  bspoly_params_d = 0
  try
    vec_bsppd = Vector(bspoly_params_d)
    bspoly_params_d = [[vec_bsppd[Int(2 * i - 1)], vec_bsppd[Int(2 * i)]] for i in 1:length(vec_bsppd)/2]
  catch err
    bspoly_params_d = bspoly_params
  end 
  initU = Vector(initial_U)
  initD = Vector(initial_D)       
  pdf_params2 = BernsteinDirichletPDFParams(initial_U=initU,
                    initial_D=initD,
                    λ_g1=λ_g1, λ_g2=λ_g2,
                    K_g=K_g, λ_q=λ_q, K_q=K_q,
                    bspoly_params=bspoly_params,
                    bspoly_params_d=bspoly_params_d,
                    θ=Vector(θ)
  )
end


qcdnum_grid = QCDNUM.GridParams(x_min=[1.0e-3, 1.0e-1, 5.0e-1], x_weights=[1, 2, 2], nx=100, qq_bounds=[1.0e2, 3.0e4], qq_weights=[1.0, 1.0], nq=50, spline_interp=3)
qcdnum_params = QCDNUM.EvolutionParams(order=2, α_S=0.118, q0=100.0, grid_params=qcdnum_grid,n_fixed_flav=5, iqc=1, iqb=1, iqt=1, weight_type=1);
    
splint_params = QCDNUM.SPLINTParams();
quark_coeffs = QuarkCoefficients();
forward_model_init(qcdnum_params, splint_params)

func_c = @cfunction(func, Float64, (Ref{Int32}, Ref{Float64}))


xmin = Float64.([1.0e-4])
iwt = Int32.([1])
ng = 1
nxin = 100
iosp = 3
nx = 10

qq = Float64.([2e0, 1e4])
wt = Float64.([1e0, 1e0])
nq = 1
nqin = 60
ngq = 2
itype = 1

as0 = 0.364
r20 = 2.0

q2c = 3.0
q2b = 25.0
q0 = 2.0
iqt = 999

def = Float64.([0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0.,     
                0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0.,      
                0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0.,      
                0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,      
                0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,      
                0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      
                0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,      
                0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      
                0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      
                0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      
                0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      
                0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]);

nfin = 0
x = 1.0e-3
q = 1.0e3
pdf = Array{Float64}(undef, 13)
qmz2 = 8315.25;

QCDNUM.qcinit(-6, " ")
nx = QCDNUM.gxmake(xmin, iwt, ng, nxin, iosp)
nq = QCDNUM.gqmake(qq, wt, ngq, nqin)
nw = QCDNUM.fillwt(itype)
QCDNUM.setord(3)
QCDNUM.setalf(as0, r20)
iqc = QCDNUM.iqfrmq(q2c)
iqb = QCDNUM.iqfrmq(q2b)
iqt = QCDNUM.iqfrmq(1e11)
QCDNUM.setcbt(0, iqc, iqb, 0)
iq0 = QCDNUM.iqfrmq(q0)
#println("OK1")
QCDNUM.evolfg(itype, func_c, def, iq0)

allx = Float64.([1.0e-3,2.0e-3,3.0e-3,1.0e-2,2.0e-2,1.0e-1,1.0])
allq = Float64.([100.,200,300,400,500,1000])
allalpha = Float64.([100.,200,300,400,500,1000])
for i in 1:6
allalpha[i] =QCDNUM.asfunc(allq[i]*allq[i])[1]
end

open("CABCHSV.info", "w") do f
write(f,"
SetDesc: \"CABCHSV PDF fixed 5-flavour fits\"
SetIndex: 11082
Authors: Francesca Capel, Ritu Aggarwal, Michiel Botje, Allen Caldwell, Richard Hildebrandt, Oliver Schulz, Andrii Verbytskyi
Reference: arXiv:2309.xxxx
Format: lhagrid1
DataVersion: 1
NumMembers: 2
Particle: 2212
Flavors: [-5,-4, -3, -2, -1, 1, 2, 3, 4, 5, 21]
OrderQCD: 2
ForcePositive: 1
FlavorScheme: fixed
NumFlavors: 5
XMin: 0.001
XMax: 1
QMin: 100
QMax: 100000
MZ: 91.1876
MUp: 0
MDown: 0
MStrange: 0
MCharm: 1.3
MBottom: 4.75
MTop: 172
AlphaS_MZ: 0.118
AlphaS_OrderQCD: 2
AlphaS_Type: ipol
AlphaS_Qs: ")
write(f,allq)
write(f,"\n")
write(f,"AlphaS_Vals: ")
write(f,allalpha)
write(f,"\n")
write(f,"AlphaS_Lambda4: 0.326
AlphaS_Lambda5: 0.226
\n")

end





open("CABCHSV_0000.dat", "w") do f
write(f,"PdfType: central
Format: lhagrid1
---")
write(f,"\n")
write(f,["$v " for v in allx]...) 
write(f,"\n")
write(f,["$v " for v in allq]...) 
write(f,"\n")
write(f,"-5 -4 -3 -2 -1 1 2 3 4 5 21")
write(f,"\n")
for j in 1:7
for i in 1:6
pdf = QCDNUM.allfxq(itype, allx[j], allq[i], 0, 1)
 write(f,["$v " for v in pdf]...) 
 write(f,"\n")
end 
end

end







Ns = 50
NN=0
#Fit results!!!
samples_data = bat_read(string("fitresults/", parsed_args["fitresults"], ".h5")).result;
sub_samples = BAT.bat_sample(samples_data, BAT.OrderedResampling(nsamples=Ns)).result;
for s in eachindex(sub_samples)

    NN=NN+1
    pdf_params_s = DirichletPDFParams(K_u=sub_samples.v.K_u[s], K_d=sub_samples.v.K_d[s], K_q=sub_samples.v.K_q[s],
                                      λ_g1=sub_samples.v.λ_g1[s], 
                                      λ_g2=sub_samples.v.λ_g2[s],
                                      K_g=sub_samples.v.K_g[s], 
                                      λ_q=sub_samples.v.λ_q[s], 
                                      θ=Vector(sub_samples.v.θ[s]))
    λ_u= pdf_params_s.θ[1]*(1+pdf_params_s.K_u)/(2-pdf_params_s.θ[1])
    
    pdf_params = pdf_params_s
    QCDNUM.evolfg(itype, func_c, def, iq0)


open(string("CABCHSV_", lpad(string(NN),3,"0"),".dat"), "w") do f
write(f,"PdfType: error
Format: lhagrid1
---")
write(f,"\n")
write(f,["$v " for v in allx]...) 
write(f,"\n")
write(f,["$v " for v in allq]...) 
write(f,"\n")
write(f,"-5 -4 -3 -2 -1 1 2 3 4 5 21")
write(f,"\n")
for j in 1:7
for i in 1:6
pdf = QCDNUM.allfxq(itype, allx[j], allq[i], 0, 1)
 write(f,["$v " for v in pdf]...) 
 write(f,"\n")
end 
end

end

    
end














end

main()
