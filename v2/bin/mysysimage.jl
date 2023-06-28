#!/usr/bin/julia

using PackageCompiler
function main()
create_sysimage(["Plots","ArgParse","LaTeXStrings", "PartonDensity","QCDNUM"], sysimage_path="sys_forplots2.so", precompile_execution_file="bin/precompile.jl")


end
main()









