#!/usr/bin/julia

using PackageCompiler
function main()
create_sysimage(["Plots","ArgParse","LaTeXStrings"], sysimage_path="sys_forplots.so", precompile_execution_file="bin/precompile.jl")


end
main()









