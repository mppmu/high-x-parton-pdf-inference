FROM fedora:39
ENV JULIA_PKGDIR=/opt/julia
ENV JULIA_DEPOT_PATH=/opt/julia   
#    export JULIA_DEPOT_PATH=/opt/myjulia   
RUN    dnf -y update && dnf -y install dnf5  && dnf -y clean all 
RUN    dnf5 -y install julia python gcc-c++ gcc-gfortran && dnf5 -y clean all 
RUN  JULIA_DEPOT_PATH=/opt/julia JULIA_PKGDIR=/opt/julia julia  -e 'import Pkg;     \
Pkg.add("SpecialFunctions");                           \
Pkg.add(url="https://github.com/bat/BAT.jl.git");      \
Pkg.add(url="https://github.com/cescalara/PartonDensity.jl.git",rev="Dec");     \
Pkg.add("ArgCheck");                                     \
Pkg.add("Colors");                                     \
Pkg.add("Plots");                                      \
Pkg.add("HDF5"); \
Pkg.add("PackageCompiler"); \
Pkg.add("ArgParse"); \
Pkg.add("LaTeXStrings"); \
Pkg.add("DensityInterface"); \
Pkg.add("QCDNUM"); \
Pkg.add("Plots"); \
Pkg.add("Random"); \
Pkg.add("Distributions"); \
Pkg.add("ValueShapes"); \
Pkg.add("ParallelProcessingTools"); \
Pkg.add("StatsBase"); \
Pkg.add("LinearAlgebra"); \
Pkg.add("SpecialFunctions"); \
Pkg.add("Printf"); \
Pkg.add("DelimitedFiles"); \
Pkg.add("LaTeXStrings"); \
Pkg.add("HypothesisTests"); \
Pkg.add("Statistics"); \
Pkg.add("Measures"); \
Pkg.add("WorkerUtilities"); \
Pkg.add("PooledArrays");   \
Pkg.add("FilePathsBase");  \
Pkg.add("SentinelArrays"); \
Pkg.add("WeakRefStrings"); \
Pkg.add("InlineStrings");  \
Pkg.add("Documenter");  \
Pkg.add("CSV");            \
Pkg.add("ArgParse");'
