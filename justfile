default:
    just --list

format:
    julia --project=@JuliaFormatter -e 'using JuliaFormatter; format(".")'

docs:
    julia --project=docs docs/make.jl

servedocs:
    julia --project=docs -e 'using LiveServer; serve(dir="docs/build")'

servedocs-continuously:
    julia --project=docs -e 'using ProbNumDiffEq, LiveServer; servedocs()'
