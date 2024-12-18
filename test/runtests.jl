using ChaoticDynamicalSystemLibrary
using Test
using SafeTestsets
using Aqua, JET

@testset "ChaoticDynamicalSystemLibrary.jl" begin
    # Write your tests here.

    include("test_chaotic_attractors.jl")

    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(ChaoticDynamicalSystemLibrary;
            ambiguities = (imported = false))
    end

    if VERSION >= v"1.7"
        @testset "Code linting (JET.jl)" begin
            JET.test_package(ChaoticDynamicalSystemLibrary; target_defined_modules = true)
        end
    end
end
