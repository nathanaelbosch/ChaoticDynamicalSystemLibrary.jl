using ChaoticDynamicalSystemLibrary
using Test
using Aqua, JET

@testset "ChaoticDynamicalSystemLibrary.jl" begin
    # Write your tests here.

    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(
            ChaoticDynamicalSystemLibrary;
            ambiguities=(imported=false),
        )
    end

    @testset "Code linting (JET.jl)" begin
        JET.test_package(ChaoticDynamicalSystemLibrary; target_defined_modules = true)
    end
end
