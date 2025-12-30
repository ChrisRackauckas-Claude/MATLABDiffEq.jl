using DiffEqBase, MATLABDiffEq, ParameterizedFunctions, Test
using JET

@testset "JET static analysis" begin
    @testset "buildDEStats type stability" begin
        # Test with full stats dictionary
        full_stats = Dict{String, Any}(
            "nfevals" => 100,
            "nfailed" => 5,
            "nsteps" => 95,
            "nsolves" => 50,
            "npds" => 10,
            "ndecomps" => 8
        )
        result = MATLABDiffEq.buildDEStats(full_stats)
        @test result isa DiffEqBase.Stats
        @test result.nf == 100
        @test result.nreject == 5
        @test result.naccept == 95

        # Test with empty stats dictionary (common case)
        empty_stats = Dict{String, Any}()
        result_empty = MATLABDiffEq.buildDEStats(empty_stats)
        @test result_empty isa DiffEqBase.Stats
        @test result_empty.nf == 0

        # Test with partial stats dictionary
        partial_stats = Dict{String, Any}("nfevals" => 42)
        result_partial = MATLABDiffEq.buildDEStats(partial_stats)
        @test result_partial.nf == 42
        @test result_partial.nreject == 0

        # JET analysis on buildDEStats
        @test_opt target_modules = (MATLABDiffEq,) MATLABDiffEq.buildDEStats(full_stats)
    end

    @testset "Algorithm struct instantiation" begin
        # Verify all algorithm types instantiate without issues
        @test MATLABDiffEq.ode23() isa MATLABDiffEq.MATLABAlgorithm
        @test MATLABDiffEq.ode45() isa MATLABDiffEq.MATLABAlgorithm
        @test MATLABDiffEq.ode113() isa MATLABDiffEq.MATLABAlgorithm
        @test MATLABDiffEq.ode23s() isa MATLABDiffEq.MATLABAlgorithm
        @test MATLABDiffEq.ode23t() isa MATLABDiffEq.MATLABAlgorithm
        @test MATLABDiffEq.ode23tb() isa MATLABDiffEq.MATLABAlgorithm
        @test MATLABDiffEq.ode15s() isa MATLABDiffEq.MATLABAlgorithm
        @test MATLABDiffEq.ode15i() isa MATLABDiffEq.MATLABAlgorithm
    end
end

f = @ode_def_bare LotkaVolterra begin
    dx = a * x - b * x * y
    dy = -c * y + d * x * y
end a b c d
p = [1.5, 1, 3, 1]
tspan = (0.0, 10.0)
u0 = [1.0, 1.0]
prob = ODEProblem(f, u0, tspan, p)
sol = solve(prob, MATLABDiffEq.ode45())

function lorenz(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz, u0, tspan)
sol = solve(prob, MATLABDiffEq.ode45())

algs = [MATLABDiffEq.ode23
        MATLABDiffEq.ode45
        MATLABDiffEq.ode113
        MATLABDiffEq.ode23s
        MATLABDiffEq.ode23t
        MATLABDiffEq.ode23tb
        MATLABDiffEq.ode15s]

for alg in algs
    sol = solve(prob, alg())
end
