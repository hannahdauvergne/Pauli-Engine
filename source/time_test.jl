###############################################################################
# Adiabaticity vs runtime sweep
# Uses ONLY time_evolution (no Hamiltonian, no eigenstates)
###############################################################################

include("Main/Main.jl")
using .Harmonic_Oscillator_solver
using Plots
using Printf
using Dates


###############################################################################
# System definition
###############################################################################

particles = [2]
system = Few_Particles_Hamonic_Oscillator(
    particles,
    10,
    "bose"
)


###############################################################################
# Discretization
###############################################################################

np = 101   # points per stroke


###############################################################################
# Forward Pauli-engine protocol
###############################################################################

function forward_protocol(np)
    oms1 = range(0.0, 0.0, np); gom1 = range(0.0, 15.0, np)
    oms2 = range(0.0, 8.0, np); gom2 = range(15.0, 15.0, np)
    oms3 = range(8.0, 8.0, np); gom3 = range(15.0, 0.0, np)
    oms4 = range(8.0, 0.0, np); gom4 = range(0.0, 0.0, np)

    gom = [gom1; gom2; gom3; gom4]
    oms = [oms1; oms2; oms3; oms4]

    return gom, oms
end


###############################################################################
# Reverse Pauli-engine protocol
###############################################################################

function reverse_protocol(np)
    oms1 = range(0.0, 8.0, np); gom1 = range(0.0, 0.0, np)
    oms2 = range(8.0, 8.0, np); gom2 = range(0.0, 15.0, np)
    oms3 = range(8.0, 0.0, np); gom3 = range(15.0, 15.0, np)
    oms4 = range(0.0, 0.0, np); gom4 = range(15.0, 0.0, np)

    gom = [gom1; gom2; gom3; gom4]
    oms = [oms1; oms2; oms3; oms4]

    return gom, oms
end


###############################################################################
# Time grid
###############################################################################

time_grid(T, np) = range(0.0, T, 4*np)


###############################################################################
# Runtime sweep
###############################################################################

Ts = [1e3, 2e3, 5e3, 1e4, 2e4, 5e4]

ΔE_irrev = Float64[]
walltime = Float64[]

gom_f, oms_f = forward_protocol(np)
gom_r, oms_r = reverse_protocol(np)

println("Starting sweep at ", now())
println("------------------------------------------------------")

for T in Ts
    ts = time_grid(T, np)

    t0 = time()

    ener_f, _, _ = time_evolution(
        system,
        gom_f,
        ts,
        type="breath",
        omegas=oms_f
    )

    ener_r, _, _ = time_evolution(
        system,
        gom_r,
        ts,
        type="breath",
        omegas=oms_r
    )

    t1 = time()

    ΔE = abs(ener_f[end] - ener_r[end])

    push!(ΔE_irrev, ΔE)
    push!(walltime, t1 - t0)

    @printf(
        "T = %8.2e | ΔE_irrev = %.5e | runtime = %.2f s\n",
        T, ΔE, t1 - t0
    )
end


###############################################################################
# Plots
###############################################################################

plot(
    Ts,
    ΔE_irrev,
    xscale=:log10,
    yscale=:log10,
    marker=:o,
    xlabel="Total runtime T",
    ylabel="Cycle irreversibility |ΔE|",
    title="Adiabaticity vs runtime",
    legend=false
)
savefig("Figs/adiabaticity_vs_runtime.pdf")

plot(
    Ts,
    walltime,
    xscale=:log10,
    marker=:o,
    xlabel="Total runtime T",
    ylabel="Wall-clock time (s)",
    title="Computational cost",
    legend=false
)
savefig("Figs/runtime_cost.pdf")


println("------------------------------------------------------")
println("Sweep finished. Results saved in Figs/")
