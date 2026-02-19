#include("Main/Main.jl")
#using .Harmonic_Oscillator_solver
using Plots

###############################################################################
# Physical system: two noninteracting bosons
###############################################################################

particles = [2]
system = Few_Particles_Hamonic_Oscillator(particles, 10, "bose")


###############################################################################
# STA scaling function b(t)
###############################################################################

function b_poly(s, bf)
    return 1 + (bf - 1) * (10*s^3 - 15*s^4 + 6*s^5)
end

function ddb_poly(s, bf, tf)
    return (bf - 1) * (60*s - 180*s^2 + 120*s^3) / tf^2
end

###############################################################################
# STA frequency from Ermakov equation
###############################################################################

function omega_STA(s, ω0, ωf, tf)
    bf = sqrt(ω0 / ωf)
    b  = b_poly(s, bf)
    ddb = ddb_poly(s, bf, tf)
    ω2 = (ω0^2) / b^4 - ddb / b
    return real(sqrt(complex(ω2)))   # allow expulsive regions
end

###############################################################################
# Discretization
###############################################################################

np = 2
s_vals = range(0, 1, np)

###############################################################################
# Initial and final trap frequencies
###############################################################################

ω0 = 0.5
ωf = 3.0

Ω0 = ω0^2 - 1
Ωf = ωf^2 - 1

###############################################################################
# Scan different stroke durations
###############################################################################

Tstrokes = range(1,1,1)

final_energy_STA = Float64[]
final_energy_linear = Float64[]

initial_energy = NaN
final_energy_theory = NaN

for tf in Tstrokes
    println("Running tf = $tf")

    # ---------------- STA protocol ----------------
    ωs = [omega_STA(s, ω0, ωf, tf) for s in s_vals]
    oms_STA = ωs .^ 2 .- 1.0

    # ---------------- Linear protocol -------------
    oms_linear = range(Ω0, Ωf, length=np)

    # Interaction fixed to zero
    gom = zeros(np)

    # Time grid
    ts = range(0.0, tf, length=np)

    # STA evolution
    en_STA, den_STA, norm_STA, fid_STA = time_evolution(
        system, gom, ts,
        type="breath", omegas=oms_STA
    )

    # Linear evolution
    en_lin, den_lin, norm_lin, fid_lin = time_evolution(
        system, gom, ts,
        type="breath", omegas=oms_linear
    )

    global initial_energy
    global final_energy_theory
    if isnan(initial_energy)
        initial_energy = en_STA[1]
        final_energy_theory = initial_energy * (ωf / ω0)
    end

    push!(final_energy_STA, en_STA[end])
    push!(final_energy_linear, en_lin[end])
end

print("Initial energy: $initial_energy")
print(om)