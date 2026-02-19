###############################################################################
# STA protocol for trap-frequency stroke (g = 0)
# Implements inverse engineering using Ermakov equation
###############################################################################

include("Main/Main.jl")
using .Harmonic_Oscillator_solver
using Plots

###############################################################################
# Physical system: two bosons, noninteracting
###############################################################################

particles = [2]
system = Few_Particles_Hamonic_Oscillator(
    particles,
    10,
    "bose"
)

###############################################################################
# STA scaling function b(t)
###############################################################################

function b_poly(s, bf)
    return 1 + (bf - 1) * (10*s^3 - 15*s^4 + 6*s^5)
end

function db_poly(s, bf, tf)
    return (bf - 1) * (30*s^2 - 60*s^3 + 30*s^4) / tf
end

function ddb_poly(s, bf, tf)
    return (bf - 1) * (60*s - 180*s^2 + 120*s^3) / tf^2
end

###############################################################################
# STA frequency from b(t)
###############################################################################

function omega_STA(s, ω0, ωf, tf)
    bf = sqrt(ω0 / ωf)
    b  = b_poly(s, bf)
    ddb = ddb_poly(s, bf, tf)

    ω2 = (ω0^2) / b^4 - ddb / b
    return real(sqrt(complex(ω2)))  # allow expulsive regions
end

###############################################################################
# Discretization
###############################################################################

np = 101                 # number of time points per stroke
s_vals = range(0, 1, np)

###############################################################################
# Initial and final frequencies
###############################################################################

ω0 = 1.0
ωf = 8.0

###############################################################################
# Scan different stroke durations
###############################################################################

Tstrokes = range(50.0, 2000.0, 5)
final_energy = Float64[]
initial_energy = NaN
final_energy_theory = NaN

for tf in Tstrokes
    println("Running STA stroke with tf = $tf")

    # Build STA frequency ramp
    oms_in = [omega_STA(s, ω0, ωf, tf) for s in s_vals]
    oms = oms_in .^ 2 .- 1.0   # omega(t)^2 - 1

    # Interaction strength fixed to zero
    gom = zeros(np)

    # Time grid
    ts = range(0.0, tf, length=np)

    enert, denst, normt, fidelityt = time_evolution(
        system, gom, ts,
        type="breath", omegas=oms
    )

    global initial_energy   # tell Julia you mean the global one
    global final_energy_theory
    if isnan(initial_energy) 
        initial_energy = enert[1]
        final_energy_theory = initial_energy * (ωf / ω0)
    end


    push!(final_energy, enert[end])
end

###############################################################################
# Plot final energy vs stroke duration
###############################################################################

p = plot(Tstrokes, final_energy,
         marker=:circle,
         label="Final energy (STA)",
         xlabel="Stroke duration",
         ylabel="Energy")

hline!(p, [final_energy_theory],
       linestyle=:dash,
       label="Theoretical final energy")

savefig("sta_trap_frequency_scan.pdf")
display(plot!())
