###############################################################################
# Load project code and required modules
###############################################################################

include("Main/Main.jl")
using .Harmonic_Oscillator_solver
using Plots

###############################################################################
# Define the physical system
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
# Control protocol (shape only, time-independent)
###############################################################################

# Stroke 1: increase interaction
oms1 = range(0.0, 0.0, np)
gom1 = range(0.0, 15.0, np)

# Stroke 2: increase trap frequency
oms2 = range(0.0, 8.0, np)
gom2 = range(15.0, 15.0, np)

# Stroke 3: decrease interaction
oms3 = range(8.0, 8.0, np)
gom3 = range(15.0, 0.0, np)

# Stroke 4: decrease trap frequency
oms4 = range(8.0, 0.0, np)
gom4 = range(0.0, 0.0, np)

# Full protocol
gom = [gom1; gom2; gom3; gom4]
oms = [oms1; oms2; oms3; oms4]

###############################################################################
# Scan different total cycle times
###############################################################################

Tcycles = range(2000.0, 20000.0, 30)
final_energy = Float64[]

initial_energy = NaN   # define it first
final_energy = Float64[]

for Tcycle in Tcycles
    println("Running Tcycle = $Tcycle")

    ts = range(0.0, Tcycle, length(gom))

    enert, denst, normt, fidelityt = time_evolution(
        system, gom, ts,
        type="breath", omegas=oms
    )

    global initial_energy   # tell Julia you mean the global one
    if isnan(initial_energy)
        initial_energy = enert[1]
    end

    push!(final_energy, enert[end])
end

###############################################################################
# Plot adiabaticity test: final energy vs cycle time
###############################################################################

p = plot(Tcycles, final_energy, marker=:o, label="Final energy")
hline!(p, [initial_energy], linestyle=:dash, label="Initial energy")

plot!(p, xlabel="Total cycle time",
         ylabel="Energy",)

savefig("Figs/adiabaticity_scan_2.pdf")
display(plot!())