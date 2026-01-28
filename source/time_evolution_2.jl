###############################################################################
# Load project code and required modules
###############################################################################

include("Main/Main.jl")                        # Load main project definitions
using .Harmonic_Oscillator_solver              # Import harmonic oscillator solver
using Plots                                    # Plotting library


###############################################################################
# Define the physical system
###############################################################################

particles = [2]                                # Number of particles (2 identical particles)
system = Few_Particles_Hamonic_Oscillator(
    particles,                                # Particle configuration
    10,                                       # Number of single-particle HO basis states
    "bose"                                    # Bosonic statistics
)


###############################################################################
# Define discretization parameters
###############################################################################

np = 101                                       # Number of time steps per stroke
gs = range(0.0,15.0,np)                        # Interaction strength grid (not directly used)
ts = range(0.0,10000.0,np)                     # Time grid for one stroke


###############################################################################
# Define Pauli-engine cycle: forward direction
# Control parameters:
#   gom(t)  = interaction strength g(t)
#   oms(t)  = omega(t)^2 - 1  (trap frequency modulation)
###############################################################################

# --- Stroke 1: Increase interaction at fixed trap frequency---
oms1 = range(0.0,0.0,np)                       # Trap frequency fixed
gom1 = range(0.0,15.0,np)                      # Interaction ramped up

# --- Stroke 2: Increase trap frequency at fixed interaction ---
oms2 = range(0.0,8.0,np)                       # Trap frequency increased
gom2 = range(15.0,15.0,np)                     # Interaction fixed

# --- Stroke 3: Decrease interaction at fixed trap frequency ---
oms3 = range(8.0,8.0,np)                       # Trap frequency fixed
gom3 = range(15.0,0.0,np)                      # Interaction ramped down

# --- Stroke 4: Decrease trap frequency at zero interaction ---
oms4 = range(8.0,0.0,np)                       # Trap frequency returned
gom4 = range(0.0,0.0,np)                       # Interaction fixed at zero


###############################################################################
# Full-cycle time grid and concatenated control parameters
###############################################################################

ts4 = range(0.0,40000.0,401)                   # Time grid for full 4-stroke cycle

gom = [gom1; gom2; gom3; gom4]                 # Full interaction protocol g(t)
oms = [oms1; oms2; oms3; oms4]                 # Full trap protocol omega(t)^2 - 1


###############################################################################
# Time evolution: individual strokes
###############################################################################

# type="breath": include a time-dependent harmonic trap term (Hpot) so that 
# modulating omegas(t) drives collective breathing-mode dynamics

ener1, dens1, norm1 = time_evolution(
    system, gom1, ts,
    type="breath", omegas=oms1
)

ener2, dens2, norm2 = time_evolution(
    system, gom2, ts,
    type="breath", omegas=oms2
)

ener3, dens3, norm3 = time_evolution(
    system, gom3, ts,
    type="breath", omegas=oms3
)

ener4, dens4, norm4 = time_evolution(
    system, gom4, ts,
    type="breath", omegas=oms4
)


###############################################################################
# Time evolution: full continuous cycle
###############################################################################

enert, denst, normt = time_evolution(
    system, gom, ts4,
    type="breath", omegas=oms
)


###############################################################################
# Concatenate piecewise results for comparison
###############################################################################

tsf  = range(0.0,40000.0,4*np)                 # Time grid for stitched trajectory
ener = [ener1; ener2; ener3; ener4]            # Energy vs time
dens = [dens1; dens2; dens3; dens4]            # Density vs space and time
norm = [norm1; norm2; norm3; norm4]            # Norm vs time


###############################################################################
# Energy comparison: piecewise vs continuous evolution
###############################################################################

plot(tsf, ener, label="")
plot!(ts4, enert, label="")
plot!(xlabel="t", ylabel="E")
display(plot!())
savefig("Figs/energy_vs_t.pdf")


###############################################################################
# Norm conservation check
###############################################################################

plot(tsf, norm, label="")
plot!(ts4, normt, label="")
plot!(xlabel="t", ylabel="Norm")
display(plot!())
savefig("Figs/norm_vs_t.pdf")


###############################################################################
# Density evolution visualization
###############################################################################

xs = collect(range(-5,5,1001))                 # Spatial grid

heatmap(xs, tsf, dens, clim=(0,2.0))            # Density for stitched evolution
plot!(xlabel="x", ylabel="t")
display(plot!())
savefig("Figs/density_vs_t.pdf")

heatmap(xs, ts4, denst, clim=(0.0,2.0))         # Density for continuous evolution
plot!(xlabel="x", ylabel="t")
display(plot!())
savefig("Figs/density_vs_t_continuous.pdf")


###############################################################################
# Reverse Pauli-engine cycle (path-dependence test)
###############################################################################

# --- Stroke 1 ---
oms1 = range(0.0,8.0,np)
gom1 = range(0.0,0.0,np)

# --- Stroke 2 ---
oms2 = range(8.0,8.0,np)
gom2 = range(0.0,15.0,np)

# --- Stroke 3 ---
oms3 = range(8.0,0.0,np)
gom3 = range(15.0,15.0,np)

# --- Stroke 4 ---
oms4 = range(0.0,0.0,np)
gom4 = range(15.0,0.0,np)


###############################################################################
# Time evolution for reversed cycle
###############################################################################

ener1, dens1, norm1 = time_evolution(system, gom1, ts, type="breath", omegas=oms1)
ener2, dens2, norm2 = time_evolution(system, gom2, ts, type="breath", omegas=oms2)
ener3, dens3, norm3 = time_evolution(system, gom3, ts, type="breath", omegas=oms3)
ener4, dens4, norm4 = time_evolution(system, gom4, ts, type="breath", omegas=oms4)


###############################################################################
# Concatenate reversed results
###############################################################################

tsfv2  = range(0.0,40000.0,4*np)
enerv2 = [ener1; ener2; ener3; ener4]
densv2 = [dens1; dens2; dens3; dens4]
normv2 = [norm1; norm2; norm3; norm4]


###############################################################################
# Diagnostics for reversed cycle
###############################################################################

plot(tsfv2, enerv2, label="")
plot!(xlabel="t", ylabel="E")
display(plot!())
savefig("Figs/energy_vs_t_reverse.pdf")

plot(tsfv2, normv2, label="")
plot!(xlabel="t", ylabel="Norm")
display(plot!())
savefig("Figs/norm_vs_t_reverse.pdf")

heatmap(xs, tsfv2, densv2)
plot!(xlabel="x", ylabel="t")
display(plot!())
savefig("Figs/density_vs_t_reverse.pdf")
