# Set up the package environment
cd(@__DIR__)
using Pkg
Pkg.activate(".")
Pkg.instantiate()

# Import packages
using DifferentialEquations
using Plots

# Define parameters
l = 1.0                  # length [m]
m = 1.0                  # mass [kg]
g = 9.81                 # gravitational acceleration [m/s^2]
b = 1.0                  # damping coefficient [kg * m^2 * s^-1]
a = 0.1                  # driving torque amplitude [Nm]
Ω = 5.0                  # driving torque frequency [rad/s]
M = t -> a * sin(Ω * t)  # driving torque [Nm]
# NOTE: The driving term is defined using an anonymous function here for compactness, but a 
# function name is still required so that the function can be evaluated at each timestep within the
# ODEFunction. Alternatively, the driving term could be defined as:
# function M(t)
#     return a * sin(Ω * t)
# end
p = [l, m, g, b, M]

# Define time range to solve
t_start = 0.0  # [s]
t_end = 10.0   # [s]
tspan = (t_start, t_end)

# Define initial conditions
θ_0 = deg2rad(45)   # initial angular deflection [rad]
ω_0 = 0.0           # initial angular velocity [rad/s]
u0 = [θ_0, ω_0]

# Define the function that stores our system of equations.
function pendulum!(du, u, p, t)

    # Rename locally for clarity
    θ = u[1]
    ω = u[2]
    l = p[1]
    m = p[2]
    g = p[3]
    b = p[4]
    M = p[5]

    # Define the differential equations
    #  du[1] = dθ/dt = ω
    #  du[2] = dω/dt = -3g/(2l) sin[θ(t)] - 3b/(ml^2)ω(t) + 3/(ml^2)M(t)
    du[1] = ω
    du[2] = -3g / (2l) * sin(θ) - 3b / (m * l^2) * ω + 3 / (m * l^2) * M(t)

end

# Build the ODE problem
prob = ODEProblem(pendulum!, u0, tspan, p)

# Solve the ODE problem by specifying a numerical method
# NOTE: by passing `save_everystep=true`, it is possible to animate the solution with 
# `animate(solution, idxs=(1,2))` (or any other indices you're interested in...)
solution = solve(prob, RK4(), save_everystep=true)

# Plot the time-domain trajectory of all state variables on the same set of axes
p1 = plot(solution,
    title="Time Domain: All states variables on one set of axes",
    label=["θ [rad]" "ω [rad/s]"]
)

# Add a trace to the plot for the driving torque
plot!(M, label="Driving Torque [Nm]")

# Plot the time-domain trajectory of the state variables on subplots
p2 = plot(solution,
    layout=(2, 1),
    label=["θ [rad]" "ω [rad/s]"],
    title=["Time Domain: Angular Deflection" "Time Domain: Angular Velocity"],
    xlabel="Time [s]",
    ylabel=["Angular Deflection" "Angular Velocity"],
    color=["blue" "red"]
)

# Plot the phase portrait of the two state variables
p3 = plot(solution,
    idxs=(1, 2),
    title="Phase Portrait: Angular Velocity vs Deflection",
    label="(θ [rad], ω [rad/s])",
    xlabel="Angular Deflection [rad]",
    ylabel="Angular Velocity [rad/s]"
)