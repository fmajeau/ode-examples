# Set up the package environment
cd(@__DIR__)      # go to the directory of this script
using Pkg         # use the package manager
Pkg.activate(".") # activate the environment defined by the toml files in this directory
Pkg.instantiate() # install missing dependencies, make sure environment is ready to use

# Import packages
using DifferentialEquations # ODEProblem, solve
using Plots # plot

# Define parameters
m = 20  # [kg]
k = 15  # [kg-m/s] spring constant
b = 2  # [kg-m] damper constant
F = 1  # [kg-m/s^2] magnitude of constant external force
p = [m, k, b, F]

# Define time range to solve
t_start = 0  # [s]
t_end = 100  # [s]
tspan = (t_start, t_end) # [s]

# Define initial conditions
x0 = 2 # [m]
v0 = 3  # [m/s]
u0 = [x0, v0]

# Define the function that stores our  system of equations.
# It will be called within every iteration of the numerical method performed by "solve"
function mass_spring_damper!(du,u,p,t)

    # Rename states for clarity
    x = u[1]  # [m]
    v = u[2]  # [m/s]

    # Rename parameters for clarity
    m = p[1]  # [kg]
    k = p[2]  # [kg-m/s]
    b = p[3]  # [kg-m]
    F = p[4]  # [kg-m/s^2]

    # Define the differential equations
    du[1] = v
    du[2] = -((k/m)*x) - ((b/m)*v) + ((1/m)*F)

    # Print to get a sense of how this function is used by solve!
    print("[time $t] du:")
    println(du)

end

# Build problems
prob = ODEProblem(mass_spring_damper!, u0, tspan, p)

# Solve the system
solution = solve(prob, RK4())

# Plot the time-domain trajectory of all state variables on the same set of axes
plot(solution, 
    title="Time Domain: All states variables on one set of axes", 
    label=["x" "v"]
    )

# Plot the time-domain trajectory of the state variables on subplots
plot(solution, 
    layout = (2,1), 
    label=["x" "v"], 
    title=["Time Domain: Position" "Time Domain: Velocity"], 
    xlabel= "Time [s]", 
    ylabel = ["Position [m]" "Velocity [m/s]"], 
    color=["blue" "red"]
    )

# Plot the phase portrait of the two state variables
plot(solution, 
    idxs=(1,2), 
    title="Phase Portrait: Velocity vs Position", 
    xlabel="Position [m]", 
    ylabel="Velocity [m/s]"
    )