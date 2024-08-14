# Set up the package environment
cd(@__DIR__)      # go to the directory of this script
using Pkg         # use the package manager
Pkg.activate(".") # activate the environment defined by the toml files in this directory
Pkg.instantiate() # install missing dependencies, make sure environment is ready to use

# Import packages
using DifferentialEquations # ODEProblem(), solve()
using Plots                 # plot()

# Define parameters
m = 20  # [kg]
k = 15  # [kg-m/s] spring constant
b = 5   # [kg-m] damper constant
F = 1   # [kg-m/s^2] magnitude of constant external force
p = [m, k, b, F] # NOTE: passed into ODEProblem; order persists in the ODEFunction

# Define time range to solve
t_start = 0  # [s]
t_end = 50   # [s]
tspan = (t_start, t_end) # NOTE: passed into ODEProblem

# Define initial conditions
x0 = 2  # [m]
v0 = 3  # [m/s]
u0 = [x0, v0] # NOTE: passed into ODEProblem; order persists in the ODEFunction and ODESolution

# Define the function that stores our system of equations.
# - It must be constructed as an ODEFunction Type
# - It will be called within every iteration of the numerical method performed by solve()
function mass_spring_damper!(du, u, p, t)

    # Rename states for clarity
    x = u[1]  # [m]
    v = u[2]  # [m/s]

    # Rename parameters for clarity
    m = p[1]  # [kg]
    k = p[2]  # [kg-m/s]
    b = p[3]  # [kg-m]
    F = p[4]  # [kg-m/s^2]

    # Define the differential equations
    #  u[1] = x, so du[1] = dx/dt
    #  u[2] = v, so du[2] = dv/dt
    du[1] = v                                   
    du[2] = -((k/m)*x) - ((b/m)*v) + ((1/m)*F)  

end

# Build the ODE problem
# Return Type: ODEProblem
prob = ODEProblem(mass_spring_damper!, u0, tspan, p)

# Build callback function, which will perturb the mass position (u[1])
function digital_controller_action(integrator)
    # Update a state variable with other states at this time step and/or stored parameters
    # NOTE: this expression is meaningless!!! It is just to illustrate an arbitrary perturbation of state u[1]
    integrator.u[1] += (integrator.p[3] * integrator.u[1])
end

# Create a periodic callback object to be passed into the solve function
# NOTE: there are many types of DiscreteCallbacks, this is defined as a periodic callback, 
#  which is run at some frequency Ts
# Return Type: DiscreteCallback
Ts = 10 # [s] callback frequency
cb = PeriodicCallback(digital_controller_action, Ts)

# Solve the ODE problem by specifying a numerical method
# Return Type: ODESolution
solution = solve(prob, RK4(), callback = cb)

# ---------------------------------------------------------------------------------
# Plot trajectory and integration step sizes
# ---------------------------------------------------------------------------------

# Create main plot 
p1 = plot(
    layout=(3,1),
    size=(600,800),
    plot_title="Mass Spring Damper w/ Periodic Forcing"
    )

# Add subplots 1 & 2: Time-domain trajectory of the state variables
plot!(solution, # NOTE: plot() automatically interprets/plots the ODESolution Type
    title=["Position of mass" "Velocity of mass"],
    label=["Position" "Velocity"], 
    xlabel="Time [s]", 
    ylabel=["Position [m]" "Velocity [m/s]"], 
    color=["blue" "red"]
    )

# Add subplot 3: Step sizes of the time integration (RK4)
# NOTE: This is primarily to illustrate how the time integrator "repeats" the 
#  previous time step when the perturbation is applied (i.e. time step looks like it  
#  is zero but really the integrator is performing a callback)
t_all = solution.t
n_steps = length(t_all)
step_size = t_all[2:n_steps] - t_all[1:n_steps-1]
plot!(
    subplot=3,
    t_all[2:n_steps], 
    step_size,
    title="Time steps taken by the time integrator",
    xlabel="Time [s]", 
    ylabel="\nStep Size [s]",
    label="Step Size",
    marker=:circle,
    markersize=3,
    yaxis=[-0.1,maximum(step_size*1.1)],
    legend=:bottomright,
)

# Add lines to all subplots to show the periodic perturbations
periodic_perturbation = collect(range(Ts,tspan[2],step=Ts))
vline!(
    [periodic_perturbation periodic_perturbation periodic_perturbation], 
    lab=["Periodic perturbation" "Periodic perturbation" "Periodic perturbation"], 
    l=(:black, :dash, 0.5)
    )
display(p1)