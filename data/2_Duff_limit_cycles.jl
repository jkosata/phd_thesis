using HarmonicBalance

@variables γ, F, α, ω0, F0, η, J, ω, t, Δω, x(t), y(t);

# a vector of expressions - these must equal to zero
diff_eq = DifferentialEquation([d(x,t,2) + γ * d(x,t) + ω0^2 * x + α*x^3+ 2*J*ω0*(x-y) - F0*cos(ω*t), 
            d(y,t,2) + γ * d(y,t) + ω0^2 * y + α*y^3 + 2*J*ω0*(y-x) - η*F0*cos(ω*t)], [x,y])

# describe each variable using one or more Fourier components
add_harmonic!(diff_eq, x, ω)
add_harmonic!(diff_eq, y, ω)

# add a Hopf frequency Δω
add_Hopf!(diff_eq, Δω)
harmonic_eq = get_harmonic_equations(diff_eq);

# remove one Hopf variable, replace by Δω, construct the appropriate Jacobian
prob = HarmonicBalance.Hopf._Hopf_Problem(harmonic_eq, Δω);

fixed_parameters = (
    ω0 => 1.4504859, # natural frequency of separate modes (in paper, ħω0 - J)
    γ => 27.4E-6,    # damping
    J => 154.1E-6,   # coupling term
    α => 3.867E-7,   # Kerr nonlinearity
    ω => 1.4507941,  # pump frequency, resonant with antisymmetric mode (in paper, ħω0 + J)
    η => -0.08,      # pumping leaking to site 2  (F2 = ηF1)
    F0 => 0.011       # pump amplitude (overriden in sweeps)
)

range = F0 => LinRange(0.002, 0.05, 500)
solutions = get_steady_states(prob, range, fixed_parameters, random_warmup=true)
HarmonicBalance.save("2_Duff_limit_cycles.jld2", solutions)

