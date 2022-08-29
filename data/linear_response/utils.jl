using Revise
using LinearAlgebra
using DifferentialEquations
using FFTW

include("../plotting.jl")

unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))

function J(harm_eq::HarmonicEquation, vars, s)

	unbr = Dict(zip(get_variables(harm_eq), HarmonicBalance._remove_brackets(get_variables(harm_eq))))	
	this_J = HarmonicBalance.LinearResponse.get_Jacobian(harm_eq.equations, vars)
	this_J = HarmonicBalance.expand_derivatives.(HarmonicBalance.substitute_all(this_J, unbr))
	
	this_J = Matrix{ComplexF64}(HarmonicBalance.substitute_all(this_J, s))

end

J(res::Result, vars, s) = J(res.problem.eom, vars, s)

""" Extracting the signal amplitude and frequency conversion factor."""
function δU(res, s)

	eqs = res.problem.eom.equations
	vars = get_variables(res.problem.eom)
	T = get_independent_variables(res.problem.eom)[1]
	
	J0 = J(res, vars, s)
	J1 = J(res,(d(vars,T)), s)	

	omega = HarmonicBalance.substitute(ω, s)
	
	Jtil = inv(J1) * J0 - im*omega*I(2)
	
	eig = eigen(Jtil)
	R = eig.vectors
	D = Diagonal(eig.values)
	xi = [1, im] / 2;
	
    function r(Ω)       
    
        u,v = R * inv(im*Ω*I(2) + D) * inv(R) * inv(J1) * xi
        
        a = sqrt(2) * norm([u,v])
        w = (imag(u) * real(v) - real(u) * imag(v)) / norm([u,v])^2 + 1/2
        a , w
    end
	
end

function δU_corr(res, s, degree=3)

	eom = get_harmonic_equations(res.problem.eom.natural_equation, degree=degree)
	
	vars = get_variables(res.problem.eom)
	T = get_independent_variables(res.problem.eom)[1]
	
	J0 = J(eom, vars, s)
	J1 = J(eom,d(vars,T), s)	
	J2 = J(eom, d(vars,T,2), s)

	omega = HarmonicBalance.substitute(ω, s)
	
	xi = [1, im] / 2;
	
    function r(Ω)       
    
        u,v =  inv(-Ω^2*J2 + (2*J2*omega*I(2) + im*J1)*Ω - J2*omega^2 - im*J1*omega+J0) * xi
        
        a = sqrt(2) * norm([u,v])
        w = (imag(u) * real(v) - real(u) * imag(v)) / norm([u,v])^2 + 1/2
        a , w
    end


end
###
# TIME-DEPENDENT SIMULATIONS
###

function solution_lab(;ω, ξ, Ω, x0 = [0.3,0.3], α=1.0, g=0.01, ω0=1.0, γ=0.001, η=0.1, cycles=5000, pts_per_cycle=30)
    T = cycles * 2pi / ω
    f(du, u, p, t) = -ω0^2  * u - α*u^3 - γ*du - η*u*u*du + g*cos(ω*t) + ξ*cos(Ω*t)
    problem = SecondOrderODEProblem(f, x0..., [0., T])
    soln = DifferentialEquations.solve(problem, saveat=(1/pts_per_cycle) * 2pi/ω, reltol=1e-8, abstol=1e-8, maxiters=1e8)
    soln.t, getindex.(soln.u, 2)
end


function Fourier_filter(response; ω_range=[0.5,1.5], ringdown_cycles=10000, pts_per_cycle)
    Δt = response[1][2] - response[1][1]
    
    #k = Int(round(ringdown_time / Δt) + 1)
    k = pts_per_cycle * ringdown_cycles + 1
    
    fft = FFTW.fft(response[2][k:end]) .* Δt # the spectrum
    
    f_axis = FFTW.fftfreq(length(fft), 1/Δt)
    
    f_range = ω_range ./ (2pi)
    
    
    relevant_indices = findall(x -> f_range[1] <= abs(x) <= f_range[2], f_axis)
    fft_filtered = fft[relevant_indices]
    
    parseval_sum = sum((abs.(fft_filtered).^2)) / (length(fft)*Δt)
    sqrt(2 * parseval_sum / (response[1][end] - response[1][k]))
    
end


###
# OLD 
###

function xi_bar(u,v, ω,γ,η)
	den = (4*γ + (u^2 + v^2)*η) * (4 * γ + 3* (u^2 + v^2) * η) + 64 * ω^2
	xi_u = 2  *(4* γ + (u + im* v)* (u - 3* im* v)* η - 8* im* ω)
	xi_v = 8 *im* γ + 2* (u + im* v) *(3 *im* u + v) * η + 16* ω
	[xi_u, xi_v] ./ den
	
end

function δU_old(res, s)
    J = res.jacobian(s)
    omega = HarmonicBalance.substitute(ω, s) 
    eig = eigen(J + im*omega*I(2))
    
    R = eig.vectors
    D = Diagonal(eig.values)
    
    sub(a) = HarmonicBalance.substitute(a,s)
    omega = sub(ω)
    xi = xi_bar( s.vals[1], s.vals[2], sub(ω), sub(γ), sub(η))
    
    function r(Ω)       
    
        u,v = R * inv(im*Ω*I(2) - D) * inv(R) * xi
        
        a = norm([u,v])
        w = (imag(u) * real(v) - real(u) * imag(v)) / a^2 + 1/2
        sqrt(2) * a , w
    end
end

