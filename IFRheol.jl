# using DifferentialEquations
using OrdinaryDiffEq
using Plots

# create a triangular distribution of filament slack w(epsilon_s)
# constraint of normalisation: int w = 1


function get_slack_dist(emin,emax,de=0.001,elimit=3)
  e=Array(0:de:elimit)
  w = e * 0  # copy and reset array

  # set triangular distribution
  emid=(emin+emax)/2

  w[e.>emin] = e[e.>emin] .-emin
  w[e.>emid] = emax .- e[e.>emid]
  w[e.>emax] .= 0

  # normalise distribution to 1
  w /= sum(w)

  return( (e=e[w.!=0], w=w[w.!=0]) )
end

# get the local slope of the non-linear stress-strain curve
# this is just counting the number of engaged springs for a certain deformation e
function local_stiffness(e,ec,k)
  return( k * sum(ec.w[ec.e .<= e]) )
end

# format for functions 
# return strain and its derivative as a function of time for the differential equation solver.
#
# This function produces a simple ramp starting at t=0
function strain_ramp(rate)
  return( t-> (e=rate*t, edot=rate) )
end

# This function produces a cycles of a certain strain amplitude and rate starting at t=0
function strain_cycle(ampl,rate)
  return( t -> begin; tm = t%(2*ampl/rate); if tm <= ampl/rate; (e=rate*tm, edot=rate); else (e=-rate*tm+2*ampl, edot=-rate) end; end;)
end


# function for single non-linear branch

function f_sb(v,p,t)
  e, edot = p.strain_function(t)
  s=v[1]
  ed=v[2]

  dsdt= local_stiffness(e-ed,p.ec,p.k) * ( edot - s / p.eta ) 
  deddt = s / p.eta
  return([dsdt,deddt]) 
end


# function f_sb(s,p,t)
#   e, edot = p.strain_function(t)

#   dsdt= local_stiffness(e,p.ec,p.k) * ( edot - s / p.eta ) 
#   return(dsdt) 
# end



function solve_sb(loading, duration, eta=30, dt=0.1)
  p = (k = 1.0, eta = eta, strain_function = loading, ec = get_slack_dist(0.5, 1.5))
  tspan = (0.0,duration)
  prob = ODEProblem(f_sb,[0,0],tspan,p,reltol=1e-8, abstol=1e-8, saveat=dt)
  sol = solve(prob,Tsit5())
  # calculate gradient
  e=[loading(t).e for t in sol.t]
  s=[e[1] for e in sol.u]
  g=[0;s[2:end]-s[1:end-1]] ./ [1;e[2:end]-e[1:end-1]] 
  return( (t=sol.t, e=e, s=s, g=g, ed=[e[2] for e in sol.u]) )
end

sol=solve_sb(0.0025,4)
plot(sol.e, sol.g)
sol=solve_sb(0.01,4)
plot!(sol.e, sol.g)
sol=solve_sb(0.05,4)
plot!(sol.e, sol.g)



# function for multiple non-linear branches

function f_mb(v,p,t)
  e, edot = p.strain_function(t)
  s=v[1]
  ed=v[2]
  # only valid if s is positive? would this work for compression?
  dsdt= ( if (e-ed) <= p.ec; 0; else p.k end; ) * ( edot - s / p.eta ) 
  deddt = s / p.eta
  return( [dsdt, deddt] ) 
end


# function for multiple non-linear branches

function f_mb(s,p,t)
  e, edot = p.strain_function(t)

  # only valid if s is positive? would this work for compression?
  dsdt= ( if e <= p.ec; 0; else p.k end; ) * ( edot - s / p.eta ) 
  return(dsdt) 
end


function solve_mb(loading, duration, eta=30, dt=0.1)
  tspan = (0.0,amplitude/rate)
  ta = Array(0:dt:amplitude/rate)
  s = zeros( length(ta) )
  
  e=[strain_ramp(t,rate).e for t in ta]
  ec = get_slack_dist(0.5, 1.5)
  
  for i in 1:length(ec.w)
    p = (k = 1.0, eta = eta, strain_function = loading, ec = ec.e[i])
    prob = ODEProblem(f_mb,0,tspan,p,reltol=1e-8, abstol=1e-8,saveat=dt)
    sol = solve(prob,Tsit5())
    s += ec.w[i]*sol.u
  end
  # calculate gradient
  g=[0;s[2:end]-s[1:end-1]] ./ [1;e[2:end]-e[1:end-1]] 
  return( (t=ta, e=e, s=s, g=g) )
end


function plot_mb(plt)
  for rate in [0.001, 0.002, 0.005, 0.01, 0.02]
    sol = solve_mb(rate,3)
    plot!(plt, sol.e, sol.g)
  end
end

sol = solve_mb(0.001,3)
plot(sol.e, sol.g)
sol = solve_mb(0.002,3)
plot!(sol.e, sol.g)
sol = solve_mb(0.005,3)
plot!(sol.e, sol.g)
sol = solve_mb(0.01,3)
plot!(sol.e, sol.g)
sol = solve_mb(0.02,3)
plot!(sol.e, sol.g)
sol = solve_mb(0.05,3)
plot!(sol.e, sol.g)
sol = solve_mb(0.1,3)
plot!(sol.e, sol.g)
       


sol = solve_sb(0.001,3)
plot(sol.e, sol.g)
sol = solve_sb(0.002,3)
plot!(sol.e, sol.g)
sol = solve_sb(0.005,3)
plot!(sol.e, sol.g)
sol = solve_sb(0.01,3)
plot!(sol.e, sol.g)
sol = solve_sb(0.02,3)
plot!(sol.e, sol.g)
sol = solve_sb(0.05,3)
plot!(sol.e, sol.g)
sol = solve_sb(0.1,3)
plot!(sol.e, sol.g)
       


using Plots
pyplot()

# sigma(t) for given tau, e0 distribution, strain rate,scale


w[e.>0.5] = e[e.>0.5] .-0.5
w[e.>1.0] = 1.5 .- e[e.>1.0]
w[e.>1.5] = 0


#
#  This function return the modulus (gradient of strain response) for a distribution of taut deformation that is flat between e_0 and e_1.
#  for a given strain rate, relaxation time of the maxwell unit
#

function keff_square(eps,edot,tau,e0,e1)
  a = 1/(e1-e0)
  ke = copy(eps)
  for i in 1:length(eps)
    e = eps[i]
    if e <= e0
      ke[i]=0
    elseif e <= e1
      ke[i] = a * edot * ( 1 - exp(-(e-e0)/(edot*tau)) )
    else
      ke[i] = a * edot * ( 1 - (exp(-(e-e0)/(edot*tau)) ) - ( 1 - exp(-(e-e1)/(edot*tau)) ) )
    end
  end
  return(ke)
end

function keff_triangular(eps,edot,tau,e0,e1)
  mid=(e0+e1)/2
  h = 2 / (e1-e0)
  #a = h/(mid-e0)    #   To be completed - scale so that the integral of the distribution is 1
  ke = (h*de) .* accumulate(+,keff_square(eps,edot,tau,e0,mid) - keff_square(eps,edot,tau,mid,e1))
  return(ke)
end


function rheol_triangular(eps,edot,tau,e0,e1,k0)
  ke = keff_triangular(eps,edot,tau,e0,e1) .+ k0
  ta = eps./edot
  de = eps[2]-eps[1]  # assumes uniform sampling here
  sigma = accumulate(+,ke.*de)
  return( (t=ta, e=eps, s=sigma, grad=ke) )
end


# r=rheol_triangular(e,0.01,30,0.2,1.2,0.01)

