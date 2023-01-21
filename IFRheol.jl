# using DifferentialEquations
using OrdinaryDiffEq
using Plots

pyplot()

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


function solve_mb(loading, duration, eta=30, dt=0.1)
  tspan = (0.0,duration)
  ta = Array(0:dt:duration)
  s = zeros( length(ta) )
  
  e=[loading(t).e for t in ta]
  ec = get_slack_dist(0.5, 1.5)
  
  for i in 1:length(ec.w)
    p = (k = 1.0, eta = eta, strain_function = loading, ec = ec.e[i])
    prob = ODEProblem(f_mb,[0,0],tspan,p,reltol=1e-8, abstol=1e-8,saveat=dt)
    sol = solve(prob,Tsit5())
    sb=[r[1] for r in sol.u]
    s += ec.w[i]*sb
  end
  # calculate gradient
  g=[0;s[2:end]-s[1:end-1]] ./ [1;e[2:end]-e[1:end-1]] 
  return( (t=ta, e=e, s=s, g=g) )
end




# A few lines to test the functions and plot the output


rate = 0.001
sol = solve_sb(strain_cycle(2,rate) ,4/rate   )
plot(sol.e, sol.s)
 rate = 0.002
sol = solve_sb(strain_cycle(2,rate) ,4/rate   )
plot!(sol.e, sol.s)
 rate = 0.005
sol = solve_sb(strain_cycle(2,rate) ,4/rate   )
plot!(sol.e, sol.s)
 rate = 0.01
sol = solve_sb(strain_cycle(2,rate) ,4/rate   )
plot!(sol.e, sol.s)
 rate = 0.02
sol = solve_sb(strain_cycle(2,rate) ,4/rate   )
plot!(sol.e, sol.s)
 



rate = 0.001
sol = solve_mb(strain_cycle(2,rate) ,4/rate   )
plot(sol.e, sol.s)
 rate = 0.002
sol = solve_mb(strain_cycle(2,rate) ,4/rate   )
plot!(sol.e, sol.s)
 rate = 0.005
sol = solve_mb(strain_cycle(2,rate) ,4/rate   )
plot!(sol.e, sol.s)
 rate = 0.01
sol = solve_mb(strain_cycle(2,rate) ,4/rate   )
plot!(sol.e, sol.s)
 rate = 0.02
sol = solve_mb(strain_cycle(2,rate) ,4/rate   )
plot!(sol.e, sol.s)
 


