using OrdinaryDiffEq
using Interpolations
using Calculus
using Plots

pyplot()

# 1. STRAIN RAMP ANALYTICAL SOLUTION 
function strain_ramp(rate)
    """
        # Aim:
        - Create a function that represents a linear strain ramp starting from time t=0.

        # Arguments
        - `rate`: The rate of strain increase per unit time.

        # Returns
        - A function of time `t` that, when called, returns a tuple with: `e` and `edot`
    """
    return( t-> (e=rate*t, edot=rate) )
end

function get_slack_dist(emin,emax,de=0.05,elimit=3)
    """
        # Aim:
        - Calculate and return a normalized triangular distribution for strain between specified limits.

        # Arguments
        - `emin`: The minimum value of the distribution range.
        - `emax`: The maximum value of the distribution range.
        - `de`: The increment between successive distribution values, with a default of 0.05.
        - `elimit`: The upper limit for the distribution values, with a default of 3.

        # Returns
        - A tuple containing two elements:
        - `e`: An array of distribution values within the specified range.
        - `w`: Corresponding weights of the distribution, normalized so their sum equals 1.

    """
    e=Array(0:de:elimit)
    w = e * 0
  
    # set triangular distribution
    emid=(emin+emax)/2
  
    w[e.>emin] = e[e.>emin] .-emin
    w[e.>emid] = emax .- e[e.>emid]
    w[e.>emax] .= 0
  
    # normalise distribution to 1
    w /= sum(w)
  
    return( (e=e[w.!=0], w=w[w.!=0]) )
end

function local_stiffness(e,ec,k)
    """
        # Aim:
        - Calculate the local stiffness for a given element in a computational model.

        # Arguments:
        - `e`: The strain value for which the local stiffness is being calculated.
        - `ec`: A struct or composite type that contains the weights (`w`) and elements (`e`) of the entire model.
        - `k`: The spring constant initialisation value.

        # Returns
        - Returns the local stiffness by counting the number of engaged springs for a certain deformation e
    """
    return( k * sum(ec.w[ec.e .<= e]) )
end

  
function f_sb(v,p,t)
    """
        # Aim:
        - Calculate the time derivatives of stress and dashpot strain in a single branch (non-linear spring) model.

        # Arguments
        - `v`: A vector where `v[1]` is the current stress, and `v[2]` is the current dashpot strain.
        - `p`: A parameter struct containing: `strain_function`, `ec`, `k`, `ks`, `eta`
        - `t`: The current time.
    
        # Returns
        - A vector containing:
            - The time derivative of stress (`dsdt`).
            - The time derivative of dashpot strain (`deddt`).
    """
    e, edot = p.strain_function(t)
    s=v[1]
    ed=v[2]

    if local_stiffness(e-ed,p.ec,p.k) == 0
        # if single branch is not recuited, only simulate the spring in parallel
    return ([edot*p.ks, 0.0])
    
    else 
    # To determine local stiffness they do e-ed which removes the plastic deformation
        s1 = s - e*p.ks
        dsdt= local_stiffness(e-ed,p.ec,p.k) * ( edot - (s1) / p.eta ) + edot*p.ks
        deddt = (s1) / p.eta
    end
    return([dsdt,deddt]) 
end

function solve_sb(loading, duration, eta=8.25e5, dt=0.01)
    """
        # Aim:
        - Solve an ODE problem for the single nramch model over a given duration.

        # Arguments
        - `loading`: A function representing the loading conditions
        - `duration`: The total time duration
        - `eta`: Viscosity parameter with a default value of `8.25e5`. This can be adjusted based on specific requirements of the simulation.
        - `dt`: Time step size for the solver, with a default value of `0.01`.
        
        # Returns
        - A tuple containing:
          - `t`: An array of time points at which the solution is computed.
          - `e`: An array of strain values corresponding to each time point.
          - `s`: An array of stress values corresponding to each time point.
          - `ed`: An array of dashpot strain values corresponding to each time point.
        
    """
    p = (k = 15e3, ks=1e3, eta = eta, strain_function = loading, ec = get_slack_dist(0.0, 2.2))
    tspan = (0.0,duration)
    prob = ODEProblem(f_sb,[0,0],tspan,p,reltol=1e-8, abstol=1e-8, saveat=dt)
    sol = solve(prob,Vern9(), dtmax=dt, maxiters=10000000)

    # Unpack solution
    t = sol.t
    e=[loading(t).e for t in sol.t]
    s=[i[1] for i in sol.u]
    ed = [i[2] for i in sol.u]
    return( (t=t, e=e, s=s, ed=ed))
end

########################################################## Plotting ###################################################
rate = 0.01
duration = 500
sol = solve_sb(strain_ramp(rate), duration)
# plot
# p1 = plot(sol.t, sol.e, xlabel="strain", ylabel="time", title="Strain ramp", legend=false)
# p2 = plot(sol.t, sol.s, xlabel="stress", ylabel="time", title="Stress response", legend=false)
# p3 = plot(sol.e, sol.s, xlabel="strain", ylabel="stress", title="Single branch", legend=false)
# plot(p1, p2, p3, layout=(1,3))
#######################################################################################################################

# 2. STRESS LOADING FEEDBACK LOOP 
function create_stress_loading_function(sol)
    """
        # Aim:
        - Create a function for interpolating stress and its time derivative based on the strain ramp 
        loading solution.

        # Arguments
        - `sol`: A solution object obtained from solving the ODE for strain ramp loading: 
                contains time points `t` and corresponding stress values `s`.
        
        # Returns
        - A function that, when given a time `t`, returns a tuple containing:
          - `s`: The interpolated stress value at time `t`.
          - `ds_dt`: The derivative of stress with respect to time at `t`.      
    """
    # Extract time and stress from sol
    ts = sol.t
    ss = sol.s 

    # Create an interpolation object for stress
    stress_interp = LinearInterpolation(ts, ss, extrapolation_bc=Flat())

    # Create an interpolation object for the derivative of stress
    # Using the derivative() function from the Interpolations package
    stress_derivative_interp = derivative(stress_interp)

    # Return a function that interpolates stress and its derivative based on time
    return t -> (s = stress_interp(t), ds_dt = stress_derivative_interp(t))
end

function f_sb_sp(v,p,t)
    """
        # Aim: 
        - Calculate the time derivatives of total strain and dashpot strain for a single branch (non-linear spring) model.

        # Arguments
        - `v`: A vector where v[1] is the current total strain, and v[2] is the current dashpot strain.
        - `p`: An object that includes the model parameters e.g. stress function,'k', and 'eta'.
        - `t`: The current time.

        # Returns
        - A vector where the first element is the time derivative of the total deformation (dedt),
            and the second element is the time derivative of the dashpot deformation (deddt).
    """
    s, sdot = p.stress_function(t)
    e=v[1]
    ed=v[2]
    
    if local_stiffness(e-ed,p.ec,p.k) == 0
        # if single branch is not recuited, only simulate the spring in parallel
        return ([sdot/p.ks, 0.0])
    else 
        # find stress in the non-linear brach with dashpot and non-linear spring (= single branch)
        s1 = s - e/p.ks

        # find change in strain of Dashpot
        deddt = s1 / p.eta

        # find total change in strain using analytical ODE solution
        """
            Δε = Δεd + Δεs                              (1) (strains of series components in the non-linear branch add up)
            Δε = σ1/η + Δσ1/k1      

            σ = σ1 + σ2                                 (2) (stresses of parallel branches add up)

            σ1 = σ - σ2                                 (3) (stress of non-linear branch can be written in terms of σ and σ2)

            Δε = (σ - σ2)/η + Δ(σ - σ2)/k1              (4) (re-write the total change in strain in terms of σ and σ2)

            σ2 = ε2 * k2                                (5) (stress in the sigle spring branch can be written in terms of ε2 and k2)

            ε = ε1 = ε2                                 (6) (strains in parallel branches are equal)

            Δε = (σ - ε*k2)/η + Δ(σ - ε*k2)/k1          (7) (re-write the total change in strain in terms of σ and σ2)

            Δε = σ/η - ε*k2/η + 1/k1 (Δσ/Δt - Δε/Δt*k2) (8) (expand) 

            Δε (1 + k2/k1) = 1/k1 * Δσ + σ/η - ε*k2/η   (9) (simplify, using Δσ = 0 at plateau and collecting like terms) 

        """
        dedt = ((1/local_stiffness(e-ed,p.ec,p.k))*sdot + s/p.eta - e*p.ks/p.eta) / (1 + p.ks / local_stiffness(e-ed,p.ec,p.k))
        return([dedt, deddt])
    end
end

function solve_sb_sp(loading, duration, eta=8.25e5, dt=0.01)
    """
        # Aim:
        - Solve a differential equation problem for a specified loading function, 
        over a given duration, loading and tolerances.
        
        # Arguments
        - `loading`: A function representing the loading conditions.
        - `duration`: The total time duration for which the problem is solved.
        - `eta`: Viscosity parameter with a default value of `8.25e5`, adjusted based on the material
        - `dt`: Time step size for the solver, with a default value of `0.05`.
    
        # Returns
        - A tuple containing time (`t`), stress (`s`), strain (`e`), and strain rate (`ed`) values as arrays.
    """
    # Initialize emin and emax which are inputs for the slack distribution
    emin = 0.0
    emax = 2.2
    
    # Create a parameter tuple 'p' including: spring stiffnesses, eta, loading, and slack distribution
    p = (k = 15e3, ks=1e3 ,eta = eta, stress_function = loading, ec = get_slack_dist(emin, emax)) 
    tspan = (0.0,duration)

    # ODE problem is defined over the time span from 0 to `duration`, with initial values `[e, ed] = [0, 0]`
    prob = ODEProblem(f_sb_sp,[0,0],tspan,p,reltol=1e-8, abstol=1e-8, saveat=dt) 
    
    # The solver used is `Vern9`, with a maximum iteration limit of 20,000,000 and a maximum time step of `dt`
    sol2 = solve(prob,Vern9(), dtmax=dt, maxiters=20000000) #solve(prob,Tsit5())
    
    # Unpack sol
    t = sol2.t
    s = [loading(t).s for t in sol2.t]
    e = [i[1] for i in sol2.u]
    ed = [i[2] for i in sol2.u]
    return((t=t, s=s, e=e, ed=ed))
end

############################################################# Plotting ###################################################
duration = 500
sol2 = solve_sb_sp(create_stress_loading_function(sol), duration)
# # plot
# p1 = plot(sol.t, sol.e, xlabel="time", ylabel="strain", title="Strain ramp", legend=false)
# p1 = plot!(sol2.t, sol2.e, xlabel="time", ylabel="strain", title="Strain ramp", legend=false)
# p2 = plot(sol.t, sol.s, xlabel="time", ylabel="stress", title="Stress response", legend=false)
# p2 = plot!(sol2.t, sol2.s, xlabel="time", ylabel="stress", title="Stress response", legend=false)
# p3 = plot(sol.e, sol.s, xlabel="strain", ylabel="stress", title="Single branch", legend=false)
# p3 = plot!(sol2.e, sol2.s, xlabel="strain", ylabel="stress", title="Single branch", legend=false)
# plot(p1, p2, p3, layout=(1,3))
##########################################################################################################################

# 3. STRAIN RAMP PID CONTROL

function pid_stain_loading(gain_Kp, gain_Kd, sol)
    """
        # Aim:
        - Creates a function for PID control based on strain loading.

        # Arguments
        - `gain_Kp`: Proportional gain for the PID controller.
        - `gain_Kd`: Derivative gain for the PID controller.
        - `sol`: A solution object typically obtained from solving a differential equation. 
         It should contain time points `t` and corresponding stress values `s`.

        # Returns
        - A function that,for time `t`, returns the rate of change of the dashpot strain (`eddot`), calculated using PID control logic.
    """
    ts = sol.t
    ss = sol.s 

    # Create an interpolation object for stress
    stress_interp = LinearInterpolation(ts, ss, extrapolation_bc=Flat())
    return((t, s, dsdt)-> (eddot= gain_Kp * (stress_interp(t)-s) - gain_Kd * dsdt))
end

function f_sb_pid(v,p,t)
    """
        # Aim:
            - Calculate the time derivatives of stress, dashpot strain, total strain, and strain rate in a 
            single branch (non-linear spring) model with PID control.

        # Arguments
            - `v`: A vector where `v[1]` is the current stress, `v[2]` is the current dashpot strain, `v[3]` is the total strain, and `v[4]` is the strain rate.
            - `p`: A parameter struct containing: `strain_function`, `ec`, `k`, `ks`, `eta`
            - `t`: The current time.

         # Returns
            - A vector containing:
                - The time derivative of stress (`dsdt`).
                - The time derivative of dashpot strain (`deddt`).
                - The strain rate (`edot`).
                - The acceleration of the dashpot strain (`eddot`).
    """
    s=v[1]
    ed=v[2]
    e = v[3]
    edot = v[4]

    if local_stiffness(e-ed,p.ec,p.k) == 0
        # if single branch is not recuited, only simulate the spring in parallel
        dsdt = edot*p.ks
        # call pid to find change in edot to change the stress correctly
        eddot = p.strain_function(t, s, dsdt)
        
        return ([dsdt, 0.0, edot, eddot])
    else 
        # To determine local stiffness they do e-ed which removes the plastic deformation
        s1 = s - e*p.ks
        dsdt= local_stiffness(e-ed,p.ec,p.k) * ( edot - (s1) / p.eta ) + edot*p.ks
        deddt = (s1) / p.eta
        # call pid to find change in edot to change the stress correctly
        eddot = p.strain_function(t, s, dsdt)
        return([dsdt, deddt, edot, eddot])
    end
end

function solve_sb_pid(loading, duration, eta=8.25e5, dt=0.01)
    """
        # Aim:
        - Solve an ODE problem under PID control.

        # Arguments
        - `loading`: A PID control function that takes `t`, `s`, and `dsdt`, and returns `eddot`.
        - `duration`: The total time for which the ODE problem is solved.
        - `eta`: Viscosity parameter with a default value of `8.25e5`, can be adjusted based on the material.
        - `dt`: Time step size for the solver, with a default value of `0.01`.
        
        # Returns
        - A tuple containing:
          - `t`: An array of time points at which the solution is computed.
          - `e`: An array of total strain values corresponding to each time point.
          - `s`: An array of stress values corresponding to each time point.
          - `ed`: An array of dashpot strain values corresponding to each time point.
    """
    edot_0 = 0.01
    p = (k = 15e3, ks=1e3, eta = eta, strain_function = loading, ec = get_slack_dist(0.0, 2.2))
    tspan = (0.0,duration)
    prob = ODEProblem(f_sb_pid,[0, 0, 0, edot_0], tspan,p,reltol=1e-8, abstol=1e-8, saveat=dt)
    sol = solve(prob, Vern9(), dtmax=dt, maxiters=10000000)
    # Unpack solution
    t = sol.t
    s = [i[1] for i in sol.u]
    ed = [i[2] for i in sol.u]
    e = [i[3] for i in sol.u]
    return( (t=t, e=e, s=s, ed=ed))
end

#################################################################### Plotting #######################################################################
gain_Kp = 1.0 # (s)
gain_Kd = 0.5 # (s)
duration = 500
sol3 = solve_sb_pid(pid_stain_loading(gain_Kp, gain_Kd, sol), duration)

p1 = plot(sol.t, sol.e, xlabel="time", ylabel="strain", title="Strain ramp", label=false)
p1 = plot!(sol2.t, sol2.e, xlabel="time", ylabel="strain", title="Strain ramp")
p1 = plot!(sol3.t, sol3.e, xlabel="time", ylabel="strain", title="Strain ramp")
p2 = plot(sol.t, sol.s, xlabel="time", ylabel="stress", title="Stress response")
p2 = plot!(sol2.t, sol2.s, xlabel="time", ylabel="stress", title="Stress response")
p2 = plot!(sol3.t, sol3.s, xlabel="time", ylabel="stress", title="Stress response")
p3 = plot(sol.e, sol.s, xlabel="strain", ylabel="stress", title="Single branch")
p3 = plot!(sol2.e, sol2.s, xlabel="strain", ylabel="stress", title="Single branch")
p3 = plot!(sol3.e, sol3.s, xlabel="strain", ylabel="stress", title="Single branch")
plot(p1, p2, p3, layout=(1,3), label=["Strain loading" "Stress loading" "PID control"], legendfontsize=5, legend=:bottomright)
#####################################################################################################################################################