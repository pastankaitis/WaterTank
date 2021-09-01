using ReachabilityAnalysis, Plots, Symbolics, MathematicalSystems, SparseArrays


#constants

const H_max = 3 # Maximum height of the tank 
const B = 3 # Base of the cylinder tank
const V_max = 3 # Maximum volume flow for the tank
const V_saf_max = 2.9 # minimum allowed water level
const V_saf_min = 0.5 # minimum allowed water level

# water tank modes

@taylorize function stable!(du, u, p, t)

    d_in, d_out, h = u

    #differential equations for the stable water tank mode
    du[1] = zero(u[1])
    du[2] = zero(u[1])
    du[3] = zero(u[1])

    return du


end

@taylorize function filling!(du, u, p, t)


    d_in, d_out, h = u
    
    #differential equations for the filling water tank mode
    du[1] = one(u[1])
    du[2] = zero(u[1])
    du[3] = d_in


    return du

end

@taylorize function emptying!(du, u, p, t)

    d_in, d_out, h = u

    #differential equations for the emptying water tank mode
    du[1] = zero(u[1])
    du[2] = one(u[1])
    du[3] = -d_out

    return du


end

@taylorize function normal!(du, u, p, t)

    d_in, d_out, h = u

    #differential equations for the normal water tank mode
    du[1] = one(u[1])
    du[2] = one(u[1])
    du[3] = d_in - d_out


    return du

end

function waterTank()

    var = @variables d_in, d_out, h

    automaton = LightAutomaton(3)
    
    #X = HPolyhedron([d_in > 0, d_out > 0, B * d_in <= V_max, B * d_out <= V_max], var) # invariants
    #stable_mode = @system(h' = stable!(h), dim:3, h ∈ X) # water tank mode stable

    X = HPolyhedron([d_in > 0, d_out > 0, B * d_in <= V_max, B * d_out <= V_max], var) # invariants
    normal_mode = @system(h' = normal!(h), dim:3, h ∈ X) # water tank normal mode

    X = HPolyhedron([d_in > 0, d_out > 0, B * d_in <= V_max, B * d_out <= V_max], var) # invariants
    filling_mode = @system(h' = filling!(h), dim:3, h ∈ X) # water tank filling mode

    X = HPolyhedron([d_in > 0, d_out > 0, B * d_in <= V_max, B * d_out <= V_max], var) # invariants
    emptying_mode = @system(h' = emptying!(h), dim:3, h ∈ X) # water tank emptying mode


    # transition "stable" -> "normal"
    #add_transition!(automaton, 1, 2, 1)
    #guard_1 = HalfSpace(h >= 0, var)
    #t1 = @ConstrainedIdentityMap(h -> h, dim: 3)

   
    # transition "normal" -> "filling"
    add_transition!(automaton, 1, 2, 1)
    guard_2 = HalfSpace(h <= V_saf_min, var)
    t1 = ConstrainedIdentityMap(3, guard_2) 

    # transition "normal" -> "emptying"
    add_transition!(automaton, 1, 3, 2)
    guard_3 = HalfSpace(h >= V_saf_max, var)
    t2 = ConstrainedIdentityMap(3, guard_3)

    # transition "filling" -> "normal"
    add_transition!(automaton, 2, 1, 3)
    guard_4 = HPolyhedron([h >= V_saf_min, h <= V_saf_max], var)
    t3 = ConstrainedIdentityMap(3, guard_4) 

    # transition "emptying" -> "normal"
    add_transition!(automaton, 3, 1, 4)
    guard_5 = HPolyhedron([h >= V_saf_min, h <= V_saf_max], var)
    t4 = ConstrainedIdentityMap(3, guard_5) 


    H = HybridSystem(automaton=automaton,
                     modes=[normal_mode, filling_mode, emptying_mode],
                     resetmaps=[t1, t2, t3, t4])


    ## initial condition in mode 1
    X0 = Hyperrectangle([0.1, 0.1, 1.5], [0.5, 0.5, 2.0])
    init = [(1, X0)]

    return InitialValueProblem(H, init)


end

prob = waterTank()

boxdirs = BoxDirections(3)


sol = solve(prob,
                tspan=(0.0, 200.0),
                alg=TMJets20(abstol=1e-5, maxsteps=10000, orderT=3, orderQ=1, disjointness=BoxEnclosure()),
                intersect_source_invariant=false,
                intersection_method=TemplateHullIntersection(boxdirs),
                clustering_method=LazyClustering(1),
                disjointness_method=BoxEnclosure())


plot(sol, vars=(0, 3), xlab="t", ylab="h", lc=:blue)






