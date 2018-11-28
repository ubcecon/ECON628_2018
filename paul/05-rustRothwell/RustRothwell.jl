"""
Module for reproducing Rust & Rothwell (1995)
"""
module RustRothwell

export panel_lag, load_rust_rothwell_data,
  probfn, payoffs, transition_prob, feasible_actions
  
using DataFrames, DelimitedFiles, ShiftedArrays, StatsModels, GLM,
  SparseArrays

"""
    read_plant_file(filename)

Reads data file from Rust & Rothwell 1995 about operation of a nuclear
power plant.

Input:
  - `filename` name of file

Output:
  - `DataFrame`

"""
function read_plant_file(filename)
  tmp = DataFrame(readdlm(filename),
                 [:name,:ID,:steamSystem,:year,:month,
                  :vintage,:age,:hrsRefuel,:hrsPlanOut,
                  :hrsForcedOut,:hrsTotal,:scramIn,:scramOut,
                  :numForcedOut])
  # put correct types on each variable
  # FIXME: is there a better way to do this?
  pd = DataFrame([String, Int64, Int64, Int64,Int64,
                  Int64,Int64, Float64, Float64,
                  Float64, Float64, Int64, Int64,
                  Int64], names(tmp), nrow(tmp))
  for v in names(pd)
    pd[v] .= tmp[v]
  end
  return(pd)
end

"""
    panel_lag(x::Symbol, data::AbstractDataFrame, id::Symbol, t::Symbol,
              lags::Integer=1)

Create lags of variables in panel data.

Inputs:
  - `x` variable to create lag of
  - `data` DataFrame containing `x`, `id`, and `t`
  - `id` cross-section identifier
  - `t` time variable
  - `lags` number of lags. Can be negative, in which cause leads will
     be created

Output:
  - A vector containing lags of data[x]. Will be missing for `id` and
   `t` combinations where the lag is not contained in `data`.

"""
function panel_lag(x::Symbol, data::AbstractDataFrame, id::Symbol, t::Symbol,
                   lags::Integer=1)
  if (!issorted(data, (id, t)))
    @warn "data is not sorted, panel_lag() will be more efficient with a sorted DataFrame"
    p = sortperm(data, (id, t))
    df = data[p,:]
  else
    p = nothing
    df = data
  end
  idlag= lag(df[id], lags)
  tlag = lag(df[t], lags)
  xlag = lag(df[x], lags)
  xlag = copy(xlag)
  xlag[ (typeof.(tlag).==Missing) .|
        (tlag .!= df[t].-lags) .|
        (idlag .!= df[id]) ] .= missing
  if (p == nothing)
    return(xlag)
  else
    xout = similar(xlag)
    xout[p] .= xlag
    return(xout)
  end
end

"""
     load_rust_rothwell_data()

Downloads, loads, and prepares Rust & Rothwell (1995) data.

"""
function load_rust_rothwell_data()
  # download data if necessary
  if !isfile("rr-data.zip")
    download("http://qed.econ.queensu.ca/jae/1995-v10.S/rust-rothwell/rr-data.zip",
             "rr-data.zip")
  end
  # unzip data if necessary
  if !isfile("rr-data/yank1.mon")
    run(`unzip rr-data.zip -d rr-data`)
  end

  files = filter(x->occursin(r"mon$", x), readdir("./rr-data"))
  data = read_plant_file("./rr-data/" * files[1])
  for file in files[2:end]
    pd = read_plant_file("./rr-data/" * file)
    data = vcat(data, pd)
  end
  # categorical variables
  categorical!(data, :steamSystem)
  data[:steamSystem] = recode(data[:steamSystem],
          1=>"BabcockWilcox", 2=>"Combustion",
          3=>"GE", 4=>"Westinghouse",
          5=>"Other1" , 6=>"Other2")
  categorical!(data,:vintage)
  data[:vintage] = recode(data[:vintage],
          0=>"preTMI", 1=>"postTMI")

  data[:hrsOperating] = similar(data[:hrsTotal])
  data[:hrsOperating] .= data[:hrsTotal] .- (data[:hrsRefuel] .+
                                             data[:hrsPlanOut] .+
                                             data[:hrsForcedOut])
  data[:hrsOperating][data[:hrsOperating].<0] .= 0.0
  data[:hrsOperating][data[:numForcedOut].==-1] .= 0.0 # exited plants

  data[:action] = ""
  portionOp = data[:hrsOperating]./data[:hrsTotal]
  data[:action][data[:hrsOperating].<=0] .= "shutdown"
  data[:action][portionOp .> 0.0] .= "run25"
  data[:action][portionOp .> 0.25] .= "run50"
  data[:action][portionOp .> 0.50] .= "run75"
  data[:action][portionOp .> 0.75] .= "run99"
  data[:action][portionOp .> 0.99] .= "run100"
  data[:action][data[:hrsRefuel] .> 0] .= "refuel"
  data[:action][data[:numForcedOut] .== -1] .= "exit"

  data[:spelltype] = "operating"
  data[:spelltype][data[:hrsRefuel] .> 0] .= "refueling"
  data[:spelltype][data[:numForcedOut] .== -1] .= "exit"

  data[:t] = data[:year]*12 + data[:month]
  sort!(data, (:ID,:t))

  data[:nppSignal] = "none"
  data[:nppSignal][data[:hrsForcedOut] .>0.0 ] .= "forcedOut"
  # indexing by boolean expression with missings creates problems,
  # surely there's a better workaround, but whatever
  missing_to_false = x -> ifelse(ismissing(x), false, x)
  idx = (data[:hrsRefuel] .> 0.0) .& missing_to_false.(panel_lag(:hrsRefuel, data, :ID, :t, 1) .> 0.0)
  data[:nppSignal][idx] .= "contRefuel"

  # shift signal by 1 period to match model timing convention (signal
  # and action at time t determines state at time t+1)
  data[:nppSignal] = panel_lag(:nppSignal, data, :ID, :t, -1)
  data[:action] = panel_lag(:action, data, :ID, :t, -1)

  # calculate spell durations
  data[:duration] = 1
  lagspell = panel_lag(:spelltype, data, :ID, :t, 1)
  for r in 2:nrow(data)
    if ( (lagspell[r]==data[r,:spelltype]) &
         (data[r,:ID]==data[r-1,:ID]) &
         (data[r,:t]==(data[r-1,:t]+1)) )
      data[r,:duration] = data[r-1,:duration]+1
    else
      data[r,:duration] = 1
    end
  end

  # create indicator for "major problem spells" (outages lasting longer than 9 months)
  nop = data[:hrsOperating].==0.0
  tmp = (nop .!= lag(nop)) .| (data[:ID] .!= lag(data[:ID]))
  tmp[typeof.(tmp).==Missing] .= false
  opid = cumsum(tmp)
  olength = similar(opid)
  for o in unique(opid)
    d = sum(opid.==o)
    olength[opid.==o] .= d
  end
  data[:majorProblemSpell] = (olength.>9) .& nop

  return(data)
end

""" 
    probfn(state, mod, transform!; coefs)
  
This is specific to this example, so I'm not writing good
documentation.

"""
function probfn(state::AbstractDataFrame, mod, transform!;
                coefs=coef(mod)) 
  if (issubset(mod.mf.terms.terms, names(state)))
    # no need to transform X's
    sc = state
  else
    sc = deepcopy(state) # create a copy, so we don't modify the input
    transform!(sc)
    #out=predict(pof, sc)
    newTerms = StatsModels.dropresponse!(mod.mf.terms)
    # create new model frame/matrix
    mf = ModelFrame(newTerms, sc; contrasts = mod.mf.contrasts)
    newX = ModelMatrix(mf).m
    exb = exp.(newX * coefs)
    p=exb./(1.0 .+ exb)
  end
  if (nrow(state)==1)
    return(p[1])
  else
    return(p)
  end
end

""" 
    transition_prob(poutage, presume, statefn, actionfn, maxduration)

   Computes transition probabilities for all states and actions.
   
   Inputs:
   - `poutage` function that returns probability of forced outage 
   - `presume` function that returns probability of resuming operation
   after refueling
   - `statefn`  function pair for converting between state indices and data
   vectors. Can be reated by `vector_index_converter`.
   - `actionfn` function pair for converting between action indices and data
   vectors. Can be reated by `vector_index_converter`.
   - `maxduration` maximum allowed duration. If
   `oldstate[:duration]+1>maxduration`, then 
   `newstate[:duration] =  maxduration` 

   Output:
   - `P[a][snew, sold]` Array of sparse matrices with P[a][snew,sold]
   equal to the probability of state `snew` given action `a` and state
   `sold`. 
    
"""
function transition_prob(poutage, presume,
                         statefn, actionfn, maxduration)
  elementtype = typeof(poutage(statefn.data(1)))
  P = Array{typeof(sparse(zeros(elementtype,1,1))),1}(undef,
                                                      actionfn.n)
  for a in 1:actionfn.n
    P[a] = sparse(zeros(elementtype,statefn.n, statefn.n))
  end
  for i in 1:statefn.n
    oldstate = statefn.data(i)
    newstate = deepcopy(oldstate)
    # assign probabilities to reachable states given each action

    # exit
    a = actionfn.index("exit")
    newstate[:spelltype] = "exit"
    P[a][statefn.index(newstate),i] = 1.0

    if (oldstate[:spelltype][1] != "exit")
      # choose to refuel
      if (oldstate[:spelltype][1] =="refueling")
        newstate[:duration] = min(oldstate[:duration][1]+1, maxduration)
      else
        newstate[:duration] = 1
      end
      newstate[:spelltype] = "refueling"
      a = actionfn.index("refuel")
      newstate[:nppSignal] = "none"
      P[a][statefn.index(newstate),i] = presume(oldstate)*(1.0-poutage(oldstate))  
      newstate[:nppSignal] = "forcedOut"
      P[a][statefn.index(newstate),i] = presume(oldstate)*poutage(oldstate)

      newstate[:nppSignal] = "contRefuel"
      P[a][statefn.index(newstate),i] = (1-presume(oldstate))

      # choose to operate
      if (oldstate[:nppSignal][1] != "contRefuel")
        newstate = deepcopy(oldstate)
        newstate[:spelltype] = "operating"
        if (oldstate[:spelltype][1] == "refueling")
          newstate[:duration] = 1
        else
          newstate[:duration] = min(oldstate[:duration][1]+1, maxduration)
        end

        for alev in ["shutdown", "run25", "run50", "run75",
                     "run99", "run100"]
          a = actionfn.index(alev)
          newstate[:nppSignal] = "none"
          P[a][statefn.index(newstate),i] = (1-poutage(oldstate))

          newstate[:nppSignal] = "forcedOut"
          P[a][statefn.index(newstate),i] = poutage(oldstate)
        end
      end
    end # oldstate != exit
  end # loop over states
  return(P)
end # function transition_prob

""" 
     payoffs(ϕ, statefn, actionfn)

Returns matrix of payoffs for each state and action.  
"""
function payoffs(ϕ, statefn, actionfn)
  u = zeros(eltype(ϕ), statefn.n, actionfn.n)
  for a in 1:actionfn.n
    u[:,a] .= ϕ[a]
  end
  action = actionfn.data(1)
  for s in 1:statefn.n
    state = statefn.data(s)
    if (state[:spelltype][1]=="exit")
      action[:action] = "exit"
      u[s,:] .= ϕ[actionfn.index(action)]
    else
      if (state[:spelltype][1]=="operating")
        for alev in ["run25", "run50", "run75",
                     "run99", "run100"]
          action[:action] = alev
          u[s, actionfn.index(action)] +=
            state[:duration][1]*ϕ[actionfn.n + 1]
        end
      end
      if (state[:nppSignal][1]=="forcedOut")
        action[:action]="refuel"
        u[s, actionfn.index(action)] += ϕ[actionfn.n+2]
        action[:action]="shutdown"
        u[s, actionfn.index(action)] += ϕ[actionfn.n+3]
        action[:action]="run100"
        u[s, actionfn.index(action)] += ϕ[actionfn.n+4]
      end
    end
  end
  return(u)
end # function payoffs()

""" 
    feasible_actions(statefn,actionfn)

Returns matrix of booleans indicating which actions are feasible in
each state.
"""
function feasible_actions(statefn, actionfn)
  feasible = fill(true, statefn.n, actionfn.n)
  for s in 1:statefn.n
    state = statefn.data(s)
    if state[:nppSignal][1]=="contRefuel"
      feasible[s,:].=false
      feasible[s,actionfn.index("refuel")] = true
    end
    if state[:spelltype][1]=="exit"
      feasible[s,:].=false
      feasible[s,actionfn.index("exit")] = true
    end
  end
  return feasible
end


end # module RustRothwell
