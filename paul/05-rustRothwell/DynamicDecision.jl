"""
Module of functions useful for estimating dynamic decision models.  
"""
module DynamicDecision

export choicevalue, vector_index_converter, loglikei_u,
  emax_mlogit, choice_prob_mlogit

using Base.Iterators, LinearAlgebra, DataFrames, Base.Threads

""" 
    emax_mlogit(x)

  Computes E[max_i(x[i] + ϵ[i])] where each ϵ[i] ~ Logistic()
"""
function emax_mlogit(x::AbstractVector{<:Number})
  # Prevent overflow
  mx = maximum(x)
  log(sum(exp.(x.-mx)))+mx + Base.MathConstants.eulergamma
end

"""
    choice_prob_mlogit(u)
  
  Compute conditional choice probabilities for a discrete choice model
  with choice specific payoffs given by u[i] + ϵ[i] where ϵ[i] ~
  Logistic() (aka Type I extreme value) 
  
  Inputs:
  
  - `u::AbstractVector` expected utility of each choice

  Output:

  - `p =` probability of each choice

"""
function choice_prob_mlogit(u::AbstractVector)
  mu = maximum(u) # prevent overflow
  eu = exp.(u.-mu)
  eu./(sum(eu))
end


"""
    choicevalue(payoff, discount, transprob, T, feasible, emax)

Computes choice specific value function for a finite horizon dynamic
discrete choice problem with Type-1 extreme value errors

Inputs:
  - `payoff` is a states by actions array of period payoffs.
  - `discount` is the discount rate
  - `T` is the time horizon
  - `transprob` is an array of arrays of transition probabilities. `transprob`
     may have either 1 or 2 dimensions. When 1 dimensional,
    `transprob[a][s1,s0]` is the probability of state `s1` given current
    state `s0` and action `a`. When 4 dimensional, `transprob[a,t][s1,s0]` is
    the probability of state `s1` given current state `s0`, action `a` and
    time `t`. In models where many transition probabilities are 0, you
    can have faster computation by storing transprob[a] as sparse.
  - `feasible` is a states by actions array indicating which actions
    are possible in which states
  - `emax` is a function that given a vector, returns 
     E[max_a(v[a] +ϵ[a])] 

Output:
  -  `v[t, s, a]` is choice specific value at time t in state s of
     action a. typeof elements of v will be the same as the return
     type of payoff

"""
function choicevalue(payoff::AbstractArray{<:Number, 2},
                     discount::Number, T::Integer,
                     transprob::AbstractArray{<:AbstractMatrix},
                     feasible::AbstractArray{<:Bool, 2};
                     emax::Function=emax_mlogit)

  # Use multiple dispatch with this helper function  to detect whether
  # transprob is time varying or not.
  function tp(p::AbstractArray{<:AbstractMatrix, 1}, a, t)
    p[a]
  end
  function tp(p::AbstractArray{<:AbstractMatrix, 2}, a, t)
    p[a,t]
  end
  states = size(tp(transprob,1,1))[1]
  actions= length(transprob)

  # careful with types to make sure autodiff works with payoff or
  # trannsprob or both 
  Veltype = eltype(payoff[:,1]'*tp(transprob,1,1))
  v = Array{Veltype, 3}(undef,T,states,actions) 
  vtp1 = Array{Veltype, 2}(undef,states,actions)
  vtp1 .= payoff
  ev = similar(vtp1, states)
  v[T,:,:] .= vtp1
  for t in (T-1):(-1):1
    @threads for s in 1:states
      ev[s] = emax(vtp1[s, feasible[s,:]] )
    end
    @threads for a in 1:actions
      v[t,:,a] = payoff[:,a] .+ (discount * (ev' * tp(transprob,a,t)))'
    end
    vtp1 .= v[t,:,:]
  end
  @threads for s in 1:states
    for a in 1:actions
      if !(feasible[s,a])
        v[:,s,a] .= -Inf
      end
    end
  end      
  return(v)
end # function choicevalue


"""
    vector_index_converter(data, vars)

Given a vector of variables, `vars` that define state, return
functions that given an index, returns a vector of values for the state
variables, and, conversely, given a vector return an index.

Inputs:

  - `data` DataFrame containing state variables.
  - `vars` Array{1} of symbols for each state
     variable. The state variables are assumed to be discrete.
     All levels of each variable should be present in `data`, but all
     combinations of values need not be
     present. `product(unique(data[vars]))` is used to determine the
     possible levels of the state variables.

Output:
  - Function to convert between indexes and state representation.
       - out.data(index) returns a single row dataframe for the
         state with a given index
       - out.index(dataframe) returns the indices for the states in
         dataframe.

"""
function vector_index_converter(data::AbstractDataFrame,
                                vars::AbstractArray{Symbol,1})
  valfn = v->unique(skipmissing(data[v]))
  tstates = collect(product(valfn.(vars)...))[:]
  states = map(s -> DataFrame(map(x->[x], collect(s)), vars), tstates)
  # states is an array of single row data frames
  stateDict = Dict(states[i]=>i for i = 1:length(states))
  function datarow(index::Integer)
    return(deepcopy(states[index]))
  end
  function index(state::AbstractDataFrame)
    sdf = state[vars]
    function anymissing(sdf)
      for v in vars
        if (ismissing(sdf[v][1]))
          return(true)
        end
      end
      return(false)
    end

    if (size(sdf)[1]==1)
      if (anymissing(sdf))
        return(missing)
      else
        return(stateDict[sdf])
      end
    else
      out = Array{Union{Missing,Int64},1}(undef,size(state)[1])
      for r in 1:size(sdf)[1]
        if (anymissing(sdf[r,:]))
          out[r] = missing
        else
          out[r]=stateDict[sdf[r,:]]
        end
      end
      return(out)
    end
  end
  return(data=datarow,index=index, n=length(states))
end # function vector_index_converter

"""
     logfinite(cut)

  Returns function that = log(x) if x>=cut, first order taylor
  expansion of log(x)around cut if x<cut
"""
function logfinite(cut::Number)
  lc = log(cut)
  dlc = 1.0/cut  
  (x)->ifelse(x>=cut, log(x), lc + (x-cut)*dlc)
end

"""
    loglikei_u(u, data, time, stateidx, actionidx, vfn, vargs...; cpfn)

  loglikehood for observation 

  Inputs:
  - `u` payoff for each choice in each state
     Must be indexable by u[stateidx[i], actionidx[i]]
  - `data` dataframe containing time, state, and action indices
  - `time` varible in `data` containing time index
  - `stateidx` varible in `data` containing state index
  - `actionnidx` varible in `data` containing action index
  - `vfn` function computes choice specific value functions. The
     arguments for vfn are vfn(u,vargs...). It should
     return an array of size 
     `[maximum(data[time]),maximum(stateidx), maximum(data[actionidx]`
     If `v = vfn(u,vargs...)`, then `v[t,s,a]` should be choice specific
     value function at time `t` in state `s` for action `a`.
  - `vargs` additional arguments for vfn.
  - `cpfn` function giving conditional choice
     probabilities. If x is a vector of choice specific values
     associated with each action in a single state, ccp(x) should
     return a vector of probabilities of each action. 
  - `logfn` function to use in place of log to avoid log(0)=-Inf

  Output: 
  - `loglikei` vector of length nrow(data) giving log likelihood of each
     observation. If data for observation i is missing we set,
     loglikei[i] = missingval
"""
function loglikei_u(u::AbstractMatrix,
                    data::AbstractDataFrame,
                    time::Symbol,
                    stateidx::Symbol,
                    actionidx::Symbol,
                    vfn::Function,
                    vargs...;
                    cpfn::Function=choice_prob_mlogit,
                    logfn::Function=logfinite(1e-300),
                    missingval=missing)
  v = vfn(u,vargs...)
  ccp = similar(v)
  @threads for t in 1:size(v)[1]
    for s in 1:size(v)[2]
      ccp[t,s,:] = cpfn(v[t,s,:])
    end
  end
  #ccp = mapslices(cpfn, v, dims=[3])
  if (ismissing(missingval))
    T = Union{Missing, eltype(ccp)}
  else
    T = eltype(ccp)
  end
  loglikei = Array{T,1}(undef,nrow(data))  
  @threads for i in 1:nrow(data)
    if (ismissing(data[time][i]) ||
        ismissing(data[stateidx][i]) ||
        ismissing(data[actionidx][i]))
      loglikei[i] = missingval
    else
      loglikei[i] = logfn.(ccp[data[time][i], data[stateidx][i],
                               data[actionidx][i]])
    end
  end
  return(loglikei)
end



end # module dynamicDecision


