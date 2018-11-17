## Questions about Julia

1. Specifying argument types. Jesse cautions against this. I think the
   reasons are: (i) choosing types is hard and an extra complication
   (ii) choosing an overly specific type unneccessarily limits the
   generality of your code.
   
   I like specifying argument types. I view specifying types as a
   first step in error checking and documentation. Knowing the most
   general abstract type that guarantees the operations in your
   function are all well-defined is, admittedly, a challenge, but I 
   think it's worth the trouble to learn. Why am I wrong? 

> JP: You are not wrong at all.
> * The  caution is against doing it for introductory users (outside of dispatch) is because it requires a lot more maturity in general programming than they might realize, causes harm if they do it wrong, and never helps performance if they do it right.  
> * For advanced users, there are plenty of good reasons to correctly put in abstract types (although it is nice not to have to do it). 
> * Putting in abstract types for functions is generally innocuous once you figure it out.  Creating your own `struct` and properly parameterizing it is very easy (even for advanced users) to make mistakes.


2. I create a DataFrame from a file, and all columns are type
   `Any`. To make the DataFrame more self documenting and for efficiency
   reasons, I think I want more specific types. Is that correct? (Maybe
   it's better to deal with types later.)

   If it is a good idea to assign types to the DataFrame, is the
   following a good way to do so?
```julia
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
```
> JP: I agree that for dataframes not generated from named tuples, adding types is a decent idea.
> * The reason is that they are always concrete, and you don't have the issues of picking the wrong abstract type.
> * You also want to avoid having columns with ``Any``
> * I will try to add something to the notes, but see  http://juliadata.github.io/DataFrames.jl/latest/man/getting_started.html#Constructing-Row-by-Row-1

3. Suppose I have a function for which I want a small portion to do
   something different depending on the type of the inputs. What's the
   preferred way to accomplish this? One option is to make two
   functions for the types of the inputs, but then 90+% of these
   functions are identical, so it would copy many lines of code and be
   difficult to maintain. Another option is to create two versions of
   a helper function that just contains the part of code that
   differs. For example, in the code below, I want to allow
   `transprob` to be either 4 or 3 dimensional. Is this good style?
   Will it be efficient? 
 
> JP: You would need to benchamark it to know (using `@btime`) but my guess is that there is a slightly better way to do it
> * The general rule is you want to let the compiler specialize a compilation at the whole function level, Otherwise, it may need to dynamically dispatch within a function.
> * It is tough to know what the compiler can figure out, though
> * Would it work if you just moved the `tp` function outside of `choicevalue`?  If so, then the compiler would use the specialized function since `transprob` is a single type.
> * There is another few hints on creating temporaries with less code using `similar` and `eltype`.  Let me try (wihtout testing it)
```julia
# transprob is time varying or not. 
function tp(p::AbstractArray{Number, 3}, s, a, t)
   p[:,s,a]
end
function tp(p::AbstractArray{Number, 4}, s, a, t)
   p[:,s,a,t]
end
function choicevalue(payoff, discount, T, transprob, feasible)
  states = size(transprob)[1]
  actions= size(transprob)[3]
  v = similar(payoff, 3)
  vtp1 = zeros(eltype(v), states, actions) # though you don't use it?  would similar work?
  for t in T:(-1):1
    for s in 1:states
      for a in 1:actions
        v[t,s,a] = payoff[s,a] + discount *
          log(sum(exp(vtp1).*feasible, dims=2))*tp(transprob,s,a,t)
      end
    end
  end
  return(v)
end
```

Your own

```julia
function choicevalue(payoff, discount, T, transprob, feasible)

  # Use multiple dispatch with this helper function  to detect whether
  # transprob is time varying or not. 
  function tp(p::AbstractArray{Number, 3}, s, a, t)
    p[:,s,a]
  end
  function tp(p::AbstractArray{Number, 4}, s, a, t)
    p[:,s,a,t]
  end
  
  states = size(transprob)[1]
  actions= size(transprob)[3]
  v = Array{typeof(payoff[1,1]), 3}(undef,T, states,actions)
  vtp1 = zeros(typeof(v[1,1,1]), states, actions)
  for t in T:(-1):1
    for s in 1:states
      for a in 1:actions
        v[t,s,a] = payoff[s,a] + discount *
          log(sum(exp(vtp1).*feasible, dims=2))*tp(transprob,s,a,t)
      end
    end
  end
  return(v)
end
```


4. A nice feature of R's formulas is that you can include
transformations of variables, like 
```r
y ~ log(x) + I(x^2) + as.factor(z)
```
Is this possible with Julia formulas? It seems like no, other than
adding interactions with * and &.
> JP: the GLM guys told me that this will be possible when https://github.com/JuliaStats/StatsModels.jl/pull/71 is merged, but it might take a little while.  Watch it!

5. Subsetting data with logicals when missings occur. I find that I
often want to subset data by a logical expression, e.g.
```julia
data[:x][data[:y].>0] = newvalue # error if any data[:y] is missing
```
> JP: I will ask on slack and get back to you.
My approach has been
```julia
missing_to_false = x -> ifelse(ismissing(x), false, x)
data[:x][missing_to_false.(data[:y].>0)] = newvalue 
```
This is fine, but somehow feels clunky. 
