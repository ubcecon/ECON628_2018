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

5. Subsetting data with logicals when missings occur. I find that I
often want to subset data by a logical expression, e.g.
```julia
data[:x][data[:y].>0] = newvalue # error is any data[:y] is missing
```
My approach has been
```julia
missing_to_false = x -> ifelse(ismissing(x), false, x)
data[:x][missing_to_false.(data[:y].>0)] = newvalue # error is any data[:y] is missing
```
This is fine, but somehow feels clunky. 