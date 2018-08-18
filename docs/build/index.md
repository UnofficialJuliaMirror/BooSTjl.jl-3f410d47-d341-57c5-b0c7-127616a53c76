<a id='BooSTjl.estimate_derivatives-Tuple{Tuple{Dict{Float64,Array{Float64,2}},Array{Float64,1},Array{Float64,1},Float64,Tuple{Float64,Float64,Int64,Array{Float64,1},Int64,Float64},Int64,Array{Float64,1}},Array{Float64,2},Int64}' href='#BooSTjl.estimate_derivatives-Tuple{Tuple{Dict{Float64,Array{Float64,2}},Array{Float64,1},Array{Float64,1},Float64,Tuple{Float64,Float64,Int64,Array{Float64,1},Int64,Float64},Int64,Array{Float64,1}},Array{Float64,2},Int64}'>#</a>
**`BooSTjl.estimate_derivatives`** &mdash; *Method*.



```
estimate_derivatives(object::Tuple{Dict{Float64,Array{Float64,2}},Array{Float64,1},Array{Float64,1},Float64,Tuple{Float64,Float64,Int64,Array{Float64,1},Int64,Float64},Int64, Vector{Float64}},x::Matrix{Float64},variable::Int64)
```

**`object`** Output object from `BooST` function

Returns the derivative of y with respect to x[:,variable].


<a target='_blank' href='https://github.com/gabrielrvsc/BooSTjl/blob/c6e9e135bc563b3003ce4fd0ce1a916333d0ae79/src/export.jl#L3-L9' class='documenter-source'>source</a><br>

