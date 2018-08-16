function estimate_derivatives(object::Tuple{Dict{Float64,Array{Float64,2}},Array{Float64,1},Array{Float64,1},Float64,Tuple{Float64,Float64,Int64,Array{Float64,1},Int64,Float64},Int64, Vector{Float64}},x::Matrix{Float64},variable::Int64)::Vector{Float64}
    M::Int64 = object[5][5]
    N::Int64 = size(x,1)
    v::Float64 = object[5][1]
    rho::Vector{Float64} = object[7]
    vrho::Vector{Float64} = v*rho
    peaux::Matrix{Float64} = Matrix{Float64}(N,M)
    derivative::Vector{Float64} = Vector{Float64}(N)
    model::Dict{Float64,Array{Float64,2}} = object[1]
    for i in 1:M
        peaux[:,i] = gradient_st(model[i],x,variable)        
    end
    derivative = peaux*vrho
    return derivative
end

function BoostMore(x::Matrix{Float64},y::Vector{Float64},object::Tuple{Dict{Float64,Array{Float64,2}},Array{Float64,1},Array{Float64,1},Float64,Tuple{Float64,Float64,Int64,Array{Float64,1},Int64,Float64},Int64, Vector{Float64}}, M::Int64, display::Bool)::Tuple{Dict{Float64,Array{Float64,2}},Array{Float64,1},Array{Float64,1},Float64,Tuple{Float64,Float64,Int64,Array{Float64,1},Int64,Float64},Int64,Vector{Float64}}

    (v::Float64,p::Float64,d_max::Int64,gamma::Vector{Float64},Mold::Int64,node_obs::Int64) = object[5]
    params::Tuple{Float64,Float64,Int64,Vector{Float64},Int64,Float64} = (v,p,d_max,gamma,M+Mold,node_obs)
    d_max = d_max - 1

    save_rho::Vector{Float64} = Vector{Float64}(Mold+M)
    save_rho[1:Mold] = object[7]
    ybar::Float64 = object[4]
    N::Int64 = length(y)
    phi::Vector{Float64} = predictBoost(object,x)
    brmse::Vector{Float64} = Vector{Float64}(Mold+M)
    save_tree::Dict{Float64,Array{Float64,2}} = object[1]
    brmse[1:Mold] = object[3]
    sizehint!(save_tree,Mold+M)
    for i in (Mold+1):(M+Mold)
        u::Vector{Float64} = y .-phi 
        (tr::Matrix{Float64},fitted::Vector{Float64}) = grow_tree(x,u,p,d_max,gamma,node_obs)::Tuple{Array{Float64,2},Array{Float64,1}}
        rho::Float64 = sum([fitted[j]*u[j] for j in 1:N])/sum([fitted[j]*fitted[j] for j in 1:N])
        [phi[j] = phi[j] + v*rho*fitted[j] for j in 1:N]
        save_tree[i] = tr
        save_rho[i] = rho
        sqerror::Vector{Float64} = [(y[j]-phi[j])^2 for j in 1:N]  
        rmse::Float64 = sqrt(mean(sqerror))
        brmse[i] = rmse
        if display == true
            println("iteration ",i, " rmse = ", rmse)
        end
    end
    nvar::Int64 = size(x,2)
    
    return save_tree, phi, brmse, ybar, params, nvar, save_rho

end



function predictBoost(object::Tuple{Dict{Float64,Array{Float64,2}},Array{Float64,1},Array{Float64,1},Float64,Tuple{Float64,Float64,Int64,Array{Float64,1},Int64,Float64},Int64, Vector{Float64}},newx::Matrix{Float64})::Vector{Float64}
    v::Float64 = object[5][1]
    y0::Float64 = object[4]
    rho::Vector{Float64} = object[7]
    model::Dict{Float64,Array{Float64,2}} = object[1]
    rhov::Vector{Float64} = rho*v
    N::Int64 = size(x,1)
    M::Int64 = object[5][5]
    fittedaux::Matrix{Float64} = Matrix{Float64}(N,M)
    fitted::Vector{Float64} = Vector{Float64}(N)
    for i in 1:M
        fittedaux[:,i] = eval_tree(newx,model[i])        
    end
    fitted = y0 + fittedaux*rhov
    return fitted
end



function Boost(x::Matrix{Float64},y::Vector{Float64},v::Float64, p::Float64, d_max::Int64, gamma::Vector{Float64}, M::Int64, display::Bool, node_obs::Int64)::Tuple{Dict{Float64,Array{Float64,2}},Array{Float64,1},Array{Float64,1},Float64,Tuple{Float64,Float64,Int64,Array{Float64,1},Int64,Float64},Int64,Vector{Float64}}

    params::Tuple{Float64,Float64,Int64,Vector{Float64},Int64,Float64} = (v,p,d_max,gamma,M,node_obs)
    d_max::Int64 = d_max-1
    N::Int64 = length(y)
    ybar::Float64 = mean(y)
    phi::Vector{Float64} = ones(N)*ybar
    nvar::Int64 = size(x,2)

    brmse::Vector{Float64} = Vector{Float64}(M)
    save_rho::Vector{Float64} = Vector{Float64}(M)
    save_tree::Dict{Float64,Matrix{Float64}} = Dict{Float64,Matrix{Float64}}()
    sizehint!(save_tree,M)

    for i in 1:M
        u::Vector{Float64} = y .-phi 
        (tr::Matrix{Float64},fitted::Vector{Float64}) = grow_tree(x,u,p,d_max,gamma,node_obs)::Tuple{Array{Float64,2},Array{Float64,1}}
        rho::Float64 = sum([fitted[j]*u[j] for j in 1:N])/sum([fitted[j]*fitted[j] for j in 1:N])
        [phi[j] = phi[j] + v*rho*fitted[j] for j in 1:N]
        save_tree[i] = tr
        save_rho[i] = rho
        sqerror::Vector{Float64} = [(y[j]-phi[j])^2 for j in 1:N]  
        rmse::Float64 = sqrt(mean(sqerror))
        brmse[i] = rmse
        if display == true
            println("iteration ",i, " rmse = ", rmse)
        end

    end
    return save_tree, phi, brmse, ybar, params, nvar, save_rho
end

