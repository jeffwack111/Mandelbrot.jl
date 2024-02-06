include("RenderFractal.jl")

# Sensitivity to initial conditions
# you arent updating your lipz shit
function delta(z0,z1,c,iter)
    delta = Float64[]
    lipsch = Float64[]
    L = lipschitz(z0,abs(z1-z0))
    for ii in 1:iter
        push!(delta,abs(z1-z0))
        push!(lipsch,L^ii)
        z0 = z0*z0 + c
        z1 = z1*z1 + c
    end
    return delta,lipsch
end

function lipschitz(z,epsilon)
    return 2*(abs(z)+epsilon)
end

