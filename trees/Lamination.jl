using CairoMakie
using ColorSchemes

include("../sequences/AngleDoubling.jl")

#The Structure of the Set of Quadratic Invariant Laminations

struct Leaf
    a::Rational
    b::Rational
    #a is less than b and both are between zero and one
    function Leaf(x::Rational,y::Rational)
        if x<=y
            new(x,y)
        else
            new(y,x)
        end
    end
end

import Base.==
function ==(x::Leaf,y::Leaf)
    return x.a==y.a && x.b==y.b
end

function plot_leaf!(leaf::Leaf,color)
    a = leaf.a
    b = leaf.b
    θ = (a+b)/2
    if abs2(a-b) != 0.25
        δ = (b-a)/2
        r = sec(2*π*δ)
        center = r*exp(2*π*1im*θ)
        radius = abs(exp(2*π*1im*a)-center)
        arc!(Point2f(real(center),imag(center)),radius,-π,π, color = color)
        #arc!(Point2f(real(center),imag(center)),radius,2*π*(θ + 0.25 + δ),2*π*(θ + 0.75 - δ))
    else
        lines!([real(exp(2*π*1im*a)),real(exp(2*π*1im*b))],[imag(exp(2*π*1im*a)),imag(exp(2*π*1im*b))],color = color)
    end
end



function image(leaf::Leaf)
    c = (leaf.a*2)%1
    d = (leaf.b*2)%1
    return Leaf(c,d)
end

function orbit(leaf::Leaf)
    items = Leaf[]

    while isempty(findall(x->x==leaf,items))
        push!(items,leaf)
        leaf = image(leaf)
    end

    preperiod = findall(x->x==leaf,items)[1] - 1

    return Sequence(items,preperiod)
    
end


function pre_images(leaf::Leaf, α)
    #assume that leaf[2] is counterclockwise of leaf[1]
    #always output pairs with the same orientation

    a = leaf.a/2
    b = leaf.b/2
    c = leaf.a/2+1//2
    d = leaf.b/2+1//2

    if α/2<=a #then α/2<a<1/2
        if b<α/2
            return [Leaf(a,d),Leaf(c,b)]
        else
            return [Leaf(a,b),Leaf(c,d)]
        end
    else #then 0<a<α/2
        if b<α/2
            return [Leaf(a,b),Leaf(c,d)]
        else
            return [Leaf(a,d),Leaf(c,b)]
        end
    end
end

#Proposition II.4.5, page 69
function saturation(alpha)
    leaves = [[Leaf(alpha/2,alpha/2+1//2)]]
    for j = 1:9
        children = copy(leaves[end])
        parents = vcat([pre_images(leaf,alpha) for leaf in children]...)
        push!(leaves,parents)
    end
    return leaves
end

function show_lamination(alpha)
    f = Figure()
    ax = Axis(f[1, 1],aspect = 1)
    xlims!(ax,-2,2)
    ylims!(ax,-2,2)
    hidedecorations!(ax)

    arc!(Point2f(0.0,0.0), 1.0, 0.0, 2*π, color=:black)
    L = lamination(alpha)
    C = colorschemes[:southwest].colors
    for G in enumerate(L)
        color = C[G[1]]
        for leaf in G[2]
            plot_leaf!(leaf,color)
        end
    end
    return f
end

function linked(X::Leaf,Y::Leaf)
    if X.b < Y.a 
        return false
    elseif Y.b < X.a
        return false
    elseif X.a < Y.a && Y.a < X.b
        if Y.b < X.b 
            return false
        else
            return true
        end
    elseif Y.a < X.a && X.a < Y.b
        if X.b < Y.b 
            return false
        else
            return true
        end
    else
        error("you missed a case")
    end
end

function conjugate(angle::Rational)
    orb = orbit(angle)
    n = period(orb)

    halfa = angle/2
    halfb = angle/2 + 1//2

    if halfa in orb.items
        theta = [halfb]
    elseif halfb in orb.items
        theta = [halfa]
    else
        error("the orbit should contain one half")
    end

    criticalleaf = Leaf(halfa,halfb)

    for k in 2:n
        a = theta[k-1]/2
        b = theta[k-1]/2 + 1//2

        if ! linked(criticalleaf,Leaf(a,orb[n-k+1]))
            push!(theta,a)
        elseif ! linked(criticalleaf, Leaf(b,orb[n-k+1]))
            push!(theta,b)
        else
            error("one leaf should be unlinked")
        end
    end

    return angle + (theta[n] - angle)//(1 - 1//(2^n))

end

function hubbardtree(angle::Rational)
    conjugateangle = conjugate(angle)
    leaf = Leaf(angle,conjugateangle)
    orb = orbit(leaf)

    C = colorschemes[:southwest].colors

    f = Figure()
    ax = Axis(f[1, 1],aspect = 1)
    xlims!(ax,-2,2)
    ylims!(ax,-2,2)
    hidedecorations!(ax)
    arc!(Point2f(0.0,0.0), 1.0, 0.0, 2*π, color=:black)
    for (ii,L) in enumerate(orb.items)
        plot_leaf!(L,C[mod1(ii,10)])
    end

    return f
end