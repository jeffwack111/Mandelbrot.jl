using CairoMakie

f = Figure()
ax = Axis(f[1, 1],aspect = 1)
xlims!(ax,-2,2)
ylims!(ax,-2,2)
hidedecorations!(ax)
arc!(Point2f(0.0,0.0), 1.0, 0.0, 2*π, color=:black)

function plot_leaf!(leaf)
    a = leaf[1]
    b = leaf[2]
    θ = (a+b)/2
    if abs2(a-b) != 0.25
        δ = (b-a)/2
        r = sec(2*π*δ)
        center = r*exp(2*π*1im*θ)
        radius = abs(exp(2*π*1im*a)-center)
        arc!(Point2f(real(center),imag(center)),radius,-π,π)
        #arc!(Point2f(real(center),imag(center)),radius,2*π*(θ + 0.25 + δ),2*π*(θ + 0.75 - δ))
    else
        lines!([real(exp(2*π*1im*a)),real(exp(2*π*1im*b))],[imag(exp(2*π*1im*a)),imag(exp(2*π*1im*b))])
    end
end

function preimages(leaf, α)
    #assume that leaf[2] is counterclockwise of leaf[1]
    #always output pairs with the same orientation

    a = leaf[1]/2
    b = leaf[2]/2
    c = leaf[1]/2+1//2
    d = leaf[2]/2+1//2

    if α/2<=a #then α/2<a<1/2
        if b<α/2
            return [(a,d),(c,b)]
        else
            return [(a,b),(c,d)]
        end
    else #then 0<a<α/2
        if b<α/2
            return [(a,b),(c,d)]
        else
            return [(a,d),(c,b)]
        end
    end
end

function laminate(alpha,depth)
    if depth == 1
        leaf = (alpha/2,alpha/2+1//2)
        plot_leaf!(leaf)
        return [leaf]
    else
        leaves = vcat([preimages(leaf,alpha) for leaf in laminate(alpha,depth-1)]...)
        for leaf in leaves
            plot_leaf!(leaf)
        end
        return leaves
    end
end

laminate(9//240,10)

f
