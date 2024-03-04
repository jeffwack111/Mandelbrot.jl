using CairoMakie

include("OldHubbardTrees.jl")

function orientedtree((IA,angles))
    (E,F,markedpoints) = embedtree((IA,angles))
    
    root = findall(x->x==1, F)[1]

    T = [[root],E[root]]
 
    nadded = 1 + length(T[2])
    n = length(E)

    while nadded < n
        parents = T[end]
        children = []
        for vertex in parents
            for u in E[vertex]
                if !(u in T[end-1])
                    push!(children,u)
                end
            end
        end
        push!(T,children)
        nadded += length(children)
    end

    
    Ngens = length(T)

    X = zeros(Float64,n)
    Y = zeros(Float64,n)
    
    for (gen,vertices) in enumerate(T)
        k = length(vertices)
        for (ii,u) in enumerate(vertices)
            X[u] = ii/(k+1)
            Y[u] = 1 - gen/(Ngens+1)
        end
    end
    pos = Point.(X,Y)
    
    f = Figure()

    ax1 = Axis(f[1, 1],aspect = 1)
    hidedecorations!(ax1)
    hidespines!(ax1)

    for (ii,p) in enumerate(E)
        for n in p
            lines!([pos[ii],pos[n]])
        end
    end

    scatter!(pos)
    tex = ["$ii" for ii in 1:n]
    text!(pos,text = tex)
    
    return f

end

function orientedtree(theta::Rational)  
    return orientedtree(angledinternaladdress(theta))
end

function arbtree((E, F,markedpoints))
    R = [collect(x) for x in E]
    z = findall(x->x==1,F)[1]
    for (point,arms) in enumerate(R)

        if point !==z && length(arms) > 1
            SECONDGLARM = globalarms(R,point)

            y = findthe(z,SECONDGLARM)
            circshift!(R[point],-y)
        end

    end
    root = findall(x->x==1, F)[1]

    T = [[root],E[root]]

    nadded = 1 + length(T[2])
    n = length(E)

    while nadded < n
        parents = T[end]
        children = []
        for vertex in parents
            for u in E[vertex]
                if !(u in T[end-1])
                    push!(children,u)
                end
            end
        end
        push!(T,children)
        nadded += length(children)
    end

    
    Ngens = length(T)

    X = zeros(Float64,n)
    Y = zeros(Float64,n)
    
    for (gen,vertices) in enumerate(T)
        k = length(vertices)
        for (ii,u) in enumerate(vertices)
            X[u] = ii/(k+1)
            Y[u] = 1 - gen/(Ngens+1)
        end
    end
    pos = Point.(X,Y)
    
    f = Figure()

    ax1 = Axis(f[1, 1],aspect = 1)
    hidedecorations!(ax1)
    hidespines!(ax1)

    for (ii,p) in enumerate(E)
        for n in p
            lines!([pos[ii],pos[n]])
        end
    end

    scatter!(pos)
    tex = ["$ii" for ii in 1:n]
    text!(pos,text = tex)
    
    return f
end

    
function arbtree(intadd::Vector{Int})
    return arbtree(hubbardtree(intadd))
end