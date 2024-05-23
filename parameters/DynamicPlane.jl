include("DynamicRays.jl")
include("RenderFractal.jl")
include("../spiders/Spiders.jl")
include("../trees/Lamination.jl")
include("../trees/ShowTree.jl")

function characteristicrays(theta::Rational) 

    c = spideriterate(theta,100)

    (OrientedHubbardTree,Boundary) = orient(theta)
    (EdgeList,Nodes) = adjlist(OrientedHubbardTree)
    CharacterSet = characteristicset(OrientedHubbardTree)
    OZ = labelonezero(OrientedHubbardTree,Boundary)

    angles = [first(anglesof(OZ,(OrientedHubbardTree,Boundary),node)) for node in CharacterSet]
    n = period(theta)
    orbit = [0.0+0.0im]
    for ii in 1:n-1
        push!(orbit,orbit[end]^2+c)
    end
    
    thetaprime = conjugate(theta)

    rays = dynamicrays(c,theta,100,10,20)
    conjrays = dynamicrays(c,thetaprime,100,10,20)

    branchrays = Dict()

    for branchangle in angles
        merge!(branchrays,dynamicrays(c,branchangle,100,10,20))
    end

    raydict = merge(branchrays,rays,conjrays)
    println("a")
    feetlist = [[raydict[angle][end] for angle in anglesof(OZ,(OrientedHubbardTree,Boundary),node)] for node in Nodes]
    println("b")
    zvalues = [sum(feet)/length(feet) for feet in feetlist]

    pos = Point.(real.(zvalue),imag.(zvalues))

    julia = inverseiterate(c,22)

    fig = showtree(OrientedHubbardTree,(EdgeList,Nodes),OZ,pos)
    ax = Axis(fig[1, 1],aspect = 1)
    r = 2
    xlims!(-r,r)
    ylims!(-r,r)
    hidedecorations!(ax)
    #hidespines!(ax)

    onlyrays = [pair[2] for pair in raydict]
    plotrays!(ax,onlyrays)

    scatter!(ax,real(julia),imag(julia),markersize = 1,color = "black")
    #scatter!(ax,real(orbit),imag(orbit))

    return fig
end