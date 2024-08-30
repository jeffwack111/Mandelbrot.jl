using GLMakie
include("AngleDoubling.jl")

#for each internal address, calculate a screen which shows the future options

function wedge!(ax,start,finish,innerrad,outerrad)
    distance = finish - start
    vertex_per_deg = 1
    nvertices = max(2, ceil(Int, rad2deg(distance) * vertex_per_deg))

    outerpoints = map(LinRange(start, finish, nvertices)) do rad
                Point2(cos(rad) * outerrad, sin(rad) * outerrad)
            end
    innerpoints = map(LinRange(finish, start, nvertices)) do rad
                Point2(cos(rad) * innerrad, sin(rad) * innerrad)
            end

    points = append!(innerpoints,outerpoints)

    poly!(ax,points)
end

function drawwedges!(ax,AIA)
    N = AIA.addr[end]
    n = 20
    for ii in N+1:N+n 
        k = denominator(AIA,ii)
        angles = []
        wedges = []
        angletext = []
        for jj in 1:k-1
            if gcd(jj,k) == 1
                push!(angles,(jj-0.5)*2*pi/(k-1))
                push!(angletext, "$jj/$k $ii")
                push!(wedges,((jj-1)*2*pi/(k-1),(jj)*2*pi/(k-1),2.0^(-ii+1),2.0^-ii))
            end
        end
        r = (2.0^-ii+2.0^(-ii+1))/2
        centers = [(r*cos(theta),r*sin(theta)) for theta in angles]
        for wedge in wedges
            wedge!(ax,wedge...)
        end
        text!(ax,centers,text = angletext, align = (:center,:center))
    end
end

function phonebook!(gl,AIA)
    label = Label(gl[2,1],repr(AIA[]))
    rowsize!(gl, 1, Aspect(1, 1.0))
    ax1 = Axis(gl[1,1],aspect = 1)
    drawwedges!(ax1,AIA[])

    hidedecorations!(ax1)
    hidespines!(ax1)
    deactivate_interaction!(ax1, :rectanglezoom)
    deactivate_interaction!(ax1, :dragpan)

    register_interaction!(ax1, :select) do event::MouseEvent, ax1
        if event.type === MouseEventTypes.leftdown
            ang = atan(event.data[2],event.data[1])
            rad = sqrt(sum([x*x for x in event.data]))
            idx = -1*Int(floor(log2(rad)))
            if idx > 0
                k = denominator(AIA[],idx)
                num = Int(mod1(ceil(ang*(k-1)/(2*pi)),k-1))
                println("$idx and $num / $k")
                if gcd(num,k) == 1
                    newaddr = push!(copy(AIA[].addr),idx)
                    newangles = push!(copy(AIA[].angles),Rational(num,k))
                    AIA[] = AngledInternalAddress(newaddr,newangles)
                    println(AIA[])
                    empty!(ax1)
                    drawwedges!(ax1,AIA[])
                    label.text = repr(AIA[])
                end
            end
        end
    end

    

end


fig = Figure()
ga = fig[1, 1] = GridLayout()

AIA = Observable(AngledInternalAddress([1],[]))

phonebook!(ga,AIA)

fig