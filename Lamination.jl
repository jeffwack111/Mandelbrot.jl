using CairoMakie

f = Figure()
ax = Axis(f[1, 1],aspect = 1)
xlims!(ax,-2,2)
ylims!(ax,-2,2)
hidedecorations!(ax)
arc!(Point2f(0.0,0.0), 1.0, 0.0, 2*π, color=:black)

function leaf(u,v)
    if u<v
        a = u
        b = v
    else
        a = v
        b = u
    end
    θ = (a+b)/2
    if abs2(a-b) != 0.25
        δ = (b-a)/2
        r = sec(2*π*δ)
        center = r*exp(2*π*1im*θ)
        radius = abs(exp(2*π*1im*a)-center)
        arc!(Point2f(real(center),imag(center)),radius,2*π*(θ + 0.25 + δ),2*π*(θ + 0.75 - δ))
    else
        lines!([real(exp(2*π*1im*a)),real(exp(2*π*1im*b))],[imag(exp(2*π*1im*a)),imag(exp(2*π*1im*b))])
    end
end

intersecting_circle(1//3,2//3)

f