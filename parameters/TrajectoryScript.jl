include("RenderFractal.jl")
include("DynamicRays.jl")

c_rabbit = -0.12256117 + 0.74486177im

z0 = landingpoint(c_rabbit,1//5,100,100,10)
theta_epsilon = 2*pi*0.2
r_epsilon = 0.01
epsilon = r_epsilon*exp(2im*pi*theta_epsilon)
z1 = z0 + epsilon

maxiter = 500

P1 = EscapeTimeProblem(z0,z->z*z+c_rabbit,escape(1000),maxiter)
P2 = EscapeTimeProblem(z1,z->z*z+c_rabbit,escape(1000),maxiter)

T0 = forwardorbit(P1)
T1 = forwardorbit(P2)

fig = Figure()
ax = Axis(fig[1,1],aspect = 1)
R = 5
ylims!(-R,R)
xlims!(-R,R)
scatter!(real.(T0),imag.(T0),marker = :circle,markersize = 10,color = 1:length(T0),colormap = [:red,:yellow])
scatter!(real.(T1),imag.(T1), marker = :cross,markersize = 12,color = 1:length(T1),colormap = [:red,:yellow])

ax = Axis(fig[2,1],yscale = log10)
rowsize!(fig.layout,2,Relative(0.2))


l = min(length(T0),length(T1))

D = abs.(T0[1:l].-T1[1:l])
plot!(D)

fig