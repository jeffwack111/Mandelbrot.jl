include("RenderFractal.jl")
include("DynamicRays.jl")

c_rabbit = -0.12256117 + 0.74486177im

z0 = landingpoint(c_rabbit,1//53,100,100,10)
theta_epsilon = 2*pi*0.5
r_epsilon = 0.01
epsilon = r_epsilon*exp(2im*pi*theta_epsilon)
z1 = z0 + epsilon

maxiter = 10

Rescape = 20

P1 = EscapeTimeProblem(z0,z->z*z+c_rabbit,escape(Rescape),maxiter)
P2 = EscapeTimeProblem(z1,z->z*z+c_rabbit,escape(Rescape),maxiter)

T0 = forwardorbit(P1)
T1 = forwardorbit(P2)

L = min(length(T0),length(T1))

T0 = T0[1:L]
T1 = T1[1:L]

fig = Figure()
ax = Axis(fig[1,1])
R = 3
ylims!(-R,R)
xlims!(-R,R)

cmap = :viridis
#cmap = [:red,:yellow]

scatter!(real.(T0),imag.(T0),marker = :circle,markersize = 10,color = 1:L,colormap = cmap, label = "Trajectory of z0")
scatter!(real.(T1),imag.(T1), marker = :cross,markersize = 12,color = 1:L,colormap = cmap, label = "Trajectory of z1")

#Colorbar(fig[1,2],colormap = cmap,limits = (0,L))#,alignmode = Mixed(right = 0))
colsize!(fig.layout, 1, Aspect(1, 1))

J = inverseiterate(c_rabbit,13)
scatter!(real(J),imag(J),markersize = 1,color = "black")
Legend(fig[1,2], ax)

ax = Axis(fig[2,1],yscale = log10,ylabel = "Δ",xlabel = "Iteration count")
rowsize!(fig.layout,2,Relative(0.2))

D = abs.(T0[1:L].-T1[1:L])
plot!(D,color = 1:L,colormap = cmap,label = "Trajectory difference")
bound = [R^n*r_epsilon for n in 1:L]
scatter!(1:L,bound,color = "black",label = "Lipschitz bound")
Legend(fig[2,2], ax)

fig