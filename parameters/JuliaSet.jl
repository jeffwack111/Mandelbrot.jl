include("RenderFractal.jl")

z_rabbit = -0.12256117 + 0.74486177im

J = julia_patch(0.0+0.0im,2.0+0.0im)
f(z) = z*z + z_rabbit
PA = jproblem_array(J,f,escape(4),100)
#C = colorschemes[:glasbey_bw_minc_20_hue_150_280_n256].colors
#colors = fill(C,size(PA))
pic = [x[1] for x in escapetime.(PA)]

fig = Figure()
scene = Scene(camera=campixel!, resolution=(1000,1000))
ax = Axis(scene,aspect = 1)
hidedecorations!(ax)
hidespines!(ax)

heatmap!(scene,mod.([x[1] for x in escapetime.(PA)],50),nan_color = RGBAf(0,0,0,1),colormap = :Set3_12)

scene