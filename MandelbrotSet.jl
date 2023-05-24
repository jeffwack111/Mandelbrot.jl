include("RenderFractal.jl")

M = mandelbrot_patch(0.38693644424000584 + 0.09839851183492665im,0.38693644424000806 + 0.09839851183492719im,6)
PA = problem_array(M,4,1000)
ET = escape_time.(PA)
C = colorschemes[:glasbey_bw_minc_20_hue_150_280_n256].colors
pic = apply_color(ET,C)