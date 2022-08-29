using Revise
using Plots
using Measures
using LaTeXStrings


dir = "../../tex/figures/"

Plots.scalefontsizes()
Plots.default(label=nothing, linewidth=0.8, fontfamily="Computer Modern", tickfontsize=6, titlefont="Computer Modern" ,titlefontsize=6,labelfontsize=6,tickfont="Computer Modern", framestyle=:box, grid=nothing)
Plots.scalefontsizes(1.125)


"""
Plots.scalefontsizes()
default(label=nothing, linewidth=0.8, fontfamily="Computer Modern", tickfontsize=5, titlefont="Computer Modern" ,titlefontsize=5,labelfontsize=5,tickfont="Computer Modern", framestyle=:box, grid=nothing)
Plots.scalefontsizes(1.05)



width_pts = 220.0  # or any other value
inches_per_points = 1.0/72.27
width_inches = width_pts *inches_per_points
width_px= width_inches*100;  # or  width_inches*DPI


function plot1D(res::Result; x::String, y::String, sol_type="stable", plot_only=["physical"])

    length(size(res.solutions)) != 1 && error("1D plots of not-1D datasets are usually a bad idea.")
    X = transform_solutions(res, x)
    Y = transform_solutions(res, y) # first transform, then filter
    
    
    for class in plot_only 
        if  class!="binary_labels" 
            Y = HarmonicBalance.filter_solutions.(Y, res.classes[class])
        end
    end
    
    
    p = Plots.plot(real.(getindex.(X, 1)) , real.(getindex.(Y,1)))
    
    Ys, Yu = HarmonicBalance.filter_solutions.(Y, res.classes[sol_type]), HarmonicBalance.filter_solutions.(Y, [.!el for el in res.classes[sol_type]])
    
    p = Plots.plot(real.(getindex.(Ys,1)))
    Plots.plot!(real.(getindex.(Yu,1)), c=1, style=:dash)
    for k in 2:length(Ys[1])
    	data = real.(getindex.(Ys, k))
    	data_u = real.(getindex.(Yu, k))
    	Plots.plot!(data, c=k)
    	Plots.plot!(data_u, c=k, style=:dash)
    end
    
    p
end
"""

