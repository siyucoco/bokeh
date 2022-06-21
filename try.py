from scipy.integrate import solve_ivp
from bokeh.layouts import column
from bokeh.models import ColumnDataSource, CustomJS, Slider
from bokeh.plotting import figure, output_file, show, row
import matplotlib.pyplot as plt
import numpy as np

def deriv(t, y):
    x = y
    xdot = .04 * x 
    return xdot

t0, tf = 0, 100
y0 = 1, 0, 0
soln = solve_ivp(deriv, (t0, tf), y0)
# print((soln.t))
# print(soln.y.size)

# for i in range(soln.y.shape[0]):
#     plt.plot(soln.t, soln.y[i], label=f'$X_{i}(t)$')
# plt.xlabel('$t$') # the horizontal axis represents the time 
# plt.legend() # show how the colors correspond to the components of X
# plt.show()

source = ColumnDataSource(data=dict(x=soln.t, y=soln.y[0]))

plot = figure(y_range=(0, 10), width=400, height=400)

plot.line('x', 'y', source=source, line_width=3, line_alpha=0.6)

time_slider = Slider(start=0, end=10, value=1, step=.01, title="time")
# freq_slider = Slider(start=0.1, end=10, value=1, step=.1, title="Frequency")
# phase_slider = Slider(start=0, end=6.4, value=0, step=.1, title="Phase")
# offset_slider = Slider(start=-5, end=5, value=0, step=.1, title="Offset")

callback = CustomJS(args=dict(source=source, time=time_slider),
                    code="""
    const data = source.data;
    const A = time.value;
    const x = data['x']
    const y = data['y']
    # const xlst = source.data.x = [];
    # const ylst = source.data.y = [];
    for (let i = 0; i < x.length; i++) {
        y[i] +=A; 
        # xlst.push(x);
        # ylst.push(y);
    }
    source.change.emit();
""")

time_slider.js_on_change('value', callback)
layout = row(
    plot,
    column(time_slider),
)

show(layout)


# source = ColumnDataSource(data=dict(x=soln.t, y=soln.y[0]))
# # print(source)
# plot = Figure(plot_width = 400, plot_height = 400)
# plot.line('x', 'y', source = source, line_width=3, line_alpha = 0.6)
# # show(plot)
# callback = CustomJS(args=dict(source= source), code="""
#     var data = source.data;
#     const N = 100, t_0 = 0, t_1 = 1, y_0 = 2
#     const h = (t_1 - t_0) / N  //time step size
#     var ts = Array.from(Array(N+1), (_, k) => k * h + t_0)
#     var ys = Array(N+1).fill(0)  //empty array for the results
#     ys[0] = y_0  //initial conditions

#     for (let i = 0; i < N; i++) {
#         ys[i + 1] =  ys[i] + f(ts[i], ys[i]) * h
#     }   
#     source.change.emit();
#     console.log(y)
#     """)

# slider = Slider(start=0, end = 4, value = 1, step=0.01, title = "ode_t")
# slider.js_on_change('value', callback)

# layout = column(slider, plot)
# show(layout)