from scipy.integrate import solve_ivp
from bokeh.layouts import column
from bokeh.models import ColumnDataSource, CustomJS, Slider
from bokeh.plotting import Figure, output_file, show
import matplotlib.pyplot as plt
import math


def deriv(t, y):
    x, y, z = y
    xdot = math.sin(x)
    ydot = 0.04 * x - 1.e4 * y * z - 3.e7 * y**2
    zdot = 3.e7 * y**2
    return xdot, ydot, zdot

t0, tf = 0, 10
y0 = 2, 0, 0
soln = solve_ivp(deriv, (t0, tf), y0)
print(soln.t)
print(soln.y[0])

# for i in range(soln.y.shape[0]):
#     plt.plot(soln.t, soln.y[i], label=f'$X_{i}(t)$')
# plt.xlabel('$t$') # the horizontal axis represents the time 
# plt.legend() # show how the colors correspond to the components of X
# plt.show()


source = ColumnDataSource(data=dict(x=soln.t, y=soln.y[0]))
# print(source)
plot = Figure(plot_width = 400, plot_height = 400, y_range=(1, 3))
plot.line('x', 'y', source = source, line_width=3, line_alpha = 0.6)
# show(plot)
ode_slider = Slider(start=0, end = 1, value = 1, step=0.001, title = "ode_t")
callback = CustomJS(args=dict(source= source, ode = ode_slider), code="""
    const data = source.data;
    var f = cb_obj.value
    var x = data['x']
    var y = data['y']
    const A = ode.value;
    for(var i =0; i<x.length; i++){
        y[i] = y[i] + A
    }
    source.change.emit();
    console.log(y)
    """)


ode_slider.js_on_change('value', callback)

layout = column(ode_slider, plot)
show(layout)