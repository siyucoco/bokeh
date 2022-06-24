import numpy as np
from scipy.integrate import solve_ivp

from bokeh.io import curdoc
from bokeh.layouts import row, column, gridplot
from bokeh.models import ColumnDataSource, ColorBar, LinearColorMapper, Slider, Div, HoverTool, Grid, LinearAxis, Tabs, Panel, Button
from bokeh.plotting import figure, show
from bokeh.palettes import Blues8

def lotkavolterra(t, z, k1, k2, k3, k4):
    x, y = z
    return [k1+k2*x*y, k3-k4*x-y]

t_start = 0.0
t_end = 15.0
N = 12 # number of data points
vec_time = np.linspace(t_start, t_end, N) # vector for time
print(vec_time.size)
# Starting values of all parameters
k1_start = 1
k2_start = 1
k3_start = 1.0
k4_start = 1.0
params = (k1_start, k2_start, k3_start, k4_start)
# print(params)

# Starting concentration of A, B, C
vec_conc_t0 = np.zeros(2)
vec_conc_t0[0] = 1.0
specie_names = ['X', 'Y']
vbar_top = (vec_conc_t0[0], vec_conc_t0[1])
specie_colors = ['darkgray', 'mediumblue']

# print(params,)
soln = solve_ivp(lotkavolterra, [t_start, t_end], [vec_conc_t0[0], vec_conc_t0[1]], args=(k1_start, k2_start, k3_start, k4_start))
# print(soln.y)
# print(soln.y[0])
int_vec_X = soln.y[0]
int_vec_Y = soln.y[1]
source = ColumnDataSource(data=dict(vec_time=vec_time, int_vec_X=int_vec_X, int_vec_Y=int_vec_Y))
source_vbar = ColumnDataSource(data=dict(specie_names=specie_names, vbar_top=vbar_top, color=specie_colors))

# Set up plot for concentrations
TOOLTIPS = [("Time (s)","@vec_time"), ("X","@int_vec_X{0,0.000}"), ("Y","@int_vec_Y{0,0.000}")]
TOOLS = "pan,undo,redo,reset,save,wheel_zoom,box_zoom"
plot_conc = figure(plot_height=450, plot_width=550, tools=TOOLS, tooltips=TOOLTIPS,
              title="Sequential reactions involving X,Y", x_range=[t_start, t_end], y_range=[0, 15])
plot_conc.line('vec_time', 'int_vec_X', source=source, line_width=3, line_alpha=0.6, line_color="darkgray",
               legend_label="x Concentration")
plot_conc.line('vec_time', 'int_vec_Y', source=source, line_width=3, line_alpha=0.6, line_color="navy",
               legend_label="y Concentration")
# plot_conc.line('vec_time', 'int_vec_C', source=source, line_width=3, line_alpha=0.6, line_color="darkorange",
#                legend_label="C Concentration")
plot_conc.xaxis.axis_label = "Time (s)"
plot_conc.yaxis.axis_label = "Concentration"
plot_conc.legend.location = "top_left"
plot_conc.legend.click_policy="hide"
plot_conc.legend.background_fill_alpha = 0.5
plot_conc.grid.grid_line_color = "silver"

# Set up widgets
k1 = Slider(title="k1"+" (initial: "+str(k3_start)+")", value=k3_start, start=2.02, end=8.0, step=0.02)
k2 = Slider(title="k2"+" (initial: "+str(k4_start)+")", value=k4_start, start=0.02, end=2.0, step=0.02)
k3 = Slider(title="k3"+" (initial: "+str(k1_start)+")", value=k1_start, start=1, end=5, step=1)
k4 = Slider(title="k4"+" (initial: "+str(k4_start)+")", value=k2_start, start=1, end=5, step=1)

def update_data(attrname, old, new):

    # Get the current slider values
    k_AB_temp = k1.value
    k_BC_temp = k2.value
    O_AB_temp = k3.value
    O_BC_temp = k4.value

    # Generate the new curve
    vec_time = np.linspace(t_start, t_end, N)  # vector for time
    # params_temp = [O_AB_temp, O_BC_temp, k_AB_temp, k_BC_temp]
    soln = solve_ivp(lotkavolterra, [t_start, t_end], [vec_conc_t0[0], vec_conc_t0[1]], args=(O_AB_temp, O_BC_temp, k_AB_temp, k_BC_temp))
    int_vec_X = soln.y[0]
    int_vec_Y = soln.y[1]
    source.data =  dict(vec_time=vec_time, int_vec_X=int_vec_X, int_vec_Y=int_vec_Y)
    specie_names = ['X', 'Y']
    specie_colors = ['darkgray', 'mediumblue']
    source_vbar.data = dict(specie_names=specie_names, color=specie_colors)

for w in [k1, k2, k3, k4]:
    w.on_change('value', update_data)

inputs_reaction = column(k1, k2, k3, k4)
tab1 =Panel(child=row(inputs_reaction, plot_conc, column( height=450)), title="Desktop")
tab2 =Panel(child=column(inputs_reaction, plot_conc, column( height=475)), title="Mobile")
tabs = Tabs(tabs = [tab1, tab2])

curdoc().add_root(tabs)
curdoc().add_root(tab1)
curdoc().title = "Ode IMPLEMENT"