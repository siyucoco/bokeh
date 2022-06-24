
from bokeh.plotting import figure, output_file, show
from bokeh.models import Panel, Tabs
import numpy as np
import math
 
 
fig1 = figure(plot_width=300, plot_height=300)
 
x = [1, 2, 3, 4, 5]
y = [5, 4, 3, 2, 1]
 
fig1.line(x, y, line_color='green')
tab1 = Panel(child=fig1, title="Tab 1")
 
fig2 = figure(plot_width=300, plot_height=300)
 
fig2.line(y, x, line_color='red')
tab2 = Panel(child=fig2, title="Tab 2")
 
all_tabs = Tabs(tabs=[tab1, tab2])
 
show(all_tabs)