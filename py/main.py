from bokeh.plotting import figure, output_file, save, show, ColumnData
# columndatascource allows tooltip and other stuffs
## show: to show the html file
from bokeh.models.tools import HoverTool
from bokeh.transform import factor_cmap
from bokeh.palettes import Blues8
from bokeh.embed import components # integrate into another file

import pandas
# read in csv
df = pandas.read_csv('cars.csv')

#Create ColumDataSource from data frame
source = ColumnDataSource(df)

output_file('index.html')

# car list
car_list = source.data['Car'].tolist()

# Add plot
p = figure(
    y_range = car_list,
    plot_width = 800,
    plot_height = 600, #pixel
    title = 'cars with top horsepower',
    x_axis_label = 'horsepower',
    tools = "pan, box_select, zoom_in, zoom_out, save, reset" # no tools at all
)

# render glyph

p.hbar(
    y = 'Car',
    right = 'Horsepower',
    left = 0,
    height = 0.4,
    fill_color = factor_cmap(
        'Car',
        palette=Blues8,
        factors = car_list
    ),
    fill_alpha = 0.9,
    source = source,
    legend='Car'
)

# Add Legend
p.legend.orientation = 'vertical'
p.legend.location = 'top_right'
p.legend.label_text_font_size = '10px'

# Add tooltips
hover = HoverTool()
hover.tooltips = """
    <div>
    <h3>@Car</h3>
    <div> <strong>Price: </stronh>@Price</div>
    <div> <strong>HP: </strong>@Horsepower</div>
    <div> <img src = "@Image" alt = " " width = "200"/> </div>
    </div>
"""
p.add_tools(hover)
# show result
# show(p)
save(p)

# print out div and script
# script, div = components(p)
# print(div)
# print(script)