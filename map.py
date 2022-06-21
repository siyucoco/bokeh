
# importing the required modules
from bokeh.plotting import gmap
from bokeh.models import GMapOptions
from bokeh.io import output_file, show
from bokeh.models import RadioGroup, CustomJS

  
# file to save the model
output_file("gfg.html")
  
# configuring the Google map
lat = 40.6259
lng = 75.3705
map_type = "roadmap"
zoom = 5
google_map_options = GMapOptions(lat = lat,
                                 lng = lng,
                                 map_type = map_type,
                                 zoom = zoom)

L = ["First", "Second", "Third"]

# the active parameter sets checks the selected value
# by default
radio_group = RadioGroup(labels=L, active=1)

radio_group.js_on_click(CustomJS(code="""
	console.log('radio_group: active=' + this.active, this.toString())
"""))

show(radio_group)

  
# generating the Google map
google_api_key = ""
title = "Bethlehem"
google_map = gmap(google_api_key,
                  google_map_options,
                  title = title)
  
# displaying the model
show(google_map)