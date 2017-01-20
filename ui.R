#devtools::install_github("jcheng5/rasterfaster")
#devtools::install_github("rstudio/leaflet@raster")
#devtools::install_github("rstudio/leaflet")

#https://gist.github.com/jcheng5/a1c48a18d0ea3623c0e6  #addRaster

require(shiny)
require(leaflet)
require(rgdal)
require(raster)
require(mosaic)

# ui.R

shinyUI(

navbarPage("MESA", id="nav",
           
           tabPanel("Interactive map",
                    fluidPage(
                      sidebarLayout(
                                 sidebarPanel(   
                                   textInput("length", label = h3("Length (mm)"), value = 10),
                                   textInput("mass", label = h3("Mass (g)"), value = 10),
                                   textInput("abs", label = h3("Surface solar absorptivity (proportion)"), value = 0.9),
                                   
                                   selectInput("shape", 
                                               label = "Choose a shape",
                                               choices = c("Cylinder","Sphere","Cube"),
                                               selected = "Cylinder"),
                                   
                                   
                                   selectInput("metric", 
                                               label = "Choose a metric",
                                               choices = c("Body temperature", "Thermal stress", "Population growth rate"),
                                               selected = "Body temperature"),
                                   textInput("monthk", label = h3("Month"), value = 1),
                                   textInput("hourk", label = h3("Hour"), value = 1),
                                   submitButton("Run")
                                   
                                   
                                   #  sliderInput("date", 
                                   #            label = "Date range:",
                                   #           min = 1960, max = 2100,  nvalue = 2010)
                                   ), #end sidebar panel
                        
                      mainPanel(
                        textOutput("focal"),
                        leafletOutput('myMap')  
                      ) #end main panel
                      ) #end sidebar layout
           ) #end fluid page
           ), #end tabPanel
           
           tabPanel("Thermal Images",
                    fluidPage(  )
           )
           ) #end navbarpage
)

