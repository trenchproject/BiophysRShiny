 # server.R

setwd("C:\\Users\\Buckley\\Google Drive\\BuckleySynch\\CAREER\\Rcode\\MESA\\biophys\\")
source("CalcTe_30May2015.R")
#setwd("C:\\Users\\Buckley\\Google Drive\\BuckleySynch\\CAREER\\Rcode\\MESA\\biophys\\")
#source("Legend.R")
  
shinyServer(function(input, output, session) {

  #Get input
  #input$shape #("Cylinder","Sphere","Cube"),
  #input$metric #"Body temperature", "Thermal stress", "Population growth rate")
 
  #------------------------------------
  # Calculate Operative Environmental Temperatures 


  setwd("C:\\Users\\Buckley\\Google Drive\\BuckleySynch\\CAREER\\Rcode\\MESA\\climate\\")
  
output$myMap = renderLeaflet({
    
   Te<- TeRaster( as.numeric(input$monthk), as.numeric(input$hourk), as.numeric(input$length), as.numeric(input$mass), as.numeric(input$abs))
   
   Te.rast= raster("Tegrid.grd")
  
  # Create the map
  map <- createLeafletMap(session, "map")
  
  # download and load data
  cols<-colorBin("RdBu", c(minValue(Te.rast), maxValue(Te.rast)), 50)
  breaks= seq(minValue(Te.rast), maxValue(Te.rast), length.out = 50 )
  
pal <- colorNumeric(palette = "YlGnBu", domain = values(Te.rast) )
  
  map=leaflet() %>% addTiles() %>% setView(-100, 35, 4) %>%
    addRasterImage(Te.rast, colors = pal, opacity = 0.8, project=FALSE) %>%
    addLegend(pal = pal, position = "bottomleft", values = values(Te.rast),
              title = "Temp", opacity = 1)
  
  map})
  
  
}) #end server
