# ui.R based on RStudio's stockVis app, Shiny lesson 6, here: http://shiny.rstudio.com/tutorial/lesson6/ 

library(shiny)
library(dygraphs)

shinyUI(fluidPage(
  titlePanel("Markov Switching Cointegrating Equation"),
  
  sidebarLayout(
    sidebarPanel(
		        helpText("Select a pair of commodities to examine. 
		        A Markov Switching Cointegrating equation will be estimated, and the filtered state 1 probability will be plotted to the right."),
    
      selectInput("symb1", label = h3("Select LHS Commodity"), 
		 choices = list("Crude Oil (West Texas)" = 1 ,
			       	"Crude Oil (Brent)" = 2,
			       	"Natural Gas" = 3,
			       	"Heating Oil (No. 2, NY Harbor)" = 4,
			       	"Gulf Coast Gasoline" = 5,
			       	"RBOB Gasoline (LA)" = 6,
			       	"Jet Fuel (Kerosene-type Gulf Coast)" = 7,
			       	"Diesel Fuel (Gulf Coast, Ultra-Low Sulfur No.2)" = 8),
		selected = 3),  
    
      selectInput("symb2", label = h3("Select RHS Commodity"), 
		 choices = list("Crude Oil (West Texas)" = 1,
			       	"Crude Oil (Brent)" = 2,
			       	"Natural Gas" = 3,
			       	"Heating Oil (No. 2, NY Harbor)" = 4,
			       	"Gulf Coast Gasoline" = 5,
			       	"RBOB Gasoline (LA)" = 6,
			       	"Jet Fuel (Kerosene-type Gulf Coast)" = 7,
			       	"Diesel Fuel (Gulf Coast, Ultra-Low Sulfur No.2)" = 8),
		selected = 1),  

      dateRangeInput("dates", 
        "Date range",
        start = "1997-01-01", 
        end = as.character(Sys.Date())),
   
br(),
br(),

sliderInput("init", "Initial State 1 Probability", min = 0, max = 1, value = 0.5, step= 0.05)

  ),
    
    mainPanel(dygraphOutput("plot"))
  )
))

