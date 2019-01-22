library(dplyr)
# library(metacell)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(plotly)
library(shinydashboard)
library(shinyWidgets)
library(shinyjs)

# dashboardPage(
#   dashboardHeader(disable = T),
#   dashboardSidebar(disable = T),
#   dashboardBody()
# )

actionLink <- function(inputId, ...) {
  tags$a(href='javascript:void',
         id=inputId,
         class='action-button',
         ...)
}

tags$head(
  tags$style(type="text/css",
             ".btn-group button{background-color: #73a5a5; border: 1px solid #92b9b9; color: white; padding: 10px 24px; cursor: pointer; float: left;}
              .btn-group:after {content: '';clear: both;display: table;}
              .btn-group button:not(:last-child) {border-right: none;}
              .btn-group button:hover {background-color: #466d6d;}"
  )
)

css <- "
  *{margin:0;box-sizing:border-box;}
  html,body{height:100%;font:14px/1.4 sans-serif;}
  input, textarea{font:14px/1.4 sans-serif;}
  .container {
  margin: 0;
  padding: 0px;
  }
  .gray {
  background: #eee;
  }
  #select-container .selectize-input {
  font-size: 20px;
  line-height: 20px;
  }
  #select-container .selectize-dropdown {
  font-size: 18px;
  line-height: 24px;
  }
  #select-container .selectize-dropdown-content {
  max-height: 150px;
  padding: 0;
  }
  #select-container .selectize-dropdown-content .option {
  border-bottom: 1px dotted #ccc;
  }

  .input-group{
  display: table;
  border-collapse: collapse;
  width:100%;
  }

  .input-group > div{
  display: table-row;
  background:#eee;
  border: 0px;
  vertical-align: top;
  float: left;
  }

  .input-group-text{
  display: table-cell;
  vertical-align: top;
  border: 0px;
  background:#eee;
  color: #777;
  padding: 0px 5px 0px 5px;
  }

  .input-group-area{
  display: table-cell;
  vertical-align: top;
  border: 0px;
  background:#eee;
  padding: 0px 5px 0px 5px;
  }

  .input-group input{
  display: inline-block;
  vertical-align: top;
  }

  .select-input .selectize-control{
  display: inline-block;
  text-align: left;
  }

  .numeric-input .form-group{
  display: inline-block;
  text-align: center;
  }

  .num10 label{ display: inline; text-align: center; align-content: center; color: #777;}
  .num10 .form-group{ display: inline-block;width = 15px; text-align: center;}
  .actionButton.prevnext{display: inline-block; text-align: center; white-space: nowrap;}
  "


wellPanel(titlePanel("Analysis Of Single Cell RNAseq Data with MetaCell package",windowTitle = "ScRNAseq Analysis with MetaCell"),
          fluidPage(
            navlistPanel(widths = c(2, 10),
              tabPanel("Annoatating MetaCells",
                       useShinyjs(),
                       extendShinyjs(text = "shinyjs.resetClick1 = function() { Shiny.onInputChange('.clientValue-plotly_click-heatmap_marks', 'null'); }",functions = c("resetClick1")),
                       extendShinyjs(text = "shinyjs.resetClick2 = function() { Shiny.onInputChange('.clientValue-plotly_click-ga_gb', 'null'); }",functions = c("resetClick2")),
                       column(5,
                              uiOutput("main_col_body"),
                              div(class='col-lg-12',plotlyOutput("mc_2d"))
                              ),
                       column(7,
                              div(style = "height:70%",
                                  tabsetPanel(
                                    tabPanel("Marker genes distribution",
                                            ###### main tab ######
                                             h3("Marker genes distribution"),
                                             fluidRow(div(class='col-lg-12',plotlyOutput("heat_marks",height = 600)))
                                              ######
                                              ),
                                    tabPanel("Two genes scatter plot",
                                             h3("Two genes scatter plot"),
                                             ###### main tab ######
                                               fluidRow(
                                                 tags$head(tags$style(type="text/css",css)),
                                                 span(class = "input-group-text",style="width:100px", tags$label("X-axis Gene",'for'="xvar")),
                                                 span(class = "input-group-area",span(class ="select-container",style="width:190px", selectInput("xvar", label =NULL,mc_genes, selected = "Acta2",multiple = FALSE, selectize = TRUE, width = "190px"))),
                                                 span(class = "input-group-text",style="width:100px", tags$label("Y-axis Gene",'for'="yvar")),
                                                 span(class = "input-group-area", span(class ="select-container",style="width:190px", selectInput("yvar", label =NULL,mc_genes, selected = "Tgfbi",multiple = FALSE, selectize = TRUE, width = "190px"))),
                                                 span(class = "input-group-text",style="width:80px;",tags$label("log2 scale",'for'="log_trans")),
                                                 span(class = "input-group-area",prettyCheckbox(inputId = "log_trans",label = NULL,status = "success",outline = TRUE,inline = TRUE,value = FALSE)),
                                                 span(class = "input-group-text",style="width:50px",tags$label("Vline",'for'="t_ga")),
                                                 span(class = "input-group-area",span(class = "numeric-input",style="width:60px",numericInput(inputId ="t_ga",label = NULL,value = 0,min = 0,max = 100,step = 0.5,width = "100%"))),
                                                 span(class = "input-group-text",style="width:50px", tags$label("Hline",'for'="t_gb")),
                                                 span(class = "input-group-area",span(class = "numeric-input",style="width:60px",numericInput(inputId ="t_gb",label = NULL,value = 0,min = 0,max = 100,step = 0.5,width = "100%")))
                                                   ),
                                               fluidRow(
                                                 div(class='col-lg-12',plotlyOutput("ga_gb_plot",height = 600))
                                             )
                                    )
                                ),
                              tags$head(tags$style(type="text/css",css)),
                              div(style = "height:30%",
                                  fluidRow(column(12,align="center",h3(textOutput("range_10"),textOutput("range_10_2")))),
                                  #fluidRow( column(2,align="left", fluidRow(
                                             # span(class = "input-group-text",style="width:100px;",tags$label("Markers only",'for'="markers_only")), 
                                             # span(class = "input-group-area",prettyCheckbox(inputId = "markers_only",label = NULL,status = "success",outline = TRUE,inline = TRUE,value = FALSE)))),
                                  fluidRow(  
                                  column(10,align="center",
                                           fluidRow(
                                             tags$head(tags$style(type="text/css",css)),
                                             tags$span(class="num10",style = "padding: 0px 20px 0px 0px;",prettyCheckbox(inputId = "markers_only",label = "Markers only",status = "success",outline = TRUE,inline = FALSE,value = FALSE)),
                                             tags$span(class="num10",numericInput(inputId = "choose_mc",label = "mc", value = 1 ,min = 1, max = ncol(mc_fp))),
                                             actionButton("top10_colorize", class = "prevnext",label = "Top 10"),
                                             #uiOutput("top10_colorize",inline = TRUE),
                                             actionButton("prev_10_colorize", class = "prevnext",label = HTML("<i class='glyphicon glyphicon-arrow-left'></i> Prev 10")),
                                             # uiOutput("prev_10_colorize",inline = TRUE),
                                             tags$span(class="num10", numericInput(inputId = "num_10_colorize",label = "current 10", value = 1,min = 1, max = ceiling(length(mc_genes)/10))),
                                             actionButton("next_10_colorize", class = "prevnext",label = HTML("Next 10 <i class='glyphicon glyphicon-arrow-right'></i>"))
                                             # uiOutput("next_10_colorize",inline = TRUE)
                                             )
                                           )
                                    ),
                                fluidRow(div(class='col-lg-11',plotlyOutput("dot_plot_genes")))
                              )
                                )
                                
                              )
                       )
              ,tabPanel("Lung Rmarkdown Guide",column(12,uiOutput("lung_guide")))
              )
            )
          )



