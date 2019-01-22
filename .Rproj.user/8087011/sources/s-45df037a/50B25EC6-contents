library(dplyr)
# library(metacell)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(shinydashboard)
library(data.table)
library(DT)
library(colourpicker)
library(shinyjs)

function(input, output, session) {
  
  local_mc_colorize <- function(mc_fp, marker_colors, sequential_coloring = FALSE) 
  {
    if (class(marker_colors)[1] != "data.frame" | length(intersect(c("gene", 
                                                                     "group", "color", "priority", "T_fold"), colnames(marker_colors))) != 
        5) {
      stop("MC-ERR marker colors parameter must be a data frame with fields gene, group, color, priority, T_fold")
    }
    marker_colors$gene = as.character(marker_colors$gene)
    marker_colors$color = as.character(marker_colors$color)
    rownames(marker_colors) = marker_colors$gene
    good_marks = intersect(rownames(marker_colors), rownames(mc_fp))
    if (length(good_marks) == 0) {
      message("no color markers are found")
      return
    }
    marker_colors = marker_colors[good_marks, ]
    color_key = as.data.frame(marker_colors)
    cl_colors = rep(NA, ncol(mc_fp))
    if (sequential_coloring) {
      for (p in sort(unique(marker_colors$priority))) {
        curr_marker_colors = marker_colors[marker_colors$priority == p, ]
        marker_fold = mc_fp[curr_marker_colors$gene,]
        marker_fold = ifelse(marker_fold > curr_marker_colors$T_fold, Marker_fold, NA)
        if (nrow(curr_marker_colors) == 1) {
          passed = is.na(cl_colors) & !is.na(marker_fold)
          hit = rep(1, sum(passed))
        }
        else {
          passed = is.na(cl_colors) & colSums(!is.na(marker_fold)) > 
            0
          hit = apply(marker_fold[, passed], 2, which.max)
        }
        cl_colors[passed] = curr_marker_colors[hit, "color"]
      }
    }
    else {
      marker_colors = marker_colors[order(marker_colors$priority),]
      marker_fold = mc_fp[marker_colors$gene, ]
      marker_fold = ifelse(marker_fold > marker_colors$T_fold,log2(marker_fold), NA)
      marker_fold = marker_fold * marker_colors$priority
      if (length(good_marks) > 1) {
        nonz = colSums(!is.na(marker_fold)) > 0
        hit = apply(marker_fold[, nonz], 2, which.max)
      }
      else {
        nonz = marker_fold > 0
        hit = rep(1, sum(nonz))
      }
      cl_colors[nonz] = marker_colors[hit, "color"]
    }
    cl_colors[is.na(cl_colors)] = "gray"
    names(cl_colors) = colnames(mc_fp)
    return(cl_colors)
  }
  
  
  # Creating reactive values (from default_mc_col table) since we will need to modify the table with events and triggers
  vals=reactiveValues()
  vals$Data <- data.table(default_mc_col)
  name2color=reactiveValues()
  name2color$Data <- setNames(unique(default_mc_col[,c("group","color")])$color,unique(default_mc_col[,c("group","color")])$group)
  color2name = reactiveValues()
  color2name$Data <- setNames(unique(default_mc_col[,c("group","color")])$group,unique(default_mc_col[,c("group","color")])$color)
  
##### XY scatter plot (gene1 vs gene2) ##### 
  # updating a data frame to use it
  select_fp <- reactive({
    log_trans <- input$log_trans
    if(log_trans){
      fp = log2(mc_fp)
    }
    else{
      fp = mc_fp
    }
    fp = as.data.frame(t(fp))
    updateNumericInput(session, inputId ="t_gb",label = NULL,value = input$t_gb,min = floor(min(fp[,input$yvar])),max = ceiling(max(fp[,input$yvar])),step = 0.5)
    updateNumericInput(session, inputId ="t_ga",label = NULL,value = input$t_ga,min = floor(min(fp[,input$xvar])),max = ceiling(max(fp[,input$xvar])),step = 0.5)
    return(fp)
  })
  
  # functions for ablines
  hline <- function(y = 0,x0 = 0, x1= 1, color = "black") {
    list(
      type = "line", 
      x0 = x0, 
      x1 = x1, 
      xref = "paper",
      y0 = y, 
      y1 = y, 
      line = list(color = color)
    )
  }
  
  vline <- function(x = 0, y0 = 0, y1=1, color = "black") {
      list(
        type = "line", 
        y0 = y0, 
        y1 = y1, 
        yref = "paper",
        x0 = x, 
        x1 = x, 
        line = list(color = color)
      )
    }
  
  # A reactive expression with the plotly - XY plot
  output$ga_gb_plot <- renderPlotly({
    xvar_name <- input$xvar
    yvar_name <- input$yvar
    fp = select_fp()
    fp$cell_type = factor(as.character(color2name$Data[mc_colors[rownames(fp)]]),levels = names(name2color$Data))
    fp$color = factor(as.character(mc_colors[rownames(fp)]),levels = names(color2name$Data))
    fp$mc_num = as.character(1:nrow(fp))
    vert <- as.double(input$t_ga)
    horiz <- as.double(input$t_gb)
    p <- plot_ly(fp[,c(xvar_name,yvar_name,"cell_type","mc_num")],source = "ga_gb", type = "scatter",mode = "markers",
                 x = ~fp[[xvar_name]], y = ~fp[[yvar_name]],color = ~cell_type,colors = levels(fp$color),
                 marker = list(size = 20,opacity = 0.5), key = ~mc_num,
                 hoverinfo = "text", text = ~paste("<b>","mc: ",fp[["mc_num"]] ,"</b><br><b>","cell type: ", fp[["cell_type"]] ,"</b><br>"))%>%
      layout(xaxis = list(title = xvar_name), yaxis = list(title = yvar_name) ,
             shapes= list(hline(y=horiz,x0 = floor(min(fp[[xvar_name]])), x1 = ceiling(max(fp[[xvar_name]]))),vline(x=vert,y0 =floor(min(fp[[yvar_name]])), y1 =ceiling(max(fp[[yvar_name]]))+1 ))) %>%
      add_annotations(text = ~mc_num,size = 0.3,color = I("black"),  showarrow = FALSE,ax = 0,ay =0) %>% hide_legend() 
    p
  })
  
  
  
##### 2d-projection of colored metacells ##### 
  output$mc_2d <- renderPlotly({
    df_2d_type =  cbind(df_2d,cell_type = color2name$Data[mc_colors[mc_mc[rownames(df_2d)]]])
    df_mc_2d_type = cbind(df_mc_2d, cell_type = color2name$Data[mc_colors[df_mc_2d$mc_num]])
    g<-ggplot(df_2d_type,aes(x=sc_x,y=sc_y,color = cell_type,key=mc_num)) + geom_point(size=1) + geom_text(data=df_mc_2d_type,aes(x=mc_x,y=mc_y,label=mc_num),color = "black",size=2) +
      theme_void() + scale_color_manual(values=as.character(name2color$Data),limits=names(name2color$Data),breaks=names(name2color$Data),labels=names(name2color$Data),name="")
    ggplotly(g, tooltip=c("key","cell_type"),source = "mc2d") %>% 
      layout( xaxis = list(showgrid = F),yaxis = list(showgrid = F),legend = list(x=max(df_2d_type$sc_x)+1)) %>% hide_legend() 
    })
  
###### Annotating interactive table ######

  ##### CSS styles for buttons and inputs ##### 
  # tags$head(
  #   tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css")
  # )
  
  tags$head(
    tags$style(type="text/css",
               "*{margin:0;box-sizing:border-box;}
               html,body{height:100%;font:14px/1.4 sans-serif;}
               input, textarea{font:14px/1.4 sans-serif;}

               .input-group{
               display: table;
               border-collapse: collapse;
               width:100%;
               }
               .input-group > div{
               display: table-row;
               background:#eee;
               padding: 0px 0px;
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
               padding: 0px 0px;
               }

               .input-group-area{
               display: table-cell;
               border: 0px;
               background:#eee;
               }

               .select-input .selectize-control{
               display: inline-block;
               text-align: left;
               width: 100%;
               padding: 0px 0px;
               }
               .numeric-input .form-group{
               display: inline-block;
               text-align: left;
               width: 100%;
               padding: 0px 0px;
               }

              .text-input .input{
               display: inline-block;
               text-align: left;
               width: 100%;
               padding: 0px 0px;
              }

              .bttn-input button{
              display: inline-block;
              text-align: center;
              background-color: #73a5a5;
              border: 1px solid #92b9b9;
              padding: 50% 0%;
              color: white;
              cursor: pointer;
              width: 100%;
              height: 100%;
              float: left;}
               "))
  
file_list_choices <- sapply(list.files("db/color_schemes/"), function(x){paste0("<option value=",x,">")})

  ##### Output full UI of the table (including buttons and titles) ##### 
  output$main_col_body=renderUI({
    fluidRow(
    fluidRow(h3("Colorize table",'style'="text-align:Left; vertical-align:top;padding:1px 10px 1px 10px")),
    fluidRow(
      column(12,align = "center",
             fluidRow(div(class='col-lg-2',style='float:left;',downloadButton('downloadData', 'Download')),
             HTML(
                 paste0("
                    <div class='col-lg-10'>
                      <div class='row'>
                          <div class='col-lg-4'>
                            <span class='input-group-addon' style='width:100%;float:left;padding: 10px 4px 10px 4px;'>Upload from computer</span>
                          </div>
                          <div class='col-lg-6'>
                            <div class='form-group shiny-input-container'>
                              <div class='input-group'>
                                <label class='input-group-btn'>
                                  <span class='btn btn-default btn-file'>
                                    Browse...
                                  <input id='upload_mc_col' name='upload_mc_col' type='file' style='display: none;' accept='.txt'/>
                                  </span>
                                </label>
                                <input type='text' class='form-control' placeholder='No file selected' readonly='readonly'/>
                              </div>
                          </div>
                        </div>
                        <div class='col-lg-2'>
                          <button class='btn btn-outline-secondary' style='float:right;' type='button' id = 'upload_mc_col_file'>Upload</button>
                        </div> 
                      </div>
                      <div class='row'>
                          <div class='col-lg-4'>
                            <span class='input-group-addon' style='width:100%;float:left;padding: 10px 4px 10px 4px;'>Upload from our db</span>
                          </div>
                          <div class='col-lg-6'>
                        	    <input list='file_list' id='chosen_mc_col' value=", names(file_list_choices)[1] ," style='width:100%;float:right;border:0px;padding: 4px 0px 4px 0px;'>
                              <datalist id='file_list'>",paste0(file_list_choices,collapse=""),"</datalist>
                          </div>
                          <div class='col-lg-2'>
                            <button class='btn btn-outline-secondary' style='float:right;' type='button' id = 'load_mc_col_file'>Upload</button>
                          </div>
                        </div>
                    </div>"))),
              #<label class='btn btn-outline-secondary' style='float:left;width:100%;border:solid 1px #eee;padding: 6px 6px 6px 6px;'>
             #<b>Browse</b>",
             # <span class='input-group-text' style='width:100%;float:right;padding: 10px 4px 10px 4px;' id='uploaded_mc_col'> No file has been chosen </span>
               
                 tags$datalist(names(file_list_choices),id="file_list"))
      ),
    fluidRow(
      column(12,           
             ##Rendering the datatable
             dataTableOutput("colors_table")
             )),
    fluidRow(
      tags$head(
        tags$style(type="text/css",
                   ".btn-styled button{background-color: #73a5a5; border: 1px solid #92b9b9; color: white; padding: 10px 24px; cursor: pointer;}
                .btn-styled:after {content: '';clear: both;display: table;}
                .btn-styled button:not(:last-child) {border-right: none;}
                .btn-styled button:hover {background-color: #466d6d;}"
        )
      ),
      column(3,
        div(class="btn-group-vertical",
          div(class="btn-styled",style="width:100%",
              actionButton(inputId = "apply_colorize",label = "Apply table colors",width = '100%'),
              actionButton(inputId = "random_colorize",label = "Apply random coloring",width = '100%'),
              actionButton(inputId = "reset_colorize",label = "Reset to default colors",width = '100%')
              )
          )
      ),

      column(3,offset = 4,
             div(class="btn-group-vertical",
                 div(class="btn-styled",style="width:100%",
                  actionButton(inputId = "Del_row_head",label = "Delete Selected rows",width = '100%'),
                  actionButton(inputId = "add_row",label = "Add new row",width = '100%')
                  )
              )
      )
      ),
    tags$script(
          HTML('$(document).on("click", "input", function () {
                var checkboxes = document.getElementsByName("row_selected");
                var checkboxesChecked = [];
                for (var i=0; i<checkboxes.length; i++) {
                if (checkboxes[i].checked) {
                  checkboxesChecked.push(checkboxes[i].value);
                }
              }
              Shiny.onInputChange("checked_rows",checkboxesChecked);})')),
    tags$script("$(document).on('click', '#colors_table button', function () {
        Shiny.onInputChange('lastClickId',this.id);
        Shiny.onInputChange('lastClick', Math.random())
                  });"),
    tags$script("$(document).on('click', '#add_row', function () {
        Shiny.onInputChange('lastClickId',this.id);
        Shiny.onInputChange('lastClick', Math.random())
                  });"),
    tags$script("$(document).on('click', '#load_mc_col_file', function () {
        var chosen_mc_col = document.getElementById('chosen_mc_col');
        Shiny.onInputChange('mc_file_to_load',chosen_mc_col.value);
        Shiny.onInputChange('lastClickId',this.id);
        Shiny.onInputChange('lastClick', Math.random())
                  });"),
    tags$script("$(document).on('click', '#upload_mc_col_file', function () {
      Shiny.onInputChange('lastClickId',this.id);
      Shiny.onInputChange('lastClick', Math.random())
                  });"),
    # tags$script("$(document).on('change','#upload_mc_col',function(){
    #   //get the file name
    #   var fileName = $(this).val();
    #   Shiny.onInputChange('mc_file_to_upload', fileName)
    #   //replace the 'Choose a file' label
    #   span = document.getElementById('uploaded_mc_col');
    #   span.innerHTML = fileName;
    # })"),
    tags$script(HTML("$(document).on('click', '#save_changes', function () {
        var list_value=[]
        for (i = 0; i < $( '.new_input' ).length; i++)
        {
          list_value.push($( '.new_input' )[i].value)
        }
        Shiny.onInputChange('newValue', list_value)
        });"))
      )})

output$downloadData <- downloadHandler(
  
  # This function returns a string which tells the client
  # browser what name to use when saving the file.
  filename = function() {
    paste0("mc_colorized_",strsplit(as.character(Sys.time())," ")[[1]][1],".txt")
  },
  
  # the argument 'file'.
  content = function(file) {
    # Write to a file specified by the 'file' argument
    write.table(vals$Data, file,sep = "\t",row.names = FALSE,quote = FALSE)
  }
)
  ##### render data table
  output$colors_table=renderDataTable({
    DT=vals$Data
    if(nrow(vals$Data) >0 ){
      DT[["color"]] <- paste0('<span style="background-color:',DT[["color"]],'">',DT[["color"]],'</span>')
      DT[["Select"]]<-paste0('<input type="checkbox" name="row_selected" value="Row',1:nrow(vals$Data),'">')
      DT[["Actions"]]<-paste0('<div class="btn-group" role="group" aria-label="delete_modify">
                                <button type="button" class="btn btn-secondary delete" style="padding:5px 5px 5px 5px;font-size=18px;" id=delete_',1:nrow(vals$Data),'>Delete</button>
                                <button type="button" class="btn btn-secondary modify" style="padding:5px 5px 5px 5px;font-size=18px;" id=modify_',1:nrow(vals$Data),'>Modify</button></div>
                              ')
      }
        else{
          DT[["color"]] <- c()
          DT[["Select"]] <- c()
          DT[["Actions"]]<- c()
        }
      datatable(DT,options = list(scrollCollapse = TRUE,paging = FALSE,scrollY = 500,scrollX = TRUE),escape=F)
      })
    
  observeEvent(input$apply_colorize,{
    color_key = as.data.frame(vals$Data)
    mc_colors <<- local_mc_colorize(mc_fp,marker_colors = color_key)
    color2name$Data <- setNames(unique(color_key[,c("group","color")])$group,unique(color_key[,c("group","color")])$color)
    name2color$Data <- setNames(unique(color_key[,c("group","color")])$color,unique(color_key[,c("group","color")])$group)
  })
  
  observeEvent(input$random_colorize,{
    mc_colors <<- colorRampPalette(c("darkgray", "burlywood1", "chocolate4","orange", "red", "purple", "blue","darkgoldenrod3", "cyan"))(ncol(mc_fp))
    mc_colors <<- setNames(mc_colors, as.character(1:length(mc_colors)))
    color2name$Data <- setNames(as.character(1:length(mc_colors)),mc_colors)
    name2color$Data <- setNames(mc_colors,as.character(1:length(mc_colors))) 
  })
  

  observeEvent(input$reset_colorize,{
    mc_colors <<- setNames(read.delim("db/mc_colors.txt",stringsAsFactors = FALSE)[[1]],colnames(mc_fp))
    vals$Data = data.table(default_mc_col)
    color2name$Data <- setNames(unique(default_mc_col[,c("group","color")])$group,unique(default_mc_col[,c("group","color")])$color)
    name2color$Data <- setNames(unique(default_mc_col[,c("group","color")])$color,unique(default_mc_col[,c("group","color")])$group)
  })
  
  ##### modals for adding / deleting / modifing rows   ##### 
  
  observeEvent(input$newValue,
               {
                 newValue=lapply(input$newValue, function(col) {
                   if (suppressWarnings(all(!is.na(as.numeric(as.character(col)))))) {
                     as.numeric(as.character(col))
                   } else {
                     col
                   }
                 })
                 DF=data.frame(lapply(newValue, function(x) t(data.frame(x))))
                 colnames(DF)=colnames(vals$Data)
                 if(input$lastClickId=="add_row"){
                   vals$Data = rbind(DF,vals$Data)
                 }
                 else{
                   vals$Data[as.numeric(gsub("modify_","",input$lastClickId))]<-DF
                 }
               }
               )
  

    output$row_addition<-renderDataTable({
      selected_row=nrow(vals$Data)+1
      old_row=vals$Data[selected_row]
      new_row=list()
      for (i in colnames(vals$Data[1]))
      {
        if (i == "color"){
          new_row[[i]]<-paste0('<input class="new_input" type="color" id=new_',i,'>')
        }
        else if (is.numeric(vals$Data[[i]]))
        {
          new_row[[i]]<-paste0('<input class="new_input" style="width:80px" type="number" id=new_',i,'>')
        }
        else
          new_row[[i]]<-paste0('<input class="new_input" type="text" id=new_',i,'>')
      }
      gene_choices <- sapply(mc_genes, function(x){   # turn choices into html
        paste0("<option value='",x,"'>")
      })
      new_row[["gene"]] = paste0('<input class="new_input" list="gene_list" id=new_gene value=',old_row[["gene"]],'>',
                                    '<datalist id="gene_list">',paste0(gene_choices,collapse=""),'</datalist>')
      tags$datalist(mc_genes,id="gene_list")
      new_row=as.data.table(new_row)
      setnames(new_row,colnames(vals$Data[1]))
      new_row
    },escape=F,options=list(dom='t',ordering=F,autoWidth = TRUE,  columnDefs = list(list(width = '70px', targets = c(4,5)))),selection="none")
    
    
    modal_add_row=modalDialog(
      fluidPage(
        h3(strong("Add new row"),align="center"),
        hr(),
        dataTableOutput('row_addition'),
        actionButton(inputId = "save_changes",label = "Save changes")
      ),
      size="l",easyClose = TRUE
    )
    observeEvent(input$save_changes,{
      removeModal()      
    })
    
    observeEvent(input$Del_row_head,{
      row_to_del=as.numeric(gsub("Row","",input$checked_rows))
      vals$Data=vals$Data[-row_to_del]}
    )
    

    output$row_modif<-renderDataTable({
      selected_row=as.numeric(gsub("modify_","",input$lastClickId))
      old_row=vals$Data[selected_row]
      row_change=list()
      for (i in colnames(old_row))
      {
        if (is.numeric(vals$Data[[i]]))
        {
          row_change[[i]]<-paste0('<input class="new_input" style="width:80px" type="number" id=new_',i,' value=',old_row[[i]],'>')
        }
        else{
          row_change[[i]]<-paste0('<input class="new_input" type="text" id=new_',i,' value=',old_row[[i]],'>')
        }
      }
      gene_choices <- sapply(mc_genes, function(x){   # turn choices into html
        paste0("<option value='",x,"'>")
      })
      tags$datalist(mc_genes,id="gene_list")
      row_change[["gene"]] = paste0('<input class="new_input" list="gene_list" id=new_gene value=',old_row[["gene"]],'>',
                                    '<datalist id="gene_list">',paste0(gene_choices,collapse=""),'</datalist>')
      row_change[["color"]] = paste0('<input class="new_input" type="color" id=new_color value=',old_row[["color"]],'>')
      old_row[["color"]] <- paste0('<span style="background-color:',old_row[["color"]],'">',old_row[["color"]],'</span>')
      row_change=as.data.table(row_change)
      setnames(row_change,colnames(old_row))
      DT=rbind(old_row,row_change)
      rownames(DT)<-c("Current values","New values")
      DT
    },escape=F,options=list(dom='t',ordering=F, autoWidth = TRUE,columnDefs = list(list(width = '70px', targets = c(4,5)))),selection="none")
    
    modal_modify=modalDialog(
      fluidPage(
        h3(strong("Row modification"),align="center"),
        hr(),
        dataTableOutput('row_modif'),
        actionButton(inputId = "save_changes",label = "Save changes")
      ),
      size="l",easyClose = TRUE
    )
    
    observeEvent(input$lastClick,
                 {
                   if (input$lastClickId%like%"delete")
                   {
                     row_to_del=as.numeric(gsub("delete_","",input$lastClickId))
                     vals$Data=vals$Data[-row_to_del]
                   }
                   else if (input$lastClickId%like%"modify")
                   {
                     showModal(modal_modify)
                   }
                   else if (input$lastClickId == "add_row"){
                     showModal(modal_add_row)
                   }
                   else if (input$lastClickId == "upload_mc_col_file"){
                     tryCatch(
                       {
                         new_vals <- read.delim(input$upload_mc_col$datapath)
                       },
                       error = function(e) {
                         # return a safeError if a parsing error occurs
                         stop(safeError(e))
                       }
                     )
                      if(!is.na(new_vals) && ncol(new_vals) == 5 && colnames(new_vals) == c("group","gene","color","prioriy","T_fold")){
                         vals$Data=data.table(new_vals)
                       }
                  
                   } 
                   else if (input$lastClickId == "load_mc_col_file"){
                     if(!is.na(input$mc_file_to_load)){
                       if(length(input$mc_file_to_load) >0 && input$mc_file_to_load %in% list.files("db/color_schemes/")){
                        fn = paste0("db/color_schemes/",input$mc_file_to_load)
                        vals$Data=data.table(read.delim(fn))
                       }
                     }
                     }
                   }
    )
    
    ##### 
    
    ##### 
    output$heat_marks <- renderPlotly({
      all_genes = union(marks_genes,as.character(vals$Data[["gene"]]))
      mc_ord = 1:ncol(mc_fp)
      cell_ord = names(mc_mc)[order(order(mc_ord)[mc_mc])]
      markers_mat = log2(mc_fp[all_genes, mc_ord])
      markers_mat = pmax(pmin(markers_mat,3),-3)
      df_text = melt(log2(mc_fp[all_genes, mc_ord]))
      colnames(df_text) = c("gene","mc","log2(mc_fp)")
      conditions.text <- paste0(paste0("mc: ",df_text$mc),"<br>",paste0("gene: ",df_text$gene),"<br>",paste0("log2 mc_fp: ", round(df_text$`log2(mc_fp)`,3)),"<br>")
      conditions.text <- matrix(unlist(conditions.text), ncol = length(mc_ord), byrow = FALSE)
      
      ######## axes ########
      x <- list(
        title = "mc",
        type = "category",
        ticks = "outside",
        nticks = length(mc_ord),
        tickfont = list(size = 7)
        )
      
      y <- list(
        title = "gene",
        tickfont = list(size = 8)
      )
      
      ax <- list(
        title = "",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE,
        ticks = ""
      )
      ####### A reactive expression with the plotly - heatmap of markers #########
      marks_heatmap <- plot_ly(
        x = mc_ord, y = all_genes,z=markers_mat, type = "heatmap", hoverinfo = 'text', 
        colors = colorRamp(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")),
        text = conditions.text
      ) %>%        layout(xaxis = ax, yaxis = y) 
      colors_matrix = matrix(data = mc_ord, nrow = 1, ncol = length(mc_ord))
      
      mc_colors_heatmap <- plot_ly(type="heatmap",x = mc_ord,z=colors_matrix,source ="heatmap_marks", hoverinfo = 'text',opacity=1,xgap = 0.5,
                                   text = matrix(paste0("mc:",mc_ord,"<br>","group:",color2name$Data[mc_colors[mc_ord]]),nrow=1,ncol=length(mc_ord)), 
                                   colors = mc_colors[mc_ord],showscale = FALSE)%>% layout(xaxis = ax,yaxis=ax)

      subplot(marks_heatmap, mc_colors_heatmap,  nrows = 2, shareX = FALSE, shareY = FALSE,heights = c(0.95,0.05)) 
      
    })
    
    
    gene_rows_dot_plot <- reactive({
      if(input$markers_only == TRUE){
        genes = as.character(vals$Data[["gene"]])
      }
      else{
        genes = rownames(mc_fp)
      }
      return(genes)
      })
    
    observeEvent(input$markers_only,{
      updateNumericInput(session, inputId = "num_10_colorize", label = "current 10" ,value = 1,min = 1, max = length(gene_rows_dot_plot()))
    })
    
    observeEvent(input$choose_mc,{
      updateNumericInput(session, inputId = "num_10_colorize", label = "current 10" ,value = 1,min = 1, max = length(gene_rows_dot_plot()))
    })

    observeEvent(input$top10_colorize, {
      if(input$num_10_colorize > 1){
        updateNumericInput(session, inputId = "num_10_colorize", label = "current 10" ,value = 1,min = 1, max = length(gene_rows_dot_plot()))
      }
    })
    
    observeEvent(input$prev_10_colorize, {
      current_10 = input$num_10_colorize
      if(current_10-1 > 0){
        updateNumericInput(session, inputId = "num_10_colorize", label = "current 10" ,value = current_10-1,min = 1, max = length(gene_rows_dot_plot()))
      }else{
        updateNumericInput(session, inputId = "num_10_colorize", label = "current 10" ,value = current_10,min = 1, max = length(gene_rows_dot_plot()))
      }
    })
    
    observeEvent(input$next_10_colorize, {
      current_10 = input$num_10_colorize
      gene_list = gene_rows_dot_plot()
      if(current_10+1 < ceiling(length(gene_list)/10)+1){
        updateNumericInput(session, inputId = "num_10_colorize", label = "current 10", value = current_10+1,min = 1, max = length(gene_list))
      }else{
        updateNumericInput(session, inputId = "num_10_colorize", label = "current 10",value = current_10,min = 1, max = length(gene_list))
      }
    })
    

    ##### horizontal scatter plot of mc footprint score (x axis) for each MetaCell for each gene (y axis), using ggplotly ##### 
    output$dot_plot_genes <- renderPlotly({
      gene_list = gene_rows_dot_plot()
      event_click_1<- event_data("plotly_click", source = "heatmap_marks")
      if(!is.null(event_click_1)){
          metacell = event_click_1$x
          updateNumericInput(session,inputId = "choose_mc",label = "mc", value = metacell ,min = 1, max = ncol(mc_fp))
          js$resetClick1()
        }
      
      event_click_2<- event_data("plotly_click", source = "ga_gb")
      if(!is.null(event_click_2)){
          metacell = event_click_2$key
          updateNumericInput(session,inputId = "choose_mc",label = "mc", value = metacell ,min = 1, max = ncol(mc_fp))
          js$resetClick2()
      }
      metacell = input$choose_mc
      
      if(!is.na(as.integer(metacell))){
        metacell = as.integer(metacell)
      }
      
      current_10 = input$num_10_colorize - 1
      range_text_10 = ifelse(current_10 > 0,sprintf("the %s'th-%s'th ",10*current_10+1,10*current_10+10),sprintf("%s'st - %s'th",10*current_10+1,10*current_10+10))
      gene_text_10 = ifelse(input$markers_only," colorizing genes"," genes")
      output$range_10 = renderText({paste0("Footprint values over all mc's for the ", range_text_10,gene_text_10)})
      output$range_10_2 = renderText({paste0("Ordered by the footprint values of mc: ",metacell)}) 
      top_genes = names(sort(mc_fp[gene_list,metacell],decreasing = TRUE)[c(1:10)+10*current_10])
      df = melt(mc_fp[top_genes,])
      colnames(df) = c("gene","mc_num","mc_fp")
      df$cell_type = color2name$Data[mc_colors[df$mc_num]]
      places =  sapply(top_genes, function(x) sum(mc_fp[x,] < mc_fp[x,metacell])+1)
      genes_places = sapply(top_genes, function(x) paste0(x," (# ",formatC(places[x],format = "f",digits =0)," out of ",length(colors),")"))
      df$genes_places = factor(genes_places[df$gene],levels = rev(as.character(genes_places)))
      x_lab = paste0("Genes (",metacell," mc_fp rank compared to others)")
      g <- ggplot(df, aes(x=genes_places, y = mc_fp,key = mc_num)) +
        geom_jitter(width = 0.3, height = 0.5, show.legend = FALSE,
                    aes(size = mc_num==metacell ,color = mc_num==metacell,
                        fill = cell_type,shape = mc_num==metacell)) +
        scale_shape_manual(values = c(22,21),limits = c(TRUE,FALSE), breaks = c(TRUE,FALSE),labels=c(TRUE,FALSE),name ="Point shape")+
        scale_size_manual(values = c(4,2),limits = c(TRUE,FALSE), breaks = c(TRUE,FALSE),labels=c(TRUE,FALSE))+
        scale_fill_manual(values=as.character(name2color$Data),limits=names(name2color$Data),breaks=names(name2color$Data),labels=names(name2color$Data),name="Cell Type")+
        scale_color_manual(values = c("black","white"),limits = c(TRUE,FALSE), breaks = c(TRUE,FALSE))+
        coord_flip() + theme_bw() + xlab(x_lab) + ylab("mc footprint")
      ggplotly(g, tooltip=c("mc_fp","cell_type","mc_num"),source = "dot_plot_color_genes") %>% hide_legend() 
      
    })
    
    getPage<-function() {
      return(includeHTML("lung_guide.html"))
    }
    output$lung_guide<-renderUI({getPage()})

  
}

