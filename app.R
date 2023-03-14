library(shiny)
library(tidyverse)
library(ggplot2)
library(magrittr)
library(survminer)
library(survival)
library(openxlsx)
#Lachesis

#create data frame




# Define UI for data upload app ----
ui <- fluidPage(
    
    tags$head(
        tags$style(HTML("hr {border-top: 1px solid #000000;}"))
    ),
    
    # App title ----
    titlePanel("Lachesis"),
    
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
        
        # Sidebar panel for inputs ----
        sidebarPanel(
            #img(src='worm.png', align = "right",height = 70, width = 200),
            
            # Input: Select a file ----
            selectInput("dataformat","Select input data format",
                        c("WormBot" = "wormbot",
                          "TidyVerse" = "tidyverseformat")),
            fileInput("file1", "Choose CSV File",
                      multiple = TRUE,
                      accept = c("text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")),
            p("Select which groupings to graph and analyze"),
            uiOutput("analysisSelect"),
            hr(),
            selectInput("isCumHz","Select graph type",
                        c("Survival" = "NULL",
                          "Cumulative events" = "event",
                          "Cumulative hazard" = "cumhaz")),
            
            
            #download graph button
            numericInput("graphX", "Graph size in inches (x-axis)",value= 10),
            numericInput("graphY", "Graph size in inches (y-axis)",value= 6),
            
            downloadButton("graph", label = "Download graph"),
            
            #download graph button
            downloadButton("stats", label = "Download stats"),
            
            selectInput("fileForm","Select image file type",
                        c("PNG"="png",
                          "TIFF"="tiff",
                          "PDF"="pdf",
                          "SVG"=".svg")),
            br(),
            #save filename
            textInput("savename","Enter Savename"),
            hr(),
            
            #save filename
            br(),            
            uiOutput("colorselect")
            
        ),
        
        # Main panel for displaying outputs ----
        mainPanel(
            
            # Output: Data file ----
            #tableOutput("contents")
            tabsetPanel(type="tabs",
                        tabPanel("Survival curve",
                                 plotOutput("plot",width = "700px", height = "500px"),
                                 plotOutput("coxplot",width = "700px", height = "500px"),
                                 numericInput("xmin", "X axis min:",value= 0),
                                 numericInput("xmax", "X axis max:",value= 30),
                                 textInput("xAxisLabel","X-axis label",value=""),
                                 textInput("yAxisLabel","Y-axis label",value = ""),
                                 selectInput("legPos","Select legend position",
                                             c("Top right"= "right",
                                               "Top left" = "left",
                                               "None" = "none")),
                                 checkboxInput("confInt","Display confidence intervals"),
                                 checkboxInput("medLine","Display median lines"),
                                 textInput("tickmarks","Tick mark intervals",value = "1")),
                        tabPanel("Summary and Statistics",
                                 h4("Summary"),
                                 tableOutput("summary"),
                                 h4("Pairwise p-values"),
                                 tableOutput("survfit"),
                                 br(),
                                 selectInput("corrmethod","Select p-value adjustment",
                                             c("None" = "none",
                                               "Bonferroni" = "bonferroni",
                                               "Holm" =  "holm",
                                               "Hommel" = "hommel",
                                               "Hochberg" = "hochberg",
                                               "Benjamini-Hochberg" = "BH",
                                               "Benjamini-Yekutieli" = "BY"))
                                 ),
                        tabPanel("Raw Data",
                                 downloadButton("tidyversedata", label = "Download TidyVerse Dataset"),
                                 tableOutput("contents")
                                 ),
                        tabPanel("Instructions",
                                 p("Lifespan plotter accepts comma delimited .csv files as input. Files can either be in WormBot format or in TidyVerse Format."),
                                 br(),
                                 p("WormBot format is shown below with # signs denoting conditional fields. 
                                    Afterwards, each row should have an age in the first column, then either
                                    a 1 for a death event or a 0 for a censored animal in the second column."),
                                 br(),
                                 tableOutput("exampleTable"), 
                                 br(),
                                 downloadButton("exampleFile", label = "Download example wormbot file"),
                                 br(),
                                 p("TidyVerse format is a .csv in which the first row is the column names and each animal is given unique row that includes a column for the age-at-death. AT LEAST
                                   one conditional column must be included with a conditional field such as 'strain' with entries for every row
                                   Please note the age-at-death column name must be formatted as 'age'. If some animals were lost from observation, include a column named'censored' 
                                   with either a '0' for a censored animal or a '1' for an observed death."))
            )
        )
    )
)

# Define server logic to read selected file ----
server <- function(input, output) {
    
    ####
    #creates and updates the table as necessary.
    create_table <- reactive({
        
        if(input$dataformat == 'wormbot'){
        dataf  <- read.csv(input$file1$datapath, header=FALSE, comment.char="%", stringsAsFactors=FALSE)
        breakpoint <- dataf[1,1]
        breaks <- grep(breakpoint,dataf$V1) #find # signs
        breaksfull <- breaks
        connames <-dataf[breaks,1]
        nWells <- length(breaks)
        breaks[(length(breaks)+1)] <- length(dataf$V1)+1 #find endpoint of document and add it to breaks       
        
        #for each well
        for (i in 1L:nWells) {
            SingleWell <- slice(dataf,breaks[i]:(breaks[i+1]-1)) #extract data
            NumRows <- length(SingleWell[,1])
            Cons <- grep("#",SingleWell$V1)
            ConTable <-slice(SingleWell,Cons[1]:length(Cons))
            Anustart <- length(Cons)
            colnames(SingleWell) <- c("age","censored") 
            SingleWell <- SingleWell[-1:-Anustart,]
            SingleWell <- as_tibble(SingleWell)
            SingleWell %>% mutate(age = as.numeric(age),censored = as.numeric(censored))
            #SingleWell <- unlist(mapply(rep, SingleWell$age, SingleWell$censored))
            SingleWell <- as_tibble(SingleWell)
            colnames(SingleWell)<- c("age","censored")
            #SingleWell$censored <- 1
            
            #for each field exported
            for (k in Cons){
                NewCol <- character()
                ConName <- str_remove(ConTable[k,1],"\\#")
                ConName <- gsub(" ", "_", ConName, fixed = TRUE)
                NewCol[1:(Anustart)] <- ConName
                NewCol[(Anustart+1):NumRows] <- ConTable[k,2]
                NewCol <- as_tibble(NewCol)
                colnames(NewCol)<-NewCol[[1,1]]
                NewCol <- NewCol[-1:-Anustart,]
                SingleWell <- add_column(SingleWell,NewCol)
            }
            
            #if first well then create totalTable variable    
            if (i == 1L) {
                totalTable <- SingleWell
            }
            #if subsequent wells append new data to totalTable    
            if (i > 1L) {
                totalTable <- bind_rows(totalTable,SingleWell)
            }  
        }
        #book-keeping to make sure that variables are kept as numeric despite data-wrangling
        totalTable <- mutate(totalTable, age = as.numeric(age),censored = as.numeric(censored))
        #print(totalTable)
        return(totalTable)
        }
        else{
            dataf  <- read.csv(input$file1$datapath, header=TRUE, stringsAsFactors=FALSE)
            if("rls" %in% colnames(dataf))
            {
                dataf<- dataf %>% rename(age = rls)
            }
            if("RLS" %in% colnames(dataf))
            {
                dataf<- dataf %>% rename(age = RLS)
            }
            if("censored" %in% colnames(dataf)==FALSE)
            {
                dataf$censored <- 1
            }
            dataf <- dataf %>% relocate(censored) %>% relocate(age)
            totalTable <- dataf %>% mutate(age = as.numeric(age), censored = as.numeric(censored))
            
        }
        
    })
    
    xmarks <- reactive({
        xint <- as.numeric(input$tickmarks)
        return(xint)
    })
    
    showcon <- reactive({
        req(input$file1)
        totalTableF <- create_table()
        conditionNames <- colnames(totalTableF)[3:length(colnames(totalTableF))]
        conditionBool <- character()
        varname <- character()
        for(i in 1:length(conditionNames)){
            
            varname[i] <- paste0('chk_',gsub(" ", "", conditionNames[i], fixed = TRUE))
            if(input[[varname[i]]]==TRUE){
                conditionBool[i] <- conditionNames[i]
            }
            else{
                conditionBool[i] <- "NULL"    
            }          
        }
        conditionSurv <- paste(conditionBool,collapse = " + ")
    })
    
    ###
    #creates analysis selector objects for each condition
    output$analysisSelect <- renderUI({
        req(input$file1)
        totalTableF <- create_table()
        conditionNames <- colnames(totalTableF)[3:length(colnames(totalTableF))]
        lapply(conditionNames, function(x){
            checkboxInput(paste0('chk_',gsub(" ", "", x, fixed = TRUE)),label = x)
        })
    }) 
    
    
    ###
    #creates color selector objects for each condition
    output$colorselect <- renderUI({
        req(input$file1)
        totalTableF <- create_table()
        fit3 <- eval(parse(text = paste("survfit(Surv(age, censored) ~ ", showcon(), ", data=totalTableF)")))
        #print(showcon())
        namesCh <- gsub("description=","",names(fit3$strata))
        namesCh <- gsub("description=","",namesCh)
        
        lapply(namesCh, function(x){
            fluidRow(
                textInput(paste0('labtxt_',gsub(" ", "", x, fixed = TRUE)),
                          label = NULL,
                          value = x),
                selectInput(paste0('lab_',gsub(" ", "", x, fixed = TRUE)), 
                            label = NULL, 
                            choices = c("black","gray","blue","green","orange","purple","red","cyan","yellow","goldenrod","darkgreen","violet","coral")),
                selectInput(paste0('line_',gsub(" ", "", x, fixed = TRUE)), 
                            label = NULL, 
                            choices = c("solid","dashed","dotted")),
                checkboxInput(paste0('omit_',gsub(" ", "", x, fixed = TRUE)), 
                            label = 'Remove Curve'),
                
                
                hr()
            )
        })
    })
    
    ###
    #show raw data
    output$contents <- renderTable({
        req(input$file1)
        totalTableF <- create_table()
        return(totalTableF)
        
    })
    
    ###
    #plot graph and update when necessary
    plotInput <-  reactive({
        req(input$file1)
        totalTableF <- create_table()
        fit3 <- eval(parse(text = paste("survfit(Surv(age, censored) ~ ", showcon(), ", data=totalTableF)")))
        #if (nrow(distinct(select(totalTableF,description))) > 1){
        namesCh <- gsub("description=","",names(fit3$strata))
        #}
        #else{
        #    namesT <- distinct(select(totalTableF,description))
        #    namesCh <- namesT$description
        #}
        
        cLabels <- character() 
        varname <- character()
        
        for(i in 1:length(namesCh)){
            varname[i] <- paste0('labtxt_',gsub(" ", "", namesCh[i], fixed = TRUE))
            cLabels[i] <- input[[varname[i]]]
        } 
        
        cPallete <- character()
        varname2 <- character()          
        for(i in 1:length(namesCh)){
            varname2[i] <- paste0('lab_',gsub(" ", "", namesCh[i], fixed = TRUE))
            cPallete[i] <- input[[varname2[i]]]
        } 
        cLine <- character()
        varname3 <- character()          
        for(i in 1:length(namesCh)){
            varname3[i] <- paste0('line_',gsub(" ", "", namesCh[i], fixed = TRUE))
            cLine[i] <- input[[varname3[i]]]
        } 
        
        cOpaque <- numeric()
        varname4 <- character()          
        for(i in 1:length(namesCh)){
            varname4[i] <- paste0('omit_',gsub(" ", "", namesCh[i], fixed = TRUE))
            if(input[[varname4[i]]] == TRUE){
                cOpaque[i] <- 0
            } else {
                cOpaque[i] <- 1
            }
        }
        
        #cLabels[cOpaque==0] <- "    "
        cPallete[cOpaque==0] <- "white"
        cLine[cOpaque==0]<-"blank"
        
        
        
        
        if(input$isCumHz=="event"){
        finalPlot <- ggsurvplot(fit3, data = totalTableF, conf.int = input$confInt,ylim = c(0,1.1),xlim=c(input$xmin,input$xmax),xlab= input$xAxisLabel,ylab = input$yAxisLabel,palette = alpha(cPallete,cOpaque),legend = legendPos(),legend.labs = cLabels, linetype = cLine, legend.title = "Legend",title= input$savename,break.x.by=xmarks(),break.y.by = 0.1,surv.median.line = displayMedian(),fun = input$isCumHz,censor = FALSE)
        }
        else if(input$isCumHz=="NULL"){
        finalPlot <- ggsurvplot(fit3, data = totalTableF, conf.int = input$confInt,ylim = c(0,1.1),xlim=c(input$xmin,input$xmax),xlab= input$xAxisLabel,ylab = input$yAxisLabel,palette = alpha(cPallete,cOpaque),legend = legendPos(),legend.labs = cLabels, linetype = cLine, legend.title = "Legend",title= input$savename,break.x.by=xmarks(),break.y.by = 0.1,surv.median.line = displayMedian(),censor = FALSE)
          
        }else{
        finalPlot <- ggsurvplot(fit3, data = totalTableF, conf.int = input$confInt,xlim=c(input$xmin,input$xmax),xlab= input$xAxisLabel,ylab = input$yAxisLabel,palette = alpha(cPallete,cOpaque),legend = legendPos(),legend.labs = cLabels, linetype = cLine, legend.title = "Legend",title= input$savename,break.x.by=xmarks(),surv.median.line = displayMedian(),fun = input$isCumHz,censor = FALSE)
          
        }
        return(finalPlot)
        
    })
    
   
    ###
    #displays graph in GUI
    output$plot <- renderPlot({
        print(plotInput())
    })
    
   
    ###
    #reactive input for median line control
    displayMedian <- reactive({
        if(input$medLine == TRUE){
            return("v")
        }
        else {
            return("none")
        }
    })
    
    ###
    #reactive function for legend position in graph
    legendPos <- reactive({
        if(input$legPos == "right"){
            val <- c(.8,.8)
        }
        else if(input$legPos=="left"){
            val <- c(.2,.8)
        }
        else {
            val <- "none"
        }
        
        return(val)
        
    })
    
    ###
    #creates statistical table
    statsInput <-  reactive({
        req(input$file1)
        totalTableF <- create_table()
        
        fitsurfdif <- eval(parse(text = paste("pairwise_survdiff(Surv(age, censored) ~ ", showcon(), ",data=totalTableF,p.adjust.method = input$corrmethod)")))
        
        #fitsurfdif <- pairwise_survdiff(Surv(age, censored) ~ description, data=totalTableF,p.adjust.method = input$corrmethod)
        fitpval <- fitsurfdif$p.value
        z <-1
        newTable <- tibble(Treatments=character(),pvalue=numeric())
        for(i in seq(1,ncol(fitpval))){
            for(k in seq(1,nrow(fitpval))){
                
                if(is.na(fitpval[k,i]) == FALSE){
                    newTable[z,1] <- paste0(colnames(fitpval)[i]," vs ",rownames(fitpval)[k])
                    newTable[z,2]<-fitpval[k,i]
                    z <- z +1
                }
            }
        } 
        
        
        
        
        return(newTable)
    })
    
    
    ###
    #creates summary table
    summaryInput <- reactive({
        req(input$file1)
        totalTableF <- create_table()
        fit <- eval(parse(text = paste("survfit(Surv(age, censored) ~ ", showcon(), ", data=totalTableF)")))
        #fit <- survfit(Surv(age, censored) ~ description, data=totalTableF)
        lok <- surv_median(fit)
        lok <- surv_median(fit)
        
        for( i in seq(1L,length(fit$n))){
            lok[i,5]<-fit$n[i]
        }
        lok[ , c(1,2,3,4,5)] <- lok[ , c(1,2,5,3,4)]
        
        
        namesCh <- gsub("description=","",names(fit$strata))
        
        cLabels <- character()
        varname <- character()
        
        for(i in 1:length(namesCh)){
            varname[i] <- paste0('labtxt_',gsub(" ", "", namesCh[i], fixed = TRUE))
            lok[i,1] <- input[[varname[i]]]
        } 
        
        colnames(lok)[3] <- "N"
        colnames(lok)[1] <- "description"
        colnames(lok)[4] <- "95% lower confidence limit"
        colnames(lok)[5] <- "95% upper confidence limit"
        
        return(lok)
    })
    
    ###
    #displays summary
    output$summary <- renderTable({
        summaryInput()
    },digits = 3)
    
    ###
    #displays  stats
    output$survfit <- renderTable({
        statsInput()
    },rownames = TRUE, digits = 3)
    
    
    ###
    #creates excel book object
    excelBook <- reactive({
        wrkBook <- createWorkbook()
        addWorksheet(wrkBook,"Summary")
        addWorksheet(wrkBook,"Statistics")
        addWorksheet(wrkBook,"Raw_Data")
        writeData(wrkBook,"Summary",summaryInput())
        writeData(wrkBook,"Statistics",statsInput())
        writeData(wrkBook,"Raw_Data",create_table())
        return(wrkBook)
    })
    
    ###
    #saves excel book object
    output$stats <- downloadHandler(
        filename = function() {
            paste(input$savename, "stats.xlsx", sep = "_")
        },
        content = function(file) {
            saveWorkbook(excelBook(),file,overwrite = TRUE)
        }
    )
    
    ###
    #saves graph handler
    output$graph <- downloadHandler(
        filename = function() {paste(input$savename, "_graph.",input$fileForm, sep = "")},
        #filename = "test.png",
        content = function(file) {
            ggsave(file, plot =  print(plotInput()), device = input$fileForm,width = input$graphX, height = input$graphY,dpi=1200)
        }
    )
    
    ###
    #outputs example table
    output$exampleTable <- renderTable({
        newTable <- tibble(column1 = c("%break","#title","#treatment",5,6,6,10,11,12,13,14,14,14,15,"%break","#title","#treatment",7,7,8,9,10,10,10,13), column2 = c("","CoolData","DMSO",1,1,0,1,1,1,1,1,1,0,1,"","CoolerData","EtOH",1,1,0,1,1,1,0,1))
        
    })
    
    exampleTableP <- reactive({
        newTable <- tibble(column1 = c("%break","#title","#treatment",5,6,6,10,11,12,13,14,14,14,15,"%break","#title","#treatment",7,7,8,9,10,10,10,13), column2 = c("","CoolData","DMSO",1,1,0,1,1,1,1,1,1,0,1,"","CoolerData","EtOH",1,1,0,1,1,1,0,1))
        colnames(newTable) <- c("%","%")
        return(newTable)
    })
    
    output$exampleFile <- downloadHandler(
        filename = "exampledataset.csv",
        content = function(file) {
            write.csv(exampleTableP(), file, row.names = FALSE, quote=FALSE)
        }
    )
    
    
}






# Create Shiny app ----
shinyApp(ui, server)