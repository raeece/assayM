library(DBI)
library(RSQLite)
library(shiny)
library(shinyFeedback)
library(DBI)
library(dplyr)
library(DT)
library(ggplot2)
library(Biostrings)
library(seqinr)
library(shinythemes)
library(data.table)


library(GenomicRanges)
library(Gviz)
library(rtracklayer)
library(png)
library(lubridate)
library(scales)


# dataset size
#select count(*) from (select distinct vmeta.strain from vmeta join variants on vmeta.strain=variants.strain);
#105807
#select count(*) from (select distinct strain from tmeta);
#205754
#select count(*) from (select distinct strain from variants);
#124564


# Create an ephemeral in-memory RSQLite database
conn <- dbConnect(RSQLite::SQLite(), "example.db")
fa <- readDNAStringSet("wuhan.fasta")
reference <- fa[[1]]


countriessql = "select distinct country from vmeta order by country;"
countriesdf <- dbGetQuery(conn, countriessql)
cselect <- c("Global", countriesdf$country)



ui <- fluidPage(
    titlePanel(fluidRow(column(
        3, img(
            height = 130,
            width = 300,
            src = "logo1.png"
        )
    ),
    column(
        9,
        h1("assayM - Monitor mutations in COVID-19 diagnostic assays"),
        h5("#last Updated 15-DEC-2020 based on 261,988 strains and 34 assays.")
    )))
    ,
    theme = shinytheme("sandstone"),
    tabsetPanel(
        id = "tabset",
        tabPanel("Common Assays for h-Cov", value = "tab1",
                 fluidRow(column(
                     12, textOutput("selected_var"), div(dataTableOutput("assaytbl"))
                 )),),
        tabPanel(
            "Check Your Assay",
            value = "tab2",
            fluidRow(column(8, shinyFeedback::useShinyFeedback())),
            fluidRow(column(8, htmlOutput("assayName"))),
            fluidRow(
                column(4, textInput("Forward", h3("Forward *"),
                                    value = "GACCCCAAAATCAGCGAAAT")),
                column(4, textInput("Probe", h3("Probe"),
                                    value = "ACCCCGCATTACGTTTGGTGGACC")),
                column(
                    4,
                    textInput("Reverse", h3("Reverse"),
                              value = "TCTGGTTACTGCCAGTTGAATCTG")
                ),
                actionButton("check", "Check")
            ),
            fluidRow(
                column(4, textInput("Forwardloc", h3("Start-End *"),
                                    value = "28287-28306")),
                column(4, textInput("Probeloc", h3("Start-End"),
                                    value = "28309-28332")),
                column(4, textInput("Reverseloc", h3("Start-End"),
                                    value = "28335-28358")),
                actionButton("checkloc", "Check")
            ),
            fluidRow(column(
                8, selectInput("dataset", "Dataset:", cselect)
            )),
            fluidRow(column(12, plotOutput('track', height = 300))),
            fluidRow(column(12, plotOutput('bars', height = 400)))
            
            
        ),
        tabPanel("Mutations catalog", value = "tab3",
                 fluidRow(column(
                     12, div(dataTableOutput("mtbl"))
                 )),),
        tabPanel("About", value = "tab4",
                 fluidRow(column(
                     12,
                     HTML(
                         "<h4>Developed by Raeece Naeem at <a href='https://pgl.kaust.edu.sa/'>Prof. Arnab Pain Group</a> ,<br/> King Abdullah University of Science and Technology.<br/>
                                 Funded by <a href='https://kaust.edu.sa/en/r3t-covid-19'>KAUST Rapid Research Response Taskforce (R3T).</a></h4>"
                     )
                 )),)
    )
    
    
    
    
    
)





# Define server logic required to draw a histogram
server <- function(input, output, session) {
    print("loading1")
    
    assay_name <- reactiveValues()
    
    output$assayName <- renderText({
        #query <- parseQueryString(session$clientData$url_search)
        paste("<h2>",
              assay_name$name,
              "</h2>",
              sep = "",
              collapse = ", ")
    })
    
    viewassay <- function(assayName) {
        assay_name$name <- assayName
        assaysql <-
            paste0("select * from assays where assay_name='", assayName, "'")
        assaylist <- dbGetQuery(conn, assaysql)
        fwd <- assaylist %>% filter(assay_direction == 'F')
        prob <- assaylist %>% filter(assay_direction == 'P')
        rev <- assaylist %>% filter(assay_direction == 'R')
        
        updateTextInput(session, "Forward", value = fwd[1, "sequence"])
        updateTextInput(session, "Probe", value = prob[1, "sequence"])
        updateTextInput(session, "Reverse", value = rev[1, "sequence"])
        
        updateTextInput(session, "Forwardloc", value = paste0(fwd[1, "start"], "-", fwd[1, "end"]))
        updateTextInput(session, "Probeloc", value = paste0(prob[1, "start"], "-", prob[1, "end"]))
        updateTextInput(session, "Reverseloc", value = paste0(rev[1, "start"], "-", rev[1, "end"]))
        
        updateTabsetPanel(session, "tabset", selected = "tab2")
        adf <-
            data.frame(start = c(fwd[1, "start"], prob[1, "start"], rev[1, "start"]),
                       end = c(fwd[1, "end"], prob[1, "end"], rev[1, "end"]))
        v[['assaydf']] <- adf
    }
    
    observe({
        queryparam <- parseQueryString(session$clientData$url_search)
        if (!is.null(queryparam[['assay']])) {
            print("query param check")
            viewassay(queryparam[['assay']])
            
        }
    })
    
    
    v <- reactiveValues()
    
    observeEvent(input$check, {
        fstart <- NA
        fend <- NA
        pstart <- NA
        pend <- NA
        rstart <- NA
        rend <- NA
        
        shinyFeedback::feedbackWarning("Forward",
                                       input$Forward == "" ,
                                       "Please enter forward")
        req(input$Forward)
        
        m1 <- matchPattern(input$Forward, reference)
        fstart = start(m1)
        fend = end(m1)
        print(input$Probe)
        if (!(input$Probe == "")) {
            m2 <- matchPattern(input$Probe, reference)
            pstart = start(m2)
            pend = end(m2)
            if (length(pstart) == 0L) {
                pstart <- NA
                pend <- NA
            }
        }
        if (!(input$Reverse == "")) {
            m3 <-
                matchPattern(reverseComplement(DNAString(input$Reverse)), reference)
            rstart = start(m3)
            rend = end(m3)
        }
        
        adf <-
            data.frame(start = c(fstart, pstart, rstart),
                       end = c(fend, pend, rend))
        
        v[['assaydf']] <- adf
        print(v[['assaydf']]$start)
        print(v[['assaydf']]$end)
        if (!is.na(fstart)) {
            updateTextInput(session, "Forwardloc", value = paste0(fstart, "-", fend))
        }
        else{
            updateTextInput(session, "Forwardloc", value = "")
        }
        if (!is.na(pstart)) {
            updateTextInput(session, "Probeloc", value = paste0(pstart, "-", pend))
        }
        else{
            updateTextInput(session, "Probeloc", value = "")
        }
        if (!is.na(rstart)) {
            updateTextInput(session, "Reverseloc", value = paste0(rstart, "-", rend))
        }
        else{
            updateTextInput(session, "Reverseloc", value = "")
        }
        #v$data <- runif(100)
    })
    
    observeEvent(input$checkloc, {
        shinyFeedback::feedbackWarning("Forwardloc",
                                       input$Forwardloc == "" ,
                                       "Please enter forward location")
        req(input$Forwardloc)
        f <- unlist(strsplit(input$Forwardloc, "-"))
        fstart <- as.integer(f[1])
        fend <- as.integer(f[2])
        p <- unlist(strsplit(input$Probeloc, "-"))
        pstart <- as.integer(p[1])
        pend <- as.integer(p[2])
        r <- unlist(strsplit(input$Reverseloc, "-"))
        rstart <- as.integer(r[1])
        rend <- as.integer(r[2])
        adf <-
            data.frame(start = c(fstart, pstart, rstart),
                       end = c(fend, pend, rend))
        v[['assaydf']] <- adf
        #v$data <- runif(100)
    })
    
    
    print("loading2")
    
    q <- NULL
    output$track <- renderPlot({
        print("rendering genome view")
        if (!is.null(v[['assaydf']])) {
            withProgress(message = '...', value = 0, {
                print(v[['assaydf']]$start)
                print(v[['assaydf']]$end)
                
                
                
                axis_track <- GenomeAxisTrack()
                seq_track <-
                    SequenceTrack(fa, options(ucscChromosomeNames = FALSE))
                data <-
                    data.frame(
                        chr = c(
                            'NC_045512',
                            'NC_045512',
                            'NC_045512',
                            'NC_045512',
                            'NC_045512',
                            'NC_045512',
                            'NC_045512',
                            'NC_045512',
                            'NC_045512',
                            'NC_045512',
                            'NC_045512'
                        ),
                        start = c(
                            266,
                            21563,
                            25393,
                            26245,
                            26523,
                            27202,
                            27394,
                            27756,
                            27894,
                            28274,
                            29558
                        ),
                        end = c(
                            21555,
                            25384,
                            26220,
                            26472,
                            27191,
                            27387,
                            27759,
                            27887,
                            28259,
                            29533,
                            29674
                        ),
                        id = c(
                            'orf1ab',
                            'S',
                            'ORF3a',
                            'E',
                            'M',
                            'ORF6',
                            'ORF7a',
                            'ORF7b',
                            'ORF8',
                            'N',
                            'ORF10'
                        ),
                        strand = c(
                            '+',
                            '+',
                            '+',
                            '+',
                            '+',
                            '+',
                            '+',
                            '+',
                            '+',
                            '+',
                            '+'
                        )
                    )
                data_g <-
                    with(data, GRanges(chr, IRanges(start, end), strand, id = id))
                data_track <-
                    AnnotationTrack(
                        data_g,
                        name = "Genes",
                        width = 1,
                        showFeatureId = T,
                        max.height = 1
                    )
                
                res <- NULL
                sql3clause = ""
                if (input$dataset != "Global") {
                    sql3clause = paste0("where vmeta.country='",
                                        input$dataset,
                                        "'")
                }
                sql3 <-
                    paste0(
                        "with recursive cnt(x) as (values(1) union all select x+1 from cnt where x<29903) select x,CASE WHEN mutations IS NULL THEN 0 ELSE mutations END mutations from cnt left outer join (select start,count(start) mutations from variants join vmeta on variants.strain=vmeta.strain ",
                        sql3clause,
                        " group by start) y on x=y.start"
                    )
                rs3 <- dbGetQuery(conn, sql3)
                incProgress(0.8, detail = "checking assay")
                
                assay <-
                    data.frame(
                        chr = c('NC_045512', 'NC_045512', 'NC_045512'),
                        start = v[['assaydf']]$start,
                        end = v[['assaydf']]$end,
                        id = c('F', 'P', 'R'),
                        strand = c('+', '+', '-')
                    )
                
                # remove rows with NAs
                assay <- assay[complete.cases(assay),]
                print(assay)
                assay_g <-
                    with(assay, GRanges(chr, IRanges(start, end), strand, id = id))
                assay_track <-
                    AnnotationTrack(
                        assay_g,
                        name = "Assay",
                        width = 1,
                        showFeatureId = T,
                        max.height = 1
                    )
                
                
                gr <-
                    GRanges(seqnames = "NC_045512",
                            ranges = IRanges(seq(1, 29903, len = 29903), width = 1))
                values(gr) <- rs3$mutations
                dtrack <-
                    DataTrack(
                        range = gr,
                        genome = "SarsCov2",
                        name = "Number of Strains",
                        type = "b"
                    )
                plotTracks(
                    list(
                        axis_track,
                        seq_track,
                        data_track,
                        assay_track,
                        dtrack
                    ),
                    chromosome = 'NC_045512',
                    from = min(assay$start - 10),
                    to = max(assay$end + 10)
                )
            })
        }
        
    })
    
    
    output$bars <- renderPlot({
        if (!is.null(v[['assaydf']])) {
            for (i in 1:nrow(v[['assaydf']])) {
                if (!is.na(v[['assaydf']][i, "start"])) {
                    st <- v[['assaydf']][i, "start"]
                    end <- v[['assaydf']][i, "end"]
                    q <-
                        c(q,
                          paste(
                              '(variants.start between  ',
                              st,
                              ' and ',
                              end,
                              ')'
                          ))
                }
            }
            qu <- paste(q, collapse = " or ")
            sql4 <-
                paste0(
                    "select count(strain),country,collection_date from (select distinct(vmeta.strain),vmeta.country,vmeta.collection_date from variants join vmeta on ( ",
                    qu,
                    " ) and vmeta.strain=variants.strain  group by vmeta.strain,vmeta.country,vmeta.collection_date) t group by country,collection_date"
                )
            rs4 <- dbGetQuery(conn, sql4)
            colnames(rs4) <- c("strains", "country", "collection_date")
            if (input$dataset != "Global") {
                rs4 <- rs4 %>% filter(country == input$dataset)
            }
            rs4$collection_date <-
                parse_date_time(rs4$collection_date, c("%Y-%m-%d", "%Y-%m"))
            rs41 <-
                rs4 %>% mutate(month_year = round_date(collection_date, "month")) %>% group_by(month_year) %>% summarise(n =
                                                                                                                             sum(strains)) %>% filter(month_year >= '2020-02-01')
            
            sql5 <-
                paste0(
                    "select count(tmeta.strain),tmeta.country,tmeta.collection_date from tmeta group by tmeta.country,tmeta.collection_date"
                )
            rs5 <- dbGetQuery(conn, sql5)
            colnames(rs5) <- c("strains", "country", "collection_date")
            if (input$dataset != "Global") {
                rs5 <- rs5 %>% filter(country == input$dataset)
            }
            rs5$collection_date <-
                parse_date_time(rs5$collection_date, c("%Y-%m-%d", "%Y-%m"))
            rs51 <-
                rs5 %>% mutate(month_year = round_date(collection_date, "month")) %>% group_by(month_year) %>% summarise(n =
                                                                                                                             sum(strains)) %>% filter(month_year >= '2020-02-01')
            rs45 <-
                rs41 %>% right_join(rs51, by = "month_year") %>% mutate(perc = n.x / n.y)
            rs45$perc[is.na(rs45$perc)] <- 0
            rs45$month_year <- as.Date(rs45$month_year)
            p <- ggplot(data = rs45, aes(x = month_year)) +
                geom_bar(aes(y = n.y), stat = "identity") +
                geom_bar(aes(y = n.x, fill = "#6495ED"), stat = "identity") +
                geom_text(aes(y = 0, label = paste0(round(
                    perc * 100, 1
                ), "%")), color = "black") +
                scale_fill_identity(name = 'Percentage of strains with mutations',
                                    guide = 'legend',
                                    labels = c('')) +
                scale_x_date(
                    date_breaks = "1 month",
                    labels = date_format("%b-%Y"),
                    limits = as.Date(c('2020-01-01', '2021-01-01'))
                ) +
                theme(legend.position = "top") + xlab('Month') + ylab('Number of strains')
            print(p)
            
            
        }
    })
    
    res <- NULL
    rs1 <- NULL
    rs3 <- NULL
    withProgress(message = '...', value = 0, {
        incProgress(0.2, detail = "fetching assays")
        print("loading3")
        sql1 <-
            "select assay_name,mismatch_strains,country,gene_target,assay_type from assay_summary;"
        rs1 <- dbGetQuery(conn, sql1)
        rs2 <-
            rs1 %>% mutate(link = paste0(
                '<a  target=_blank href=?assay=',
                assay_name,
                '>',
                assay_name,
                '</a>'
            ))
        incProgress(0.8, detail = "fetching variants")
        print("loaded3")
    })
    
    
    ##########
    buttonInput <- function(FUN, len, id, ...) {
        inputs <- character(len)
        for (i in seq_len(len)) {
            inputs[i] <- as.character(FUN(paste0(id, i), ...))
        }
        inputs
    }
    
    vals <- reactiveValues()
    
    vals$Data <- data.table(
        Assay = rs1$assay_name,
        "Mutant strains" = rs1$mismatch_strains,
        Country = rs1$country,
        "Gene target" = rs1$gene_target,
        "Assay Type" = rs1$assay_type,
        Action = buttonInput(
            FUN = actionButton,
            len = 27,
            id = 'button_',
            label = "View",
            onclick = 'Shiny.setInputValue(\"lastClick\",  this.id,{priority:\"event\"})'
        )
    )
    
    
    output$assaytbl <- DT::renderDataTable({
        DT = vals$Data
        datatable(DT, escape = F, selection = 'single')
    })
    observeEvent(input$lastClick, {
        selectedRow <- as.numeric(strsplit(input$lastClick, "_")[[1]][2])
        print(vals$Data[selectedRow, 1])
        viewassay(vals$Data[selectedRow, 1])
    })
    
    
    
    output$mtbl <- renderDT({
        withProgress(message = '...fetching mutations', value = 50, {
            msql <-
                "select vmeta.strain,country,collection_date,start location,case mutation when 'S' then 'SNP' when 'I' then 'Insertion' when 'D' then 'Deletion' END  mutation_type,ref,allele from variants join vmeta on vmeta.strain=variants.strain;"
            mrs <- dbGetQuery(conn, msql)
            DT::datatable(
                mrs,
                filter = 'top',
                options = list(searching = TRUE, paging = TRUE)
            )
        })
    })
    
    
    
    
    
    
}

# Run the application
shinyApp(ui = ui, server = server)
