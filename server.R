server <- function(input, output){
  
  
  options(shiny.maxRequestSize=50*1024^2) 
  
  
  output$brushtop5 <- renderTable({
    
    validate(
      need(input$brushtop5, "Draw a rectangle around data points for further information")
    )
    
    brushedPoints(top5_df_brush, input$brushtop5)
  }, striped = T)
  
  
  shinyjs::disable("download_res")
  shinyjs::disable("download_top5")
  
  observeEvent(analyzed_df(), {
    shinyjs::enable("download_res")
    shinyjs::enable("download_top5")
  })
  
  
  # output$downloadData_ui <- renderUI({
  #    req(analyzed_df())
  #    downloadButton("downloadData", "Download results", class="down_but")
  #  })  #THIS ALSO WORKS TO SHOW BUTTON AFTER ANALYSIS
  
  
  ################################################################################################################################
  # Define conditional dynamic file upload prompted when user selects "Custom" as reference data
  output$ui_sel_ref <- renderUI ({
    
    if (input$sel_reference == "Custom"){
      fileInput("ref_file", "Upload custom reference file",
                multiple = F, 
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv"))
    }
  })
  
  
  output$ui_sel_ref_annot <- renderUI ({
    
    if (input$sel_reference == "Custom"){
      fileInput("annot_file", "Upload custom annotation file",
                multiple = F, 
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv"))
    }
  })
  
  
  output$sample_data_file <- renderImage({
    
    list(src = "data/cluster_expr_IMG.png",
         alt = "Sample SCseq data",
         width=500)
    
  }, deleteFile = F)
  
  output$sample_reference_file <- renderImage({
    
    list(src = "data/ref_data_IMG.png",
         alt = "Sample reference data",
         width=500)
    
  }, deleteFile = F)
  
  output$sample_annotation_file <- renderImage({
    
    list(src = "data/ref_annot_IMG.png",
         alt = "Sample annotation data",
         width=500)
    
  }, deleteFile = F)
  
  ################################################################################################################################
  # Read uploaded differential expression file
  de_data <- reactive({
    
    inFile <- input$de_file
    
    
    
    if(is.null(inFile) & input$example_data == T){
      
      
      req(input$run)
      
      as_tibble(read.csv("data/Trimmed_cluster_signatures.csv",
                         check.names = F,
                         strip.white = T))
      
      
      
    } else {
      
      validate(
        need(input$de_file != "", "Please upload a data set or use example data")
      )
      
      # Make sure the file type is correct
      validate(
        need(tools::file_ext(inFile$name) %in% c(
          'text/csv',
          'text/comma-separated-values',
          'text/tab-separated-values',
          'text/plain',
          'csv',
          'tsv'
        ), "Wrong File Format try again!"))
      

      
      dat <- as_tibble(
        read.csv(inFile$datapath, check.names=TRUE, strip.white = TRUE)
      )
      
      # Make sure the column names are proper for correct subsetting
      validate(
        need(
          {if(sum(Reduce("|", lapply(c("logfc", "gene", "cluster"), grepl, colnames(dat), ignore.case=T))) == 3) TRUE else FALSE},
          "Formatting error: Make sure your dataset contains at least three columns named 'logfc', 'gene', and 'cluster' (capitalization is not important)"
      )
      )
      
      req(input$run)
      
      dat
      
    }
    
  }) # close de_data reactive object
  
  
  ################################################################################################################################
  # Read reference dataset
  
  ref_data <- reactive({
    
    if(input$sel_reference == "ImmGen"){
      
      reference <- as_tibble(read.csv("data/ImmGenData20180420_18_48_38.csv", check.names=FALSE, strip.white = TRUE))
      
      # Calculate row means for each gene (mean expression across 209 different immgen cell types)
      gene_avg <- rowMeans(reference[,2:dim(reference)[2]])
      
      
      # Calculate the ratio of gene expression in a given cell type compared 
      # to the average of the whole cohort. Calculate log (natural) fold change and store it in immgen_dat2
      reference_ratio <- log(sweep(reference[,2:dim(reference)[2]], 1, FUN="/", gene_avg))
      
      
      # Combine gene names and the log fold change in one data frame
      reference_log <- cbind(reference[,1], reference_ratio)
      
      reference_log
      
      
    } else if (input$sel_reference == "Custom"){
      
      in_refFile <- input$ref_file
      
      reference <- as_tibble(read.csv(in_refFile$datapath, check.names=FALSE, strip.white = TRUE))
      
      # # Can expand this code here further by allowing the upload of sample annotations to display details of reference samples
      # reference_annotation <- read.delim("popinfo.txt", header = F)
      # names(reference_annotation) <- c("Short.Name", "Long.Name", "Details", "Laboratory", "n?", "read?")
      # reference_annotation <- reference_annotation[,1:3]
      
      # Calculate row means for each gene (mean expression across 209 different immgen cell types)
      gene_avg <- rowMeans(reference[,2:dim(reference)[2]])
      
      
      # Calculate the ratio of gene expression in a given cell type compared 
      # to the average of the whole cohort. Calculate log (natural) fold change and store it in immgen_dat2
      reference_ratio <- log(sweep(reference[,2:dim(reference)[2]], 1, FUN="/", gene_avg))
      
      
      # Combine gene names and the log fold change in one data frame
      reference_log <- cbind(reference[,1], reference_ratio)
      
      reference_log
      
    } else {NULL}
    
    
  })
  
  ################################################################################################################################
  # Read immgen annotation file for explanations of cell types
  
  
  reference_annotation <- reactive({
    
    if(input$sel_reference == "ImmGen"){
      
      ref_annotation <- as_tibble(read.csv("data/ImmGen_annotations.csv", header = T))  # THE ORDER DOESN'T HAVE TO MATCH
      ref_annotation
      
    } else if(input$sel_reference == "Custom"){
      
      annotFile <- input$annot_file
      
      ref_annotation <- as_tibble(read.csv(annotFile$datapath, check.names=FALSE, strip.white = TRUE))
      ref_annotation
      
    }
  })
  
  ################################################################################################################################
  # Define a reactive cluster object that will store cluster information
  clusters <- reactive({
    
    
    gtools::mixedsort(
      levels(
        as.factor(
          pull(de_data(), grep("cluster", x = colnames(de_data()), ignore.case = T, value = T)
          )
        )
      )
    )
    
  }) # close clusters reactive object
  
  
  
  ################################################################################################################################
   # Compare de_data against reference file
  analyzed_df <- reactive({
    
    
    
    gene_column <<- grep("gene", colnames(de_data()), ignore.case = T, value = T)
    logFC_column <<- grep("logfc", colnames(de_data()), ignore.case = T, value = T)
    cluster_column <<- grep("cluster", colnames(de_data()), ignore.case = T, value = T)
    ref_gene_column <<- grep("gene", colnames(ref_data()), ignore.case = T, value = T)
    
    
    master_df <- data.frame()
    
    withProgress(message = 'Analysis in progress', value = 0, {
      
      for (i in clusters()) {
        
        # Increment the progress bar, and update the detail text.
        incProgress(1/length(clusters()), detail = paste("Analyzing cluster", i))
        
        sel_clst <- de_data() %>%
          filter(!!rlang::sym(cluster_column) == i) %>%
          select_(.dots = c(gene_column, logFC_column))
        
        
        # Merge SCseq cluster log FC value with immgen log FC for shared genes
        merged <- merge(sel_clst, ref_data(), by.x = gene_column, by.y = ref_gene_column)
        
        
        # Calculate a scoring matrix by multiplying log changes of clusters and immgen cells
        reference_scoring <- data.frame(apply(merged[,3:dim(merged)[2]],2,function(x){x*merged[,2]}), check.names = FALSE)
        
        # Calculate the aggregate score of each immgen cell type by adding
        score_sum <- colSums(reference_scoring)
        
        df <- data.frame(reference_score_sum = score_sum)
        
        df <- rownames_to_column(df, var="reference_id")
        
        # FIX THIS LATER AFTER UPDATING THE IMMGEN/CUSTOM REFERENCE UPLOADING CODE ABOVE SECTION
        
        if(input$sel_reference == "ImmGen"){
          
          df <- left_join(df, reference_annotation(), by=c("reference_id" = "short_name"))
          
          
          
          df$ref_cell_type <- c(rep("pro/pre-B", length(1:9)),
                                rep("B cell", length(10:24)),
                                rep("DC", length(25:47)),
                                rep("pDC", length(48:51)),
                                rep("Macrophage", length(52:70)),
                                rep("Monocyte", length(71:78)),
                                rep("Granulocyte", length(79:84)),
                                rep("pre-T cell", length(85:93)),
                                rep("T cell", length(94:126)),
                                rep("NKT", length(127:133)),
                                rep("T cell (act)", length(134:149)),
                                rep("gdT", length(150:171)),
                                rep("NK cell", length(172:183)),
                                rep("Epith/Endoth", length(184:196)),
                                rep("Stem cell", length(197:209)))
          
        } else if (input$sel_reference == "Custom" & !is.null(input$annot_file)){
          
          df <- left_join(df, reference_annotation(), by=c("reference_id" = "short_name"))
          
          
        } else if(input$sel_reference == "Custom" & is.null(input$annot_file)){
          
          df$ref_cell_type <- rep("NA_ref_cell_type", dim(ref_data())[2]-1)
          
        }
        
        
        

        df$cluster <- i
        
        # Add confidence-of-prediction calculations here and append to the df
        # Calculate the mean and standard deviation of the aggregate scores per reference cell type
        mean_score_sum <- mean(df$reference_score_sum)
        score_sum_sd <- sd(df$reference_score_sum)
        
        # Calculate the distance of the identity score from population mean (how many std devs apart?)
        df$z_score <- (df$reference_score_sum - mean_score_sum)/score_sum_sd
        
        # Calculate the proportion of the genes changing in the same direction between unknown cluster and reference cell type
        df$percent_pos_correlation <- {
          
          ngenes <- dim(reference_scoring)[1]
          
          pos_corr_vector <- numeric()
          
          for(i in 1:dim(reference_scoring)[2]){
            
            # Calculate number of genes positively correlated (upregulated or downregulated in both unk cluster and reference)
            pos_cor <- ( sum(reference_scoring[, i] > 0) / ngenes ) * 100
            
            pos_corr_vector <- c(pos_corr_vector, pos_cor)
            
          } #close for loop
          
          pos_corr_vector
          
        } # close expression 
        
        
        
        master_df <- rbind(master_df,df)
        
        
        
      } # close for loop that iterates over clusters
      
    })
    
    
    master_df
    
  }) # close analyzed_df reactive expression
  
  
  
  ################################################################################################################################
  # Generate plotting area dynamically for individual cluster plots
  
  
  # Insert the right number of plot output objects into the web page (https://gist.github.com/wch/5436415/)
  output$plots <- renderUI({
    plot_output_list <- lapply(clusters(), function(i) {
      plotname <- paste("plot", i, sep="")
      plotOutput(plotname, height = 500, width = 1800, brush = "brush") # optimize plotting area
    }
    ) # close lapply 
    
    # Convert the list to a tagList - this is necessary for the list of items to display properly.
    do.call(tagList, plot_output_list)
    
  }) # close output$plots renderUI
  
  
  
  ################################################################################################################################
  # Prepare individual plots
  
  
  
  # # Call renderPlot for each one. Plots are only actually generated when they
  # # are visible on the web page.
  
  observe({
    
    withProgress(message = 'Analyzing', value = 0, {
      
      for (i in clusters()) {
        
        # Increment the progress bar, and update the detail text.
        incProgress(1/length(clusters()), detail = paste("Cluster", i))
        
        # Need local so that each item gets its own number. Without it, the value
        # of i in the renderPlot() will be the same across all instances, because
        # of when the expression is evaluated.
        local({
          
          df_plot <- analyzed_df() %>%
            filter(cluster == i)
          
          score_mean <- mean(df_plot$reference_score_sum)
          score_sd <- sd(df_plot$reference_score_sum)
          
          
          
          my_i <- i
          plotname <- paste("plot", my_i, sep="")
          
          output[[plotname]] <- renderPlot({
            
            df_plot_brushed <<- df_plot
            
            p <- ggdotplot(df_plot, x = "reference_id", y="reference_score_sum", 
                           fill = "ref_cell_type", xlab=F, ylab="Reference identity score",
                           font.y = c(14, "bold", "black"), size=1, x.text.angle=90,
                           title = paste("Cluster:",my_i), font.title = c(15, "bold.italic"),
                           font.legend = c(15, "plain", "black"))+
              theme(axis.text.x = element_text(size=10, vjust=0.5, hjust=1))+
              geom_hline(yintercept=score_mean)+
              annotate("rect", xmin = 1, xmax = length(df_plot$reference_id),
                       ymin = score_mean-score_sd, ymax = score_mean+score_sd,
                       fill = "gray50", alpha = .1)+
              annotate("rect", xmin = 1, xmax = length(df_plot$reference_id),
                       ymin = score_mean-2*score_sd, ymax = score_mean+2*score_sd,
                       fill = "gray50", alpha = .1)
            
            
            
            # Old iteration using ggdotchart function. It reorders X axis.
            # p <- ggdotchart(df_plot, x = "reference_id", y="reference_score_sum", 
            #                 group = "ref_cell_type", color = "ref_cell_type", xlab=F, ylab="Reference identity score",
            #                 font.y = c(14, "bold", "black"),
            #                 dot.size = 3, title = paste("Cluster:",my_i), font.title = c(15, "bold.italic"),
            #                 font.legend = c(15, "plain", "black"))+
            #   theme(axis.text.x = element_text(size=10, vjust=0.5, hjust=1))
            
            
            print(p)
            
          }) # close renderPlot
          
        }) # close local
        
      } # close for loop 
      
    })  # close withProgress
    
    
  }) #close observe
  
  
  ################################################################################################################################
  # Prepare top5 summary plots
  
  top_df <- reactive({
    
    
    top5_df <- analyzed_df() %>%
    group_by(cluster) %>%    #cluster
    top_n(5, wt = reference_score_sum) %>%
    arrange(as.numeric(cluster), cluster, desc(reference_score_sum))
  
  
  ordered_cluster_levels <- gtools::mixedsort(levels(as.factor(top5_df$cluster)))
  
  
  top5_df$cluster <- factor(top5_df$cluster, levels = ordered_cluster_levels)
  
  
  top5_df$index <- 1:nrow(top5_df)
  
  top5_df <- select(top5_df, cluster,
                    ref_cell_type,
                    reference_id,
                    long_name,
                    description,
                    reference_score_sum,
                    index, everything())
  
  
  
  top5_df_brush <<- top5_df
  
  top5_df
  
  })
  
  
  output$top5 <- renderPlot({
    
   top_plot <- top_df()
    
    ggdotplot(top_plot, x="index", y="reference_score_sum", 
              fill = "cluster", size=2, x.text.angle=90, 
              font.legend = c(15, "plain", "black")) +
      scale_x_discrete(labels=top_plot$reference_id)+
      theme(axis.text.x = element_text(vjust=0.5, hjust=1))
    
  })
  
  
  ################################################################################################################################
  # Download results 
  
  
  output$download_res <- downloadHandler(
    filename = "Identity_scores_all.csv",
    content = function(file) {
      write.csv(analyzed_df(), file, row.names = FALSE)
    }
  )
  
  
  output$download_top5 <- downloadHandler(
    filename = "Identity_scores_top5.csv",
    content = function(file) {
      
    write.csv(top_df(), file, row.names = FALSE)
    }
  )
  
  
} # close server function


# shinyApp(ui=ui, server=server)