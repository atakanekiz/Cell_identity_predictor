ui <- fluidPage(
  
  shinyjs::useShinyjs(), # needed for download button to work
  
  # tags$head(tags$style(".down_but{background-color:#39F407;} .down_but{color:#14150F;}")), THIS ALSO WORKS TO CONDITIONALLY SHOW DOWNLOAD BUTTON
  
  tags$head(includeHTML(("data/google_analytics.html"))),
  
  titlePanel("Cell identity predictor"),
  
  sidebarLayout(
    
    sidebarPanel(width = 2,
                 
                 # Cluster markers (Diff exp analysis, columns must have "gene", "avg_logFC", and "cluster" data)
                 fileInput("de_file", "Upload cluster diff. exp. data",
                           multiple = F,
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                 
                 # horizontal line
                 tags$hr(), 
                 
                 checkboxInput(inputId = "example_data", 
                               label = "Use example data (Ekiz HA and Huffaker TB, 2018)", 
                               value = F),
                 
                 # horizontal line
                 tags$hr(), 
                 
                 radioButtons("sel_reference", 
                              label = "Select reference data", 
                              choices = c("ImmGen", "Custom"), 
                              selected = "ImmGen"), 
                 
                 
                 conditionalPanel("input.sel_reference == 'Custom'",
                                  uiOutput("ui_sel_ref")),
                 conditionalPanel("input.sel_reference == 'Custom'",
                                  uiOutput("ui_sel_ref_annot")),
                 
                 actionButton("run", label="Analyze"),
                 
                 br(),br(), br(),
                 
                 downloadButton("download_res", "Download results"),
                 
                 br(),br(),
                 
                 downloadButton("download_top5", "Download top5 hits")
                 
                 # uiOutput("downloadData_ui") #THIS ALSO WORKS TO SHOW BUTTON AFTER ANALYSIS
                 
                 
                 
                 
                 
    ),
    
    
    
    mainPanel(
      
      
      
      tabsetPanel(
        
        tabPanel("Top-5 hits",
                 
                 fluidRow(
                   # h4("You can draw a rectangle around points for further information"),
                   plotOutput("top5", brush = "brushtop5"), height="600px"),
                 tableOutput("brushtop5")
        ),
        
        
        tabPanel("Individual Clusters",
                 
                 # This is the dynamic UI for the plots generated for each cluster
                 
                 fluidRow(
                   #   verbatimTextOutput("brush"),    )
                   h4("Please wait while the images are being generated. This can take some time depending on the cluster numbers.\n"),
                   p("In the plots below, horizontal lines indicate the mean identity score for each experimental cell cluster, and the shaded bands indicate +/- 1 and 2 standard deviations from the mean.\n\n"),
                   uiOutput("plots")

                 )       
        ),
        
        

        tabPanel("How to use this program",
                 
                 h3("Summary"),
                 p("Understanding the biological identity of cell clusters in single cell RNA sequencing (SCseq) experiments can be challenging due to overlapping gene expression profiles. An accurate assessment of the cluster identity requires analyzing multiple differentially expressed genes as opposed to a low throughput approach of examining a few selected lineage-specific markers. This program compares user-provided SCseq cluster gene signatures with", a('Immunological Genome Project (ImmGen)', href= 'https://www.immgen.org'), "immune cell datasets -- or with a user-provided reference dataset --, and calculates an", strong('identity score (IS)'), "for each SCseq cluster. To calculate the IS, our algorithm considers both upregulated and downregulated genes and the extent of differential expression, thus provide a detailed analysis of the cell phenotypes. For further information, please read the sections below."),
                 
                 br(),br(),
                 
                 p(strong(span(style="color:green", "To compare your SCseq clusters against ImmGen cell types, simply upload the genes differential expressed in clusters and click 'Analyze'"))),
                 br(),
                 p(strong(span(style="color:blue", "If you'd like to compare SCseq data with a custom reference, please select the appropriate radio button and upload reference gene expression data"))),
                 br(),
                 p(strong(span(style="color:gray", "You can also run the program by using example data that our lab has generated (Ekiz HA and Huffaker TB, in submission, 2018). These data were obtained by SCseq of murine tumor-infiltrating immune cells. Briefly, B16F10-OVA murine melanoma cells were subcutaneously injected into either wild-type mice or mice in which microRNA-155 is deleted only in the T cells via CD4-Cre. 9 or 12 days after injection, tumors were excised and disrupted mechanically and enzymatically with Accumax. CD45+ tumor infiltrating live immune cells were stained and sorted via FACS and were subjected to SCseq via 10X Genomics 3' platform. Data were analyzed by using Seurat R package, and differential gene expression in clusters were calculated after aggregating all samples in the experiment (d9/d12 WT and d9/d12 miR-155 T cell conditional knockout)"))),
                 
                 br(),br(),
                 
                 
                 h3("About your data"),
                 tags$ul(
                   tags$li(" Provide the expression data from SCseq cell clusters as a csv file."),
                   tags$li(" This file", strong("must have at least three columns named Gene, Cluster and logFC"), ". Capitalization does not matter and extra characters before or after these column names are also allowed."),
                   tags$li(" The data file can be as big as you want and can also contain other columns, although this will not be utilized by the algorithm."),
                   tags$li(" Each row of this data frame represents differentially expressed genes in cluster. You can generate this file manually, or by using", a("Seurat package's",href="https://satijalab.org/seurat/"), code("FindAllMakers()"), "function")
                   
                 ),
                 br(),br(),
                 p(strong("Sample SCseq differential expression data (generated via Seurat)")),
                 imageOutput("sample_data_file"),
                 
                 
                 
                 h3("About reference dataset"),
                 tags$ul(
                   tags$li("This program is pre-loaded with ImmGen microarray data for easily investigating immune cells in SCseq experiments"),
                   tags$li("If one would like to provide a custom reference dataset, any microarray or RNAseq data can be used, as long as data are normalized together."),
                   tags$li("Custom reference datasets can be uploaded as a csv file and should contain gene expression data from known cell types"),
                   tags$li("Reference dataset should have genes in rows and cell types in columns"),
                   tags$li(strong("The first column of the reference data frame must have gene names and must have string 'gene' in column name"),".Capitalization does not matter and characters before and after 'gene' is allowed"),
                   tags$li("This file can contain biological and technical replicates. In this case, IS is calculated and plotted for each replicate separately")
                 ),
                 br(),br(),
                 p(strong("Sample reference gene expression data")),
                 imageOutput("sample_reference_file"),
                 
                 
                 h3("About (optional) reference annotation data"),
                 tags$ul(
                   tags$li("If a custom reference dataset is used, although not necessary, an annotation file (in csv format) can be uploaded to obtain more informative plots"),
                   tags$li(strong("Annotation file must contain the columns named as 'short_name', 'long_name', 'description', and 'ref_cell_type'.")),
                   tags$li(strong("Data under 'short_name' column must EXACTLY match the column names of the reference gene expression data.")),
                   tags$li("Annotation file can have other columns as well, which will be shown as a tabular format along with other columns after calculations")
                 ),
                 br(),br(),
                 p(strong("Sample reference annotation data")),
                 imageOutput("sample_annotation_file"),
                 
                 
                 h3("How is cell identity score calculated?"),
                 p("The algorithm works first by calculating gene expression signatures in the reference file, and then comparing the experimental cluster signatures with these reference signatures. Specifically the following processes are performed:"),
                 h5(strong("Pre-processing of reference dataset")),
                 tags$ul(
                   tags$li("The algorithm uses pre-loaded ImmGen signatures or a user-uploaded reference expression file (see below for the details of reference file format)."),
                   tags$li("The reference expression file should contain normalized gene expression values (microarray or RNAseq) in rows and cell types in columns."),
                   tags$li("For each gene found in the reference file, algorithm calculates the ratio of gene expression in each cell type (i.e. individual columns) to the average gene expression value of the whole data set (i.e. all columns)."),
                   tags$li("Ratio is log transformed to obtain positive values indicating upregulation and negative values indicating downregulation . Therefore, pre-processing of the reference file results in a data frame that features logFC values of each gene for each cell type.")
                 ),
                 h5(strong("Analysis of experimental cluster signatures against reference file")),
                 tags$ul(
                   tags$li("Uploaded cluster signature file should contain information about genes, logFC in expression levels and the cluster information (see below for allowed file format)."),
                   tags$li("Algorithm finds the common genes between experimental data and the reference dataset."),
                   tags$li("For each shared gene, logFC values of differentially expressed in clusters and are multiplied with the logFC values in the reference dataset. This way if a gene is upregulated or downregulated in both a given cluster and a cell type in the reference dataset (i.e. positive correlation), the multiplication will result in a positive number (i.e. multiplication of two positive or two negative numbers resulting in a positive number). Alternatively, if a gene is differentially regulated in opposite directions in cluster and reference cell type (i.e. negative correlation), multiplication of logFC values will result in a negative number."),
                   tags$li("Multiplied logFC values per each gene are added up, resulting in an aggregate IS across the dataset for each cluster in the experimental data. This way, genes that have similar expression patterns in the experimental cell cluster and the reference cell type contribute to a higher IS, whereas genes with opposite expression patterns will result in a lower IS."),
                   tags$li(strong("This way, each cluster in the SCseq experiment is analyzed against each known cell type in the reference file and scored for its overall similarity. A higher identity score corresponds to a positive correlation in gene expression profiles and hence higher chances of a cluster belonging to the particular cell type found in the reference file")),
                   tags$li("For each cell cluster in the experiment, aggregate IS of reference cell types are plotted in dot plots which shows reference cell types in the x-axis, and aggregate IS in the -axis."),
                   tags$li("A summary plot showing the reference cell types correspoint to the highest 5 IS score is also shown.")
                 ),
                 
                 
                 h3("How confident is the prediction?"),
                 p("Since our algorithm compares unknown cluster gene signatures with signatures of reference samples, the confidence in the identity prediction can only be as high as the biological overlap between the experimental sample and the reference datasets. In an ideal scenario where gene expression data are available from all possible cell types under various experimental conditions, we envision that this algorithm can describe the identity of the unknown cell clusters in a highly accurate manner. However, in the absence of such data, our algorithm is still useful in characterizing unknown cluster identities by using a multiparametric approach. To help assess the confidence of the predictions, we utilize two different metrics:"),
                 tags$ul(
                   tags$li(strong("Deviation from mean:"), "For each unknown cluster, this metric is calculated by", strong(em("i) ")), "Averaging the identity scores from the whole reference dataset", strong(em("ii) ")), "Subtracting the identity scores of individual reference from this average", strong(em("iii) ")), "Dividing the difference by the standard deviation of the identity score across the whole reference dataset. This way, the deviation from mean is reported in the units of ", em("standard deviation."), "The higher deviation indicates a more pronounced distinction from the rest of the reference cell types, hence a higher confidence in prediction."),
                   tags$li(strong("Frequency of positively correlated genes:"), "This metric calculates the percentage of signature genes from unknown clusters that are regulated in a similar fashion in the reference cell types (i.e. upregulated or downregulated in both datasets). To exemplify this, please consider the following example: Imagine an unknown cluster, 'Cluster_1', which is defined by 100 differentially expressed genes (90 of which are upregulated, 10 are downregulated). If a reference cell type, 'Ref_cell_A', also has these 100 genes regulated in the same direction (90 up, 10 down), the frequency of positively correlated genes would be 100%. On the other hand, if 'Ref_cell_B' had all of these genes upregulated, then this frequency would be 90% (due to the 10 genes exhibiting divergent regulation). This way, the investigator can assess the overlap of signature genes that show similar patterns of expression.")
                 ), # close tags$ul
                 p("Although these metrics can be helpful in determining the confidence of prediction, they are not suitable to serve as a benchmark for the performance of our algorithm. A 'low-confidence' prediction can be explained by various factors including"),
                 tags$em(
                   tags$ul(
                     tags$li("A lack of comparable cell types in the reference dataset"),
                     tags$li("Using suboptimal number of metagene dimensions during clustering step of data analysis (which can result in  impure clusters composed of multiple cell types"),
                     tags$li("Presence of outliers in experimental data and clustering artifacts.")
                   )),
                 
                 p("Alternatively, a 'low-confidence' prediction can point to the discovery of new and interesting cell populations, and identify contaminating cell populations in the experiment. For instance, our experience with tumor-infiltrating lymphocyte single cell RNA-sequencing in murine B16F10 melanoma model suggests that, contaminating melanoma cells are loosely identified as ImmGen blood stem cells due to their proliferative capacity and overlapping gene expression profiles (Ekiz HA, manuscript in preparation)."),
                 
                 h3("Contact us"),
                 p("For questions and comments:"),
                 p("atakan.ekiz@hci.utah.edu"),
                 p("Ryan O'Connell Lab"),
                 p("Department of Pathology"),
                 p("University of Utah"),
                 p("Salt Lake City, Utah")
                 
                 
        )
        
        
        
      )
    ) 
    
    
    
  )
)

