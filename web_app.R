```r
library(pheatmap)
library(shiny)
library(RPostgreSQL)

drv = dbDriver("PostgreSQL") 

# Establish database connection
con1 = dbConnect(
 ###xxx###
)

load(file = "/srv/shiny-server/pst_gt_heatmap/mat_type1.rdata")

# Plot genotype heatmap
gt_fig <- function(chr, start, end, continent, country, location, nuclear, bioproject) {
  search <- paste0(
    "SELECT * FROM (SELECT * FROM pst239_gt_24_9_27 WHERE POS>",
    start * 1000,
    " AND POS<",
    end * 1000,
    " AND CHROM=",
    chr,
    " ORDER BY RANDOM() LIMIT 1000 ) AS subquery ORDER BY ID"
  )
  
  gt_data <- dbGetQuery(con1, search)
  
  if (length(continent) == 0) {
    continent <- "none"
  }
  if (length(country) == 0) {
    country <- "none"
  }
  if (length(location) == 0) {
    location <- "none"
  }
  if (length(nuclear) == 0) {
    nuclear <- "none"
  }
  if (length(bioproject) == 0) {
    bioproject <- "none"
  }
  
  if ("All" %in% continent) {
    select1 <- 4:242
  } else {
    select1 <- which(mat_type[, 2] %in% continent) + 3
  }
  
  if ("All" %in% country) {
    select2 <- 4:242
  } else {
    select2 <- which(mat_type[, 4] %in% country) + 3
  }
  
  if ("All" %in% location) {
    select3 <- 4:242
  } else {
    select3 <- which(mat_type[, 6] %in% location) + 3
  }
  
  if ("All" %in% nuclear) {
    select4 <- 4:242
  } else {
    select4 <- which(mat_type[, 7] %in% nuclear) + 3
  }
  
  if ("All" %in% bioproject) {
    select5 <- 4:242
  } else {
    select5 <- which(mat_type[, 8] %in% bioproject) + 3
  }
  
  print(select1)
  print(select2)
  print(select3)
  print(select4)
  print(select5)
  
  if (
    "none" %in% continent &
      "none" %in% country &
      "none" %in% location &
      "none" %in% nuclear &
      "none" %in% bioproject
  ) {
    all_select <- 4:242
  } else {
    all_select <- Reduce(union, list(select1, select2, select3, select4, select5))
  }
  
  print(head(gt_data[, 1:3]))
  print(all_select)
  
  dt = as.matrix(gt_data[, all_select])
  print(head(dt))
  
  dt <- apply(dt, 2, as.numeric)
  rownames(dt) <- c(rep(" ", nrow(dt)))
  rownames(dt)[round(seq(1, nrow(dt), length.out = 4))] <- paste0(
    " --",
    gt_data[round(seq(1, nrow(dt), length.out = 4)), 3] / 1000,
    "kb"
  )
  
  return(dt)
  
  # pheatmap(
  #   dt,
  #   cluster_rows = FALSE,
  #   cluster_cols = TRUE,
  #   show_colnames = F,
  #   annotation_col = annotation,
  #   annotation_colors = annotation_color,
  #   color = c("2" = '#EDE38B', "1" = '#800080', "0" = '#FC8902')
  # )
}

anno_select <- function(mat_type, anno) {
  # print(anno)
  
  if (length(anno) == 0) {
    anno <- "none"
  }
  
  if ("All" %in% anno) {
    a <- mat_type[, 2:7]
  } else if (!"none" %in% anno) {
    a <- as.data.frame(mat_type[, which(colnames(mat_type) %in% anno)])
    rownames(a) <- rownames(mat_type)
    colnames(a) <- colnames(mat_type)[which(colnames(mat_type) %in% anno)]
    print(a)
  } else {
    a <- NA
  }
  
  return(a)
}

# Define UI for the application
ui <- fluidPage(
  tags$head(
    tags$title("Pst-genotype heatmap"),
    tags$style(HTML(
      "
      .custom-title {
        text-align: center;  /* Center text */
        background-color: #48A5DC;  /* Background color */
        font-weight: bolder;  /* Bold text */
        font-style: italic;  /* Italic text */
        color: white;  /* Text color */
        padding: 10px; /* Add padding */
        font-size: 36px;
      }
    "
    ))
  ),
  
  # Application title
  tags$div(class = "custom-title", div("Pst-genotype heatmap")),
  
  # Sidebar with plotting controls
  tabsetPanel(
    id = "tabset",
    
    tabPanel(
      "Plotting!",
      sidebarLayout(
        sidebarPanel(
          
          # selectInput("sample", "Select samples", unique(mat_type[,1]), multiple = TRUE),
          numericInput("chr", "Chr", value = 1),
          
          numericInput("start", "Start(kb)", value = 1),
          
          numericInput("end", "End(kb)", value = 100),
          
          selectInput(
            "continent",
            "continent",
            c("All", unique(mat_type[, 2])),
            multiple = TRUE
          ),
          
          selectInput(
            "country",
            "country",
            c("All", unique(mat_type[, 4])),
            multiple = TRUE
          ),
          
          selectInput(
            "location",
            "location",
            c("All", unique(mat_type[-which(is.na(mat_type[, 6])), 6])),
            multiple = TRUE
          ),
          
          selectInput(
            "nuclear",
            "haplotype",
            c("All", unique(mat_type[, 7])),
            multiple = TRUE
          ),
          
          selectInput(
            "bioproject",
            "bioproject",
            c("All", unique(mat_type[, 8])),
            multiple = TRUE
          ),
          
          selectInput(
            "anno",
            "annotation",
            c("All", colnames(mat_type)[2:length(colnames(mat_type))]),
            multiple = TRUE
          ),
          
          actionButton("simulate", "Go!"),
          downloadButton("download2", "Download.png")
        ),
        
        mainPanel(
```
