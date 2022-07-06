# =========== #
# Spatial App #
# =========== #
library("Seurat")
library("Matrix")
library("reticulate")
library("stringr")
library("ggplot2")
library("shinycssloaders")
library("shinyWidgets")
library("DT")
library("qs")
library("RColorBrewer")
library("patchwork")
library("shinyjs")
options(warn=-1)

# Load object
obj = qs::qread("data/all_merge.qs")
obj@active.assay = "Spatial"
# spo = qs::qread("data/all_obj_list.qs")

# Load Gene Information
gene_info = read.table("data/gene_info.txt", sep="\t", stringsAsFactors = F, header = T)
human_genes = sort(unique(gene_info$human))
human_gene_names = human_genes
human_logic = paste0("input.gene == '", human_genes,"' || ", collapse = '')
human_logic = substr(human_logic, 1, nchar(human_logic)-3)
human_logic2 = paste0("input.all_gene_human == '", human_genes,"' || ", collapse = '')
human_logic2 = substr(human_logic2, 1, nchar(human_logic2)-3)
all_genes = unique(c(gene_info$mzebra, gene_info$human))
gene_names <- rownames(obj@assays$Spatial)

sample_pt_size = c(4, 3, 3, 4, 3, 3, 3, 4, 3, 4, 4, 4, 3)
names(sample_pt_size) = levels(obj$sample)

# sample_checkbox_order = c("C2a", "B2a", "C2b", "B2b", "C2c", "B2c", "C2d", "B2d", "C1a", "C1b", "C1c", "C1d", "B1c")

num_clusters = max(as.vector(obj$all_umap_cluster2))
clusters = 1:num_clusters
cluster_choices = c(clusters)

bhve_samples = c("b1", "b2")
ctrl_samples = c("c1", "c2")

cur.num.ortho.cichlid = 0

# UI ===========================================================================
# Define UI for app
ui <- fluidPage(
  
  tags$head(
    tags$style(
      HTML(
        ".checkbox-inline { 
                    margin-left: 0px;
                    margin-right: 10px;
          }
         .checkbox-inline+.checkbox-inline {
                    margin-left: 0px;
                    margin-right: 10px;
          }
        "
      )
    )
  ),
  
  # App title
  # titlePanel("Paint Expression of Spots"),
  
  # Sidebar layout with input and output definitions
  # sidebarLayout(
  div(
    # Sidebar panel for inputs
    # sidebarPanel(
    wellPanel(
      tags$style(type="text/css", '#leftPanel { width:200px; float:left;}'),
      id = "leftPanel",
      
      # selectizeInput(inputId = "gene", label = "Gene Name", choices = NULL, selected = "fosb", multiple = TRUE, options = list(maxOptions = 1000)),
      
      selectizeInput(inputId = "all_gene_cichlid", label = "Cichlid Gene Name", choices = NULL, selected = "celsr1a", multiple = TRUE, options = list(maxOptions = 1000)),
      textOutput(outputId = "closest_human_gene"),
      selectizeInput(inputId = "all_gene_human", label = "Human Gene Name", choices = NULL, selected = "Celsr1", multiple = TRUE, options = list(maxOptions = 1000)),
      conditionalPanel(condition = "output.show_otho_cond_panel", id = "ortho.cond.panel",
                       radioButtons(inputId = "cichlid.ortho", label = "MZ Orthologs",
                       # checkboxGroupInput(inputId = "cichlid.ortho", label = "MZ Orthologs",
                                          choices = "junk",
                                          selected = "junk"),
                                          # inline = T),
                       tags$head(tags$style('#info{background-color: white;
                                    font-family: "Helvetica Neue",Helvetica,Arial,sans-serif;
                                    # border-style: solid;
                                    # border-wdith: thin;
                                    border-radius: 3px;
                                   }'))
      ),
      
      checkboxGroupInput(inputId = "samples", label = "Samples", 
                         choiceNames  = str_to_title(names(sample_pt_size)),
                         choiceValues = names(sample_pt_size),
                         selected = names(sample_pt_size), 
                         inline = T),
      checkboxGroupInput(inputId = "umap_add", label = "Disply w/ UMAP", 
                         choiceNames  = c("Comp. by Cluster", "Comp. by Sample"),
                         choiceValues = c("cluster_comp", "sample_comp"),
                         selected = c("cluster_comp", "sample_comp"), 
                         inline = F),
      sliderInput(inputId = "my.ncol", "Number of Columns", value = 4, min = 1, max = 13),
      sliderInput(inputId = "pt.size.multiplier", "Point Size", value = 0.5, min = 0, max = 1),
      checkboxInput(inputId = "toInfo", label = "Display Gene Info", value = F, width = NULL),
      
      # Gene Info Panel
      # conditionalPanel(condition = paste0("input.toInfo && input.gene != null && ", human_logic),
      #                  selectizeInput(inputId = "mz", label = "Orthologous MZ Genes", choices = NULL, multiple = TRUE, options = list(maxOptions = 1000))
      # ),
      conditionalPanel(condition = "input.toInfo && (input.all_gene_cichlid != null || input.all_gene_human != null)",
                       strong("Gene Info"),
                       textOutput(outputId = "info", container = pre),
                       tags$head(tags$style('#info{background-color: white;
                                    font-family: "Helvetica Neue",Helvetica,Arial,sans-serif;
                                    # border-style: solid;
                                    # border-wdith: thin;
                                    border-radius: 3px;
                                   }'))
      )
      # # conditionalPanel(condition = "input.all_gene_human != null && cur.num.ortho.cichlid > 1",
      
    ),
    
    # Main panel for displaying outputs
    # mainPanel(
    div(
      style = "flex-grow:1; resize:horizontal; overflow: hidden",
      tabsetPanel(id = "tabs", type = "tabs",
                  tabPanel("Spatial", value="sp_plot", plotOutput("sp_plot", width = "100%", height="100vh") %>% withSpinner(color="#FFA500")),
                  tabPanel("UMAP", value="umap_plot", plotOutput("umap_plot", width = "100%", height="100vh") %>% withSpinner(color="#FFA500"))
                  # tabPanel("UMAP w/ Info", value="umap_info_plot", plotOutput("umap_info_plot", width = "100%", height="100vh") %>% withSpinner(color="#FFA500"))
                  # tabPanel("By Fish", value="fishplot", plotOutput("fishplot", width = "100%", height="500px") %>% withSpinner(color="#0dc5c1")),
                  # tabPanel("BHVE vs CTRL", value="bcplot", plotOutput("bcplot", width = "100%", height="500px") %>% withSpinner(color="#0dc5c1")),
                  # tabPanel("BarPlot Summary", value="barplot", plotOutput("barplot", width = "100%", height="500px")  %>% withSpinner(color="#0dc5c1")),
                  # tabPanel("BarPlot %", value="barplotAvg", plotOutput("barplotAvg", width = "100%", height="500px")  %>% withSpinner(color="#0dc5c1")),
                  # tabPanel("Summary", value="summary", DT::dataTableOutput("summary") %>% withSpinner(color="#0dc5c1"))
      ),
      downloadButton(outputId = "down", label = "Download the plot")
      
      # conditionalPanel(condition = 'input.tabs == "summary || input.tabs == "barplot || input.tabs == "barplotAvg',
      #                  downloadButton(outputId = "sumDown", label = "Download Summary")
      # )
      
    )
  )
)

# Server =======================================================================
# Define server logic
server = function(input, output, session) {
  # updateSelectizeInput(session, "gene", choices = all_genes, server = TRUE)
  updateSelectizeInput(session, "all_gene_human",   choices = human_gene_names, server = TRUE)
  updateSelectizeInput(session, "all_gene_cichlid", choices = gene_names, server = TRUE)
  
  # Update MZ choices based on input Human gene - reactive
  # mz_choices = reactive({
  #   gene_info$mzebra[which(gene_info$human == input$gene)]
  # })
  # observe({
  #   updateSelectizeInput(session, "mz", choices = mz_choices(), server = TRUE, selected = mz_choices()[1])
  # })
  
  output$show_otho_cond_panel = reactive({ length(input$all_gene_human) > 0 })
  outputOptions(output, "show_otho_cond_panel", suspendWhenHidden = FALSE)
  
  # Update num MZ orthologs choices based on input Human gene - reactive
  observeEvent(input$all_gene_human, {
    if (length(input$all_gene_human) > 0) {
      myValues = gene_info$mzebra[which(gene_info$human == input$all_gene_human)]
      myENS = gene_info$ens[which(gene_info$human == input$all_gene_human)]
      myENS[which(!is.na(myENS))] = paste0(" (", myENS[which(!is.na(myENS))], ")")
      myENS[which(is.na(myENS))] = ""
      pat.mart.agree = gene_info$human_pat[which(gene_info$human == input$all_gene_human)] == gene_info$human_mart[which(gene_info$human == input$all_gene_human)]
      print(pat.mart.agree)
      myNames = paste0(myValues, myENS)
      myNames[which(pat.mart.agree)] = paste0(myNames, " ✓")
      updateRadioButtons(session = getDefaultReactiveDomain(), inputId = "cichlid.ortho", choiceNames = myNames, choiceValues = myValues, selected = findClosestCichlid(input$all_gene_human))
    }
  }, ignoreNULL = FALSE, ignoreInit = T)
  
  findClosestCichlid = function(human_gene) {
    this_cichlid_genes = gene_info[which(gene_info$human == human_gene),1]
    if (length(this_cichlid_genes) > 1) {
      upper_cichlid_gene <- this_cichlid_genes[which(startsWith(tolower(this_cichlid_genes), tolower(human_gene)))]
      if (length(upper_cichlid_gene) == 1) {
        cichlid_gene <- upper_cichlid_gene
      } else {
        cur.num.ortho.cichlid <<- length(upper_cichlid_gene)
        cichlid_gene <- this_cichlid_genes[1]
      } # end bad multiple
    } else {
      cichlid_gene <- this_cichlid_genes
    }
    
    return(cichlid_gene)
  }
  
  findClosestHuman = function(cichlid_gene) {
    this_human_genes = gene_info[which(gene_info$mzebra == cichlid_gene),2]
    if (length(this_human_genes) > 1) {
      upper_human_gene <- this_human_genes[which(startsWith(tolower(this_human_genes), tolower(cichlid_gene)))]
      if (length(upper_human_gene) == 1) {
        human_gene <- upper_mouse_gene
      } else {
        human_gene <- this_human_genes[1]
      } # end bad multiple
    } else {
      human_gene <- this_human_genes
    }
    
    return(human_gene)
  }
  
  
  findOvlPCells <- function(genes) {
    clean_genes <- c()
    pos_cells <- c()
    for (gene in genes) {
      gene_lower <- tolower(gene)
      gene_upper <- toupper(gene)
      gene_title <- str_to_title(gene)
      if (gene_lower %in% gene_names) {
        gene <- gene_lower
      } else if (gene_upper %in% gene_names) {
        gene <- gene_upper
      } else if (gene_title %in% gene_names) {
        gene <- gene_title
      }
      clean_genes <- c(gene, clean_genes)
      expr <- FetchData(object = obj, vars = gene, slot = "counts")
      pos_cells <- c(colnames(bb)[which(x = expr > 0)], pos_cells)
    }
    
    # Find the overlap of the positive cells
    num_genes <- length(genes)
    counts <- table(pos_cells)
    ovlp_cells <- rownames(counts[counts >= num_genes])
    obj <- SetIdent(obj, cells=ovlp_cells, value="overlapping_cells")
    obj <- SetIdent(obj, cells=setdiff(WhichCells(obj), ovlp_cells), value="non_overlapping_cells")
    obj$ovlp <- obj@active.ident
    return(obj$ovlp)
  }
  
  createBarPlot <- function(gene, average) {
    summary = createSummary(gene)
    
    df = data.frame()
    print(head(summary))
    for (cluster in clusters) {
      df = rbind(df, t(c("bhve", 
                         cluster, 
                         sum(summary[which(summary$Sample %in% bhve_samples & summary$Cluster == cluster),4]),
                         sum(summary[which(summary$Sample %in% bhve_samples & summary$Cluster == cluster),3]))))
      df = rbind(df, t(c("ctrl", 
                         cluster, 
                         sum(summary[which(summary$Sample %in% ctrl_samples & summary$Cluster == cluster),4]),
                         sum(summary[which(summary$Sample %in% ctrl_samples & summary$Cluster == cluster),3]))))
    }
    colnames(df) <- c("condition", "cluster_num", "cells_expr", "cells_in_category")
    df[,3] = as.numeric(as.vector(df[,3]))
    df[,4] = as.numeric(as.vector(df[,4]))
    df$pct = df$cells_expr / df$cells_in_category * 100
    df$value = df$cells_expr
    if(average)
      df$value = df$pct
    
    print(head(df))
    
    my_title <- paste("Number of Cells Expressing", paste(gene, collapse = ' and '), "per Cluster")
    my_ylab <- "Number of Cells"
    if (average == TRUE) {
      my_title <- paste("% Cells Expressing", paste(gene, collapse = ' and '), "per Cluster")
      my_ylab <- "% Cells"
    }
    p <- ggplot(df, aes(fill=condition, x=cluster_num, y=value)) +
      geom_bar(position="dodge", stat="identity") +
      theme_minimal() +
      ggtitle(my_title) +
      xlab("Cluster") +
      ylab(my_ylab) +
      # scale_x_continuous(breaks = clusters) +
      geom_text(aes(label=value), vjust=1.6, color="black", position = position_dodge(0.9), size=3.5)
    theme_minimal()
    p
  }
  
  createSplPlot <- function(cluster, resolution, dims, npcs, geneMode, origClust) {
    obj <- createSpltObj(cluster, resolution, dims, npcs, geneMode)
    if (origClust) {
      print("Displaying two plots")
      p1 <- DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 1.5) + ggtitle("New Clusters")
      Idents(obj) <- obj$orig.cluster 
      p2 <- DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 1.5) + ggtitle("Old Clusters")
      plot_grid(p1,p2)
    } else {
      DimPlot(obj, reduction = "umap", split.by = "cond", label = TRUE, pt.size = 1.5)
    }
  }
  
  createSplFeaturePlot <- function(gene, cluster, resolution, dims, npcs, geneMode) {
    obj <- createSpltObj(cluster, resolution, dims, npcs, geneMode)
    FeaturePlot(obj, features = c(gene), split.by = "cond", reduction = "umap", pt.size = 1.5, label=TRUE, order = TRUE)
  }
  
  createSpatialPlot = function(gene, split, samples) {
    # obj@active.assay <- "SCT"
    if (split == "cond") {
      FeaturePlot(obj, features = c(gene), split.by = "cond", reduction = "umap", pt.size = 1.5, label=TRUE, order = TRUE)
    } else {
      this.obj <- obj[,which(obj$sample %in% samples)]
      plist = SpatialFeaturePlot(this.obj, images = samples, features = gene, ncol = 4, combine = F)
      names(plist) = samples
      plist2 = list()
      gene.values = FetchData(object = this.obj, vars = gene)
      max.val = max(gene.values)
      min.val = min(gene.values)
      for (s in names(plist)) { plist[[s]]$layers[[1]]$aes_params$point.size.factor = sample_pt_size[s]*(input$pt.size.multiplier*2) }
      for (s in names(plist)) { plist2[[s]] = plist[[s]] + ggplot2::scale_fill_gradientn(colors = rev(brewer.pal(11, "Spectral")), limits = c(min.val, max.val)) }
      wrap_plots(plist2, ncol = input$my.ncol, guides = "collect") & theme(legend.position = "bottom")
      # & plot_annotation(theme = theme(plot.background = element_rect(fill ="#b5b5b5")))
    }
  }
  
  createUmapPlot = function(gene, split, samples, cluster_comp = F, sample_comp = F) {
    # obj@active.assay <- "SCT"
    this.obj <- obj[,which(obj$sample %in% samples)]
    plist = list()
    plist[[1]] = FeaturePlot(this.obj, features = gene, order = T, pt.size = 0.9*(input$pt.size.multiplier*2)) + theme_void() + coord_fixed() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    if (cluster_comp) {
      cluster.df = data.frame(table(this.obj$all_umap_cluster2))
      cluster.df$pos = data.frame(table(this.obj$all_umap_cluster2[which( this.obj@assays$Spatial@counts[gene,] > 0 )]))[,2]
      cluster.df$pct = (cluster.df$pos/cluster.df$Freq) * 100
      cluster.df$round_pct = format(round(cluster.df$pct, 1), nsmall = 1)
      cluster.df$round_pct[which(cluster.df$pct == 0)] = 0
      plist[[2]] = ggplot(cluster.df, aes(x = Var1, y = pct, fill = Var1)) + geom_bar(stat  = "identity") + xlab("Cluster") + ylab("% of Spots in Cluster") + scale_x_discrete(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + ggtitle("Comp. by Cluster") + theme_classic() + NoLegend() + geom_text(aes(label=round_pct), vjust=1.25, hjust = 0.5, color="black", position = position_dodge(0.9), size=3)
    }
    if (sample_comp) {
      sample.df = data.frame(table(this.obj$sample))
      sample.df$pos = data.frame(table(this.obj$sample[which( this.obj@assays$Spatial@counts[gene,] > 0 )]))[,2]
      sample.df$pct = (sample.df$pos/sample.df$Freq) * 100
      sample.df$round_pct = format(round(sample.df$pct, 1), nsmall = 1)
      sample.df$round_pct[which(sample.df$pct == 0)] = 0
      plist[[length(plist)+1]] = ggplot(sample.df, aes(x = Var1, y = pct, fill = Var1)) + geom_bar(stat  = "identity") + xlab("Sample") + ylab("% of Spots in Sample") + scale_x_discrete(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + ggtitle("Comp. by Sample") + theme_classic() + NoLegend() + geom_text(aes(label=round_pct), vjust=1.25, hjust = 0.5, color="black", position = position_dodge(0.9), size=3)
    }
    
    if (! (cluster_comp && sample_comp)) {
      wrap_plots(plist, ncol = 1)
    } else {
      # (plist[[1]] | (plist[[2]] / plist[[3]])) 
      wrap_plots(plist, ncol = 1)
    }
    
  }
  
  createOvlpPlot <- function(gene, split) {
    if (length(gene) > 1) {
      
      obj$ovlp <- findOvlPCells(gene)
      
      Idents(obj) <- obj$ovlp
      obj <- obj
      if (split == "sample") {
        Idents(obj) <- obj$sample
        only_b1_c1 <- obj[,WhichCells(obj,idents = c("b1", "c1"))]
        obj <- only_b1_c1
      }
      
      Idents(obj) <- obj$ovlp
      DimPlot(obj, reduction="umap", group.by = "ident", split.by=split, pt.size=2, order=TRUE)
      
    } else if (length(gene) == 1) {
      Idents(obj) <- obj$seurat_clusters
      obj <- obj
      if (split == "sample") {
        Idents(obj) <- obj$sample
        only_b1_c1 <- obj[,WhichCells(obj,idents = c("b1", "c1"))]
        obj <- only_b1_c1
        Idents(obj) <- obj$seurat_clusters
      }
      # obj@active.assay <- "RNA"
      FeaturePlot(obj, features = c(gene), split.by = split, reduction = "umap", pt.size = 1.5, label=TRUE, order = TRUE) 
    }
  }
  
  createSummary <- function(gene) {
    obj$ovlp <- findOvlPCells(gene)
    Idents(obj) <- obj$ovlp
    
    sample_clust = paste("overlapping_cells", obj$sample, obj$seurat_clusters)
    ovlp_sample_clust = paste(obj$ovlp, obj$sample, obj$seurat_clusters)
    pos_levels = expand.grid("overlapping_cells", unique(obj$sample), unique(obj$seurat_clusters))
    pos_levels = paste(pos_levels[,1], pos_levels[,2], pos_levels[,3])
    
    ovlp_sample_clust = factor(ovlp_sample_clust, levels = pos_levels)
    ovlp_sample_clust = ovlp_sample_clust[which(! is.na(ovlp_sample_clust) )]
    
    df = as.data.frame(table(ovlp_sample_clust))
    df_all = as.data.frame(table(sample_clust))
    df$all = df_all[match(df[,1], df_all[,1]), 2]
    df$pct = df[,2] / df$all * 100
    df2 = data.frame(do.call('rbind', strsplit(as.character(df[,1]), ' ', fixed=TRUE)))
    df = cbind(df2[,2:3], df[,c(3,2,4)])
    colnames(df) = c("Sample", "Cluster", "# Cells in Category", "# Cells Expressing", "% Cells Expressing")
    
    return(df)
  }
  
  geneInfo = function(gene) {
    str = ""
    if (length(gene) > 0) {
      if (gene %in% gene_names) {
        human = gene_info$human[which(gene_info$mzebra == gene)]
        human_description = gene_info$human_description[which(gene_info$mzebra == gene)]
        mzebra_description = gene_info$mzebra_description[which(gene_info$mzebra == gene)] 
        str = paste0(str, "\nClosest Human Gene:\n", human, "\n")  
        str = paste0(str, "\nHuman Description:\n\"", human_description, '"\n')
        str = paste0(str, "\nMZ Description:\n\"", mzebra_description, '"\n') 
        
      } else {
        human_description = gene_info$human_description[which(gene_info$human == gene)]
        if (length(human_description) > 1) {
          str = paste0(str, "\n# of Close MZ Genes:\n", length(human_description), "\n")
        } else {
          mzebra = gene_info$mzebra[which(gene_info$human == gene)]
          mzebra_description = gene_info$mzebra_description[which(gene_info$human == gene)]
          str = paste0(str, "Closest MZ Gene:\n", mzebra, "\n")
          str = paste0(str, "\nMZ Description:\n\"", mzebra_description, '"\n')    
        }
        str = paste0(str, "\nHuman Description:\n\"", human_description[1], '"\n')
      }
    }
    return(str)
  }
  
  # geneParser = function() {
  #   if (length(input$gene) < 1) { return(NULL) }
  #   mz_gene = input$gene
  #   if (! input$gene %in% gene_names)
  #     mz_gene = input$mz
  #   return(mz_gene)
  # }
  
  # numCichlidOrtho = function() {
  #   if (length(input$all_gene_human) > 0) {
  #     return(length(which(gene_info$human == gene)))
  #   }
  #   return(0)
  # }
  # output$ortho.cichlid = numCichlidOrtho()
  
  output$info = renderText({
    
    if (length(input$all_gene_cichlid) > 0) {
      geneInfo(input$all_gene_cichlid)
    } else if (length(input$all_gene_human) > 0) {
      gene <- input$all_gene_cichlid
      geneInfo(input$all_gene_human)
    }
    
  })
  
  output$closest_cichlid_gene <- renderText({
    if (length(input$all_gene_human) > 0) {
      str <- "Closest Cichlid Gene: "
      str <- paste(str, as.character(findClosestCichlid(input$all_gene_human)))
      return(str)
    }
  })
  
  output$closest_human_gene <- renderText({
    if (length(input$all_gene_cichlid) > 0) {
      str <- "Closest Human Gene: "
      str <- paste(str, as.character(findClosestHuman(input$all_gene_cichlid)))
      return(str)
    }
  })
  
  # Normal Plot (All Samples)
  output$sp_plot <- renderPlot({
    
    if (length(input$all_gene_cichlid) > 0) {
      gene <- input$all_gene_cichlid
      suppressMessages(createSpatialPlot(gene, "split", input$samples))
    } else if (length(input$all_gene_cichlid) == 0 && length(input$all_gene_human) > 0 && length(input$cichlid.ortho) > 0 && input$cichlid.ortho %in% gene_info$mzebra[which(gene_info$human == input$all_gene_human)]) {
      print(input$cichlid.ortho)
      gene = input$cichlid.ortho
      suppressMessages(createSpatialPlot(gene, 'split', input$samples))
    }

  })
  
  output$umap_plot <- renderPlot({
    if (length(input$all_gene_cichlid) > 0) {
      gene <- input$all_gene_cichlid
      suppressMessages(createUmapPlot(gene, "split", input$samples, cluster_comp = "cluster_comp" %in% input$umap_add, sample_comp = "sample_comp" %in% input$umap_add))
    } else if (length(input$all_gene_cichlid) == 0 && length(input$all_gene_human) > 0 && length(input$cichlid.ortho) > 0 && input$cichlid.ortho %in% gene_info$mzebra[which(gene_info$human == input$all_gene_human)]) {
      print(input$cichlid.ortho)
      gene = input$cichlid.ortho
      suppressMessages(createUmapPlot(gene, "split", input$samples, cluster_comp = "cluster_comp" %in% input$umap_add, sample_comp = "sample_comp" %in% input$umap_add))
    }
  })
  
  # Plots cells that are overlapping if >1 genes provided (BHVE vs CTRL)
  output$bc_ovlp_plot <- renderPlot({
    gene = geneParser()
    if (length(gene) > 0) {
      createOvlpPlot(gene, "cond")
    }
  })
  
  output$barplot <- renderPlot({
    gene = geneParser()
    if (length(gene) > 0) {
      createBarPlot(gene, FALSE)
    }
  })
  
  output$barplotAvg <- renderPlot({
    gene = geneParser()
    if (length(gene) > 0) {
      createBarPlot(gene, TRUE)
    }
  })
  
  output$splClust <- renderPlot({
    if (length(input$cluster) > 0 || length(input$gene) > 0) {
      gene = geneParser()
      createSplPlot(input$cluster, input$resolution, input$dims, input$npcs, input$geneMode, input$origClust)
    }
  })
  
  output$pntSplClust <- renderPlot({
    if (length(input$cluster) > 0 && length(input$gene) > 0) { 
      gene = geneParser()
      createSplFeaturePlot(gene, input$cluster, input$resolution, input$dims, input$npcs, input$geneMode)
    }
  })
  
  output$summary <- DT::renderDataTable({
    gene = geneParser()
    if (length(gene) > 0) {
      DT::datatable(createSummary(gene))
    }
  })
  
  # Download Functionality
  output$sumDown <- downloadHandler(
    filename = function() {
      gene = geneParser()
      paste0(paste0(c("summary", input$cluster, gene), collapse = "_"), ".tsv")
    }, 
    content = function(filename) {
      gene = geneParser()
      summary <- createSummary(gene)
      write.table(summary, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  )
  
  output$degClustCond <- downloadHandler(
    filename = function() {
      gene = geneParser()
      paste0(paste0(c("degByClustCond", input$cluster, gene), collapse = "_"), ".tsv")
    }, 
    content = function(filename) {
      if (input$tabs == "splClust" || input$tabs == "pntSplClust") {
        obj <- createSpltObj(gene, input$cluster, input$resolution, input$dims, input$npcs, input$geneMode)
        deg <- data.frame()
        obj$clust.cond <- paste0(obj$seurat_clusters, obj$cond)
        obj_num_clusters <- as.numeric(tail(levels(obj@meta.data$seurat_clusters), n=1))
        for (i in 0:obj_num_clusters) {
          newRow <- FindMarkers(obj, ident.1 = paste0("CTRL_", i), ident.2 = past0("BHVE_", i))
          newRow$gene <- rownames(newRow)
          newRow$cluster <- i
          deg <- rbind(deg, newRow)
        }
        write.table(deg, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
      }
      
    }
  )
  output$degClust <- downloadHandler(
    filename = function() {
      gene = geneParser()
      paste0(paste0(c("degByClust", input$cluster, gene), collapse = "_"), ".tsv")
    }, 
    content = function(filename) {
      if (input$tabs == "splClust" || input$tabs == "pntSplClust") {
        gene = geneParser()
        obj <- createSpltObj(gene, input$cluster, input$resolution, input$dims, input$npcs, input$geneMode)
        deg <- FindAllMarkers(obj)
        deg$gene <- rownames(deg)
        write.table(deg, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
      }
      
    }
  )
  
  output$down <- downloadHandler(
    filename = function() {
      gene = geneParser()
      if (length(gene) > 1) {
        paste(paste(gene, collapse = '_'), ".png", sep="")
      } else {
        paste(gene, ".png", sep="")
      }
    }, 
    content = function(filename) {
      print(filename)
      print(input$gene)
      if (input$tabs == "bc_ovlp_plot") {
        png(filename = filename, width = 900, height = 500, type="cairo")
        p <- createOvlpPlot(gene, "cond")
      } else if (input$tabs == "b1c1_ovlp_plot") {
        png(filename = filename, width = 900, height = 500, type="cairo")
        p <- createOvlpPlot(gene, "sample")
      } else if (input$tabs == "allplot") {
        height <- 500 * length(input$gene)
        png(filename = filename, width = 900, height = height, type="cairo")
        p <- createPlot(gene, "sample", c("b1", "b2", "b3", "b4", "b5", "c1", "c2", "c3", "c4", "c5"))
      } else if (input$tabs == "bcplot")  {
        height <- 500 * length(input$gene)
        png(filename = filename, width = 900, height = height, type="cairo")
        p <- createPlot(gene, "cond", c())
      } else if (input$tabs == "barplot") {
        png(filename = filename, width = 900, height = 500, type="cairo")
        p <- createBarPlot(gene, FALSE)
      } else if (input$tabs == "barplotAvg") {
        png(filename = filename, width = 900, height = 500, type="cairo")
        p <- createBarPlot(gene, TRUE)
      } else if (intput$tabs == "splClust") {
        png(filename = filename, width = 900, height = 500, type="cairo")
        p <- createSplPlot(input$cluster, input$resolution, input$dims, input$npcs, input$geneMode, input$origClust)
      } else if (intput$tabs == "pntSplClust") {
        png(filename = filename, width = 900, height = 500, type="cairo")
        p <- createSplPlot(input$cluster, input$resolution, input$dims, input$npcs, input$geneMode)
      }
      print(p)
      dev.off()
    }
  )
  
}

shinyApp(ui = ui, server = server)