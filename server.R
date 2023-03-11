# Copyright (C) 2023 Wenjie SUN <sunwjie@gmail.com>
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without
# modifications, as long as this notice is preserved.

library(shiny)
library(magrittr)
library(data.table)
library(ggplot2)
library(plotly)

## load data
## all data in in the ./data/ folder
meta_d = fread("./data/meta_data.tsv")
exp_d = fread("./data/expression.tsv")
pseudotime_lum_d = fread("./data/luminal_pseudotime.tsv")
pseudotime_bas_d = fread("./data/basal_pseudotime.tsv")
edges_d = fread("./data/edges_d.tsv")
color_setting_d = fread("./data/color_scheme.csv")

## color configuration
color_setting = color_setting_d$color
names(color_setting) = color_setting_d$name

## define ggplot theme
theme0 <- theme_classic() + theme(
    text = element_text(size = 15),
    line = element_line(size = 1),
    axis.line = element_line(size = 1),
    axis.ticks.length = unit(3, units = "mm"),
    axis.text.x = element_text(
        margin = margin(t = 2, unit = "mm")
        , angle = 60, vjust = 1, size = 15, hjust = 1),
    axis.text.y = element_text(margin = margin(r = 3, l = 5, unit = "mm")),
    legend.position = "right",
) 
theme1 = theme0 + theme(
    axis.text.x = element_text( 
        margin = margin(t = 2, unit = "mm")
        , angle = 0, vjust = 1, size = 12, hjust = 0.5)
)

## function to preprocess gene expression for plotting
get_gene_expression = function(exp_d, gene, meta_d) {
    exp_sub = exp_d[, c("rn", gene), with=F]
    merge(meta_d, exp_sub, by.x = "Row.names", by.y = "rn")
}

## cell cluster ordering
level_order = c(
    "E13_Mammary Epithelial cells",
    "E14_Mammary Epithelial cells",
    "E15_Mammary Epithelial Basal-like cells",
    "E15_Mammary Epithelial Hybrid cells",
    "E15_Mammary Epithelial Luminal-like cells",
    "P0_Mammary Epithelial Basal cells",
    "P0_Mammary Epithelial Luminal Progenitors cells",
    "P0_Mammary Epithelial Luminal HRhigh cells"
    )

## pesudotime plot for basal and luminal trajectories
check_gene_both = function(gene, exp_d, meta_d, cell_ord1, cell_ord2, level_order, color_setting) {

    ## example
    # gene = "Anxa1"
    # cell_ord = pesudotime_lum
    # meta_d = d
    # cell_ord1 = pesudotime_lum
    # cell_ord2 = pesudotime_bas
    

    setkey(exp_d, "rn")
    setkey(meta_d, "Row.names")
    library(stringr)
    library(ggsci)
    x = exp_d[cell_ord1][, gene, with=F][[1]]
    d1 = data.table(Pseudotime = seq_along(x), expression = x, cell_type = meta_d[cell_ord1, CellType])
    x = exp_d[cell_ord2][, gene, with=F][[1]]
    d2 = data.table(Pseudotime = seq_along(x), expression = x, cell_type = meta_d[cell_ord2, CellType])

    d1[, trajectory := "Luminal"]
    d2[, trajectory := "Basal"]
    d3 = rbind(d1, d2)
    d3[trajectory == "Basal", Pseudotime := - Pseudotime - 20]
    d3[trajectory == "Luminal", Pseudotime := Pseudotime + 20]

    d3$cell_type = factor(d3$cell_type, levels = level_order)

    # constants
    y_axis_begin  <- min(d3$expression) %>% floor
    y_axis_end    <- max(d3$expression) %>% ceiling
    total_ticks   <- 4

    tick_frame_y <- data.frame(ticks = seq(y_axis_begin, y_axis_end, by = 1), zero=0) 


    lab_frame_y <- data.frame(
        lab = c(min(tick_frame_y$ticks), max(tick_frame_y$ticks)), zero = 0
        )

    tick_sz_y <- 10

    g = ggplot(d3, aes(x = Pseudotime, y = expression)) + theme_classic() +
        # y axis
        geom_segment(x = 0, xend = 0, y = -10, yend = tail(tick_frame_y$ticks, 1), size = 0.5) +
        # y ticks
        geom_segment(data = tick_frame_y, aes(x = zero, xend = zero + tick_sz_y, y = ticks, yend = ticks)) +
        geom_point(aes(color=cell_type)) + labs(title = gene) + scale_color_npg() +
        scale_color_manual(values = color_setting) +
        geom_smooth(aes(linetype=trajectory), color = 'black', method = "gam", formula = y ~ s(x, bs = "cs")) +
        geom_text(data=tick_frame_y, aes(y=ticks, x=zero, label=ticks), family = 'Times', vjust=0.5, hjust = 1.5) +
        guides(colour = guide_legend(override.aes = list(size=5)))

    g = g + theme1 + theme(
        axis.text.y = element_blank(), 
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.line.x = element_line(arrow = arrow(ends = "both", length = unit(0.1, "inches")))) +
    labs(x = "Pseudotime", y = "Expression", linetype = "Trajectory", color = "Cell Type") 

    g
}

shinyServer(function(input, output) {

    ## choose a gene name; if empty, default is the cluster
    select_d = reactive({

        col_var = input$col
        
        if (col_var == "" | !col_var %in% names(exp_d)) {
            col_var = "CellType"
            d = meta_d
        } else {
            validate(
                need(col_var %in% colnames(exp_d), "The gene you searched is not in the dataset.")
            )

            d = get_gene_expression(exp_d, col_var, meta_d)
        }
        list(
            col_var = col_var,
            d = d
            )
    })

    ## 3d trajectory plot
    output$trajectory_3d = renderPlotly({

        ## prepare data
        d_l = select_d()
        d = d_l$d
        d$lum_score %<>% scale
        d$bas_score %<>% scale
        d$PC_2%<>% scale
        value = d[[d_l$col_var]]
        d = d[!grepl("Skin", CellType)]
        super_cell_d = d[, .(lum_score = median(lum_score), bas_score = median(bas_score), PC_2 = median(PC_2)), by = CellType]

        ## initialize figure
        if (is.character(value)) {
            fig <- plot_ly(colors=color_setting) 
        } else {
            fig <- plot_ly()
        }

        ## draw cell
        fig = add_trace(fig,
            x = d$lum_score,
            y = d$bas_score, 
            z = d$PC_2, 
            data = d,
            type = 'scatter3d', mode = 'markers', 
            hoverinfo = "skip",
            color = value,
            marker = list(
                size = 5,
                opacity = 0.5
                ),
            showlegend = T
        )

        ## add edges
        for (i in 1:nrow(edges_d)) {
            node = edges_d[i, ] %>% unlist
            x = super_cell_d[CellType %in% node, lum_score]
            y = super_cell_d[CellType %in% node, bas_score]
            z = super_cell_d[CellType %in% node, PC_2]

            fig = add_trace(fig, x = x, y = y, z = z, 
                type = "scatter3d",
                mode = "lines",
                line = list(color = 'rgb(0, 0, 0)', wdith=2), 
                hoverinfo="skip",
                showlegend=F
            )
        }
         
        ## add super cell
        fig = add_trace(fig,
            x = super_cell_d$lum_score,
            y = super_cell_d$bas_score,
            z = super_cell_d$PC_2,
            data = super_cell_d,
            type = 'scatter3d', mode = 'markers', 
            text = super_cell_d$CellType,
            hoverinfo = "text",
            color = super_cell_d$CellType,
            marker = list(
                size = 10, 
                line = list(color = 'rgba(0, 0, 0, 0.8)', width = 2) 
                ),
            showlegend = F
        )

        ## title & legends
        fig <- fig %>% layout(scene = list(xaxis = list(title = 'Luminal Score'),
                             yaxis = list(title = 'Basal Score'),
                             zaxis = list(title = 'PC 2')))
        fig

    })

    ## 2d trajectory plot
    select_pseudotime_gene = reactive({
        input_gene = input$gene
        validate(
            need(input_gene %in% colnames(exp_d), "The gene you searched is not in the dataset.")
        )
        input_gene
    })
    output$pseudotime = renderPlot({
        gene = select_pseudotime_gene()
        check_gene_both(gene, exp_d, meta_d, pseudotime_lum_d$Row.names, pseudotime_bas_d$Row.names, level_order, color_setting)
    })

})
