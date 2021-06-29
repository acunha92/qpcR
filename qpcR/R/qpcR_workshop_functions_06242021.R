##Setup working directory


##Alexander Cunha
##qpcR Functions
##6/29/2021

##Summary function
##file = .csv raw data cols[Gene, Sample, Ct, Replicate]
pcR.summary <- function(file,
                        replicates) {data.og <- read.csv(file = file,
                                                         header = T,
                                                         stringsAsFactors = F);
                        data.og$Ct <- as.numeric(data.og$Ct);
                        data.summary <- tidyr::pivot_wider(data.og,
                                                           id_cols = c(Sample, Gene),
                                                           names_from = Replicate,
                                                           values_from = Ct);
                        data.mutate <- dplyr::mutate(data.summary,
                                                     Mean.Ct = rowMeans(data.summary[,(3:(replicates+2))]),
                                                     StdDev = apply(data.summary[,(3:(replicates+2))], 1, sd)) |> dplyr::arrange(Sample, Gene);
                        write.csv(data.mutate, file = paste0("qpcR_summary_", Sys.Date(),".csv"), row.names = F)}

##Plot Mean Ct Values
##summary.file = .csv produced from pcR.summary function, .csv only
pcR.meanplot <- function(summary.file,
                         title = as.character(summary.file)) {summary.data <- read.csv(file = summary.file,
                                                                                  header = T,
                                                                                  stringsAsFactors = F);
                         mean.plot <- ggplot2::ggplot(data = summary.data,
                                                      ggplot2::aes(x = Sample,
                                                                   y = Mean.Ct,
                                                                   fill = Gene)) +
                           ggplot2::geom_bar(stat = "identity",
                                             position = "dodge",
                                             color = "black",
                                             size = 0.9) +
                           ggplot2::geom_errorbar(ggplot2::aes(ymin = Mean.Ct - StdDev,
                                                               ymax = Mean.Ct + StdDev),
                                                  size = 0.9,
                                                  width = 0.2,
                                                  position = ggplot2::position_dodge(0.9)) +
                           ggplot2::ylab("Mean Ct Value") +
                           ggplot2::labs(title = title) +
                           ggplot2::theme_bw() +
                           ggplot2::scale_y_continuous(limits = c(0,40), n.breaks = 4) +
                           ggplot2::scale_fill_brewer(palette = 'Set2') +
                           ggplot2::theme(axis.text = ggplot2::element_text(size = "12", face = "bold", color = "black"),
                                          axis.title = ggplot2::element_text(size = "18", face = "bold", color = "black"));
                         ggplot2::ggsave(plot = mean.plot,
                                         file = paste0("meanCTplot_",Sys.Date(),".pdf"),
                                         device = "pdf",
                                         width = 14,
                                         height = 8.5,
                                         units = "in")}

##DeltaCt and Relative Expression Summary
##summary.file = .csv produced from pcR.summary function, .csv only
##goi = list of genes of interest
##hkg = string of single housekeeping gene to compare to, output table specifies hkg of interest
# pcR.deltaCT <- function(summary.file,
#                         goi,
#                         hkg) {summary.data <- read.csv(file = summary.file,
#                                                        header = T,
#                                                        stringsAsFactors = F);
#                         data.pivot <- tidyr::pivot_wider(summary.data,
#                                                          id_cols = Sample,
#                                                          names_from = Gene,
#                                                          values_from = Mean.Ct);
#                         cross <- purrr::cross2(.x = goi,
#                                                .y = hkg);
#                         cross.frame <- as.data.frame(do.call(rbind, cross));
#                         delta.ct <- purrr::map2_df(.x = cross.frame[,1],
#                                                    .y = cross.frame[,2],
#                                                    .f = ~dplyr::mutate(data.pivot,
#                                                                        deltaCT_ = data.pivot[.x] - data.pivot[.y]));
#                         delta.ct$HKG <- paste(hkg);
#                         write.csv(delta.ct,
#                                   file = paste0('deltaCT_',hkg,'_',Sys.Date(),'.csv'),
#                                   row.names = F)}

##Relative Expression Added to DeltaCt?
pcR.expression <- function(summary.file,
                                goi,
                                hkg) {summary.data <- read.csv(file = summary.file,
                                                               header = T,
                                                               stringsAsFactors = F);
                                data.pivot <- tidyr::pivot_wider(summary.data,
                                                                 id_cols = Sample,
                                                                 names_from = Gene,
                                                                 values_from = Mean.Ct);
                                cross <- purrr::cross2(.x = goi,
                                                       .y = hkg);
                                cross.frame <- as.data.frame(do.call(rbind, cross));
                                delta.ct <- purrr::map2_df(.x = cross.frame[,1],
                                                           .y = cross.frame[,2],
                                                           .f = ~dplyr::mutate(data.pivot,
                                                                               deltaCT = data.pivot[.x] - data.pivot[.y])) |> dplyr::mutate(RelExpr = 2^-(deltaCT),
                                                                                                                                            HKG = paste(hkg)) |> dplyr::select(-all_of(goi), -all_of(hkg)) |> as.matrix();
                                utils::write.csv(x = delta.ct,
                                                 file = paste0('qpcR_expression_',as.character(hkg),'_',Sys.Date(),'.csv'),
                                                 row.names = F)}


##Plot Relative Expression
##expression.file = .csv produced from pcR.expression function, .csv only
pcR.expressionplot <- function(expression.file,
                               title = as.character(expression.file)) {expression.data <- read.csv(file = expression.file,
                                                                                                   header = T,
                                                                                                   stringsAsFactors = F);
                               expression.melt <- expression.data |> dplyr::select(Sample, starts_with('RelExpr')) |> reshape::melt();
                               expression.plot <- ggplot2::ggplot(data = expression.melt,
                                                                  ggplot2::aes(x = Sample,
                                                                               y = value,
                                                                               fill = variable)) +
                                 ggplot2::geom_bar(stat = "identity",
                                                   position = "dodge",
                                                   color = "black",
                                                   size = 0.9) +
                                 ggplot2::geom_abline(slope = 0, intercept = 1, linetype = 'dashed', color = 'darkred') +
                                 ggplot2::ylab(paste0("Relative Expression to ",expression.data$HKG)) +
                                 ggplot2::labs(title = title) +
                                 ggplot2::theme_classic() +
                                 ggplot2::scale_fill_brewer(palette = 'Set2') +
                                 ggplot2::theme(axis.text = ggplot2::element_text(size = "12", face = "bold", color = "black"),
                                                axis.title = ggplot2::element_text(size = "18", face = "bold", color = "black"));
                               ggplot2::ggsave(plot = expression.plot,
                                               file = paste0("pcRexpressionplot_",Sys.Date(),".pdf"),
                                               device = "pdf",
                                               width = 14,
                                               height = 8.5,
                                               units = "in")}
