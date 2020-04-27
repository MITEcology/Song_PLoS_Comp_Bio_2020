library(tidyverse)
library(factoextra)
library(ggfortify)
library(here)
library(patchwork)

load("data.RData")


# Separabilioty in environment-independent approach-----------------------------------------------------------
data_info_ori <- data_info %>%
  filter(type == 'ori') %>% 
  select(var1, var2, var3, interaction_type) %>% 
  na.omit()

res_pca_ori <- prcomp(select(data_info_ori, -interaction_type), scale = TRUE)
p1 <- fviz_pca_ind(res_pca_ori,
                      col.ind = data_info_ori$interaction_type,
                      geom.ind = "point",
                      addEllipses = TRUE, # Concentration ellipses
                      ellipse.level=0.68,
                      legend.title = "Groups",
                      repel = TRUE,
                      title = 'Environmental-independent approach'
)

# Separabilioty in environment-dependent approach-----------------------------------------------------------
data_info_ori <- data_info %>%
  filter(type == 'ori') %>% 
  select(Temp_var, var1, var2, var3, interaction_type) %>% 
  na.omit()

res_pca_ori <- prcomp(select(data_info_ori, -interaction_type), scale = TRUE)
p2 <- fviz_pca_ind(res_pca_ori,
                      col.ind = data_info_ori$interaction_type,
                      geom.ind = "point",
                      addEllipses = TRUE, # Concentration ellipses
                      ellipse.level=0.68,
                      legend.title = "Groups",
                      repel = TRUE,
                      title = 'Environmental-dependent approach'
)

# Scalability -------------------------------------------------------------
data_info_ori <- data_info %>%
  filter(type == 'ori') %>% 
  select(Temp_var, var1, var2, var3, web_type) %>% 
  na.omit()

res_pca_ori <- prcomp(select(data_info_ori, -web_type), scale = TRUE)
p3 <- fviz_pca_ind(res_pca_ori,
                      col.ind = data_info_ori$web_type,
                      geom.ind = "point",
                      addEllipses = TRUE, # Concentration ellipses
                      ellipse.level=0.68,
                      legend.title = "Groups",
                      repel = TRUE,
                      title = 'Scalability'
)


# Specificity -------------------------------------------------------------
set.seed(1010)
data_info_er <- data_info %>%
    filter(type == 'er') %>% 
    group_by(ID) %>% 
    sample_n(1) %>% 
    ungroup() %>% 
    select(Temp_var, var1, var2, var3, interaction_type) %>% 
    na.omit()
  
res_pca_er <- prcomp(select(data_info_er, -interaction_type), scale = TRUE)
p4 <-   fviz_pca_ind(res_pca_er,
                  col.ind = data_info_er$interaction_type,
                  geom.ind = "point",
                  addEllipses = TRUE, # Concentration ellipses
                  ellipse.level=0.68,
                  legend.title = "Groups",
                  repel = TRUE,
                  title = 'Specificity'
  )

(p1 + p2) / (p3 + p4) + plot_annotation(tag_levels = 'A') & theme(aspect.ratio = 1)
ggsave('figure.pdf', height = 8, width =10)
