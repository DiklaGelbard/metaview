# library(metacell)
# scdb_init("db/",force_reinit = T)
# default_mc_id <- "lung_kinetics_sorted"
# default_mc2d_id <- "lung_kinetics_sorted"
# default_mat_id <- "lung_kinetics"
# current_mc_id <- paste0(default_mc_id,"_new")
# current_mc2d_id <- paste0(default_mc2d_id,"_new")

default_mc_col = read.delim("db/color_schemes/lung_mc_colorize.txt")

df_2d = read.delim("db/df_2d.txt")
df_mc_2d =read.delim("db/df_mc_2d.txt")

mc_mc = read.delim("db/mc_mc.txt")[[1]]
names(mc_mc) = rownames(read.delim("db/mc_mc.txt"))
marks_genes = read.delim("db/marks_genes.txt",stringsAsFactors = FALSE)[[1]]
mc_fp = as.matrix(read.delim("db/mc_fp.txt"))
colnames(mc_fp) = 1:ncol(mc_fp)
mc_genes = sort(rownames(mc_fp))
mc_colors = setNames(read.delim("db/mc_colors.txt",stringsAsFactors = FALSE)[[1]],colnames(mc_fp)) 

# mc = scdb_mc(default_mc_id)
# mc_genes = sort(rownames(mc@mc_fp))
# mc_cols = mc@colors
# mc_cols[is.na(mc@colors)] = "gray"
# mat = scdb_mat(default_mat_id)
# mc2d = scdb_mc2d(default_mc2d_id)


