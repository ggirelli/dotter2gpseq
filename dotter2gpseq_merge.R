#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 2.1.0
# Description:
# 	merge output of dotter2gpseq.py and add dataset and cell_type information.
# 
# Changelog:
#  v2.1.0 180209 - changed terminology in output, copy instead of allele.
#  v2.0.1 - fixed probe name assignment.
#  v2.0.0 - requires only one metadata table.
#  v1.1.0 - now supporting set/probe label columns in metadata.
#  v1.0.0 - first implementation.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

suppressMessages(library(argparser))

# INPUT ========================================================================

# Create arguent parser
parser = arg_parser('Description: merge dotter2gpseq.py output, add dataset and
cell type information. The software looks for dotter2gpseq.py output in
subfolders of the specified input directories. These subfolders must be named as
dataset_series, with series being in XXX format with leading zeros. Please, note
that the channel is enforced as lower-case by the merge operation.

Example 1: output in current directory.
./merge_data.R -m meta.tsv -i HAP1/dots_auto IMR90/dots_auto

Example 2: output to "/home/user/out" directory.
./merge_data.R -m meta.tsv -i HAP1/dots_auto IMR90/dots_auto
	-o /home/user/out',
	name = 'dotter2gpseq_merge.R'
)

# Define mandatory arguments
parser = add_argument(parser, arg = '--meta', short = '-m', nargs = 1,
	help = paste0('Metadata table.',
		' Needed columns: dataset, series, cell_line, set_label, probe_label.'))
parser = add_argument(parser, arg = '--indir', short = '-i', nargs = Inf,
	help = 'List of input folders with dotter2gpseq output as subfolders.')
parser = add_argument(parser, arg = '--outdir', short = '-o', nargs = 1,
	help = 'Output folder, created if missing. Default to current one.',
	default = ".")
parser = add_argument(parser, arg = '--aspect', short = '-a', nargs = 3,
	help = paste0('Physical size of Z, Y and X voxel sides.',
		' Default: 300.0 130.0 130.0'), default = c(300.0, 130.0, 130.0))

# Parse arguments
p = parse_args(parser)

# Attach argument values to variables
attach(p['' != names(p)])

# Additional checks
if( is.na(meta) ) {
	stop("at least one matadata table must be provided.")
}
if( 1 == length(indir) ) { if( is.na(indir) ) {
	stop("at least one input directory must be provided.")
}}
if( is.na(outdir) ) {
	stop("at least one output folder must be provided.")
}
# RUN ==========================================================================

# Read metadata table
md = read.delim(meta, as.is = T, header = T)

# Check that all columns are present
l0 = lapply(c("dataset", "series", "cell_line", "set_label", "probe_label"),
	FUN = function(x) {
	if( !x %in% colnames(md) ) {
		stop(paste0("metadata missing the '", x, "' column.",
			"\nFile at: ", meta))
	}
})

# Iterate through metadata
l2 = by(md, paste0(md$dataset, "~", md$series), FUN = function(subt) {

	# Extract dataset and series information -------------------------------
	dataset = subt$dataset[1]
	series = sprintf("%03d", subt$series[1])
	flag = paste0(dataset, "_", series)
	cell_type = subt$cell_line[1]
	label = subt$set_label[1]

	# Log current status
	cat(paste0("Working on ", flag, ".\n"))

	# Identify input folder and skip if missing ----------------------------
	ipaths = paste0(indir, "/", flag)
	l = lapply(ipaths, FUN = function(ipath) {
		out = list()
		out$good = T
		out$partial = F

		if ( !dir.exists(ipath) ) {
			out$error = paste0("Warning: cannot find folder for ", flag,
				". Skipped.\nFolder at: ", ipath, "\n")
			out$good = F
			return(out)
		}

		# Identify input files -------------------------------------------------
		flist = list.files(ipath)
		out$nuclei = flist[grepl("nuclei.out", flist)]
		out$dots = flist[
			grepl("wCentr.out", flist) & ! grepl("noAllele", flist)]

		if( 0 == length(out$nuclei) ) {
			out$error = paste0("Warning: cannot find nuclei ",
				"information in ", flag, ". Skipping ", dataset, ".\n",
				"Folder: ", ipath)
			out$good = F
			out$partial = T
		} else {
			out$nuclei = read.delim(paste0(ipath, "/", out$nuclei),
				as.is = T, header = T)
		}

		if( 0 == length(out$dots) ) {
			out$error = paste0("Warning: cannot find dot information in ",
				flag, ". Skipping ", dataset, ".\n",
				"Folder: ", ipath)
			out$good = F
			out$partial = T
		} else {
			out$dots = read.delim(paste0(ipath, "/", out$dots),
				as.is = T, header = T)
			out$dots$Channel = tolower(out$dots$Channel)
		}

		return(out)

	})

	good = unlist(lapply(l, FUN = function(x) { x$good }))
	if ( !any(good) ) {
		partial = unlist(lapply(l, FUN = function(x) { x$partial }))
		if ( any(partial) ) {
			# Print error due to partial information present
			cat(l[[which(partial)[1]]]$error)
		} else {
			# Print error due to not-found information
			cat(paste0("Warning: cannot find information on dataset ",
				dataset, " in any of the input directories.\n"))
		}
		# Stop evaluating current dataset
		return(NULL)
	}

	lid = which(good)[1]
	dots = l[[lid]]$dots
	nuclei = l[[lid]]$nuclei
	print(l)

	# Add dataset, series and cell_type information ------------------------
	dots$dataset = rep(dataset, nrow(dots))
	dots$label = rep(label, nrow(dots))
	dots$probe_label = subt$probe_label[match(dots$Channel, subt$channel)]
	dots$cell_type = rep(cell_type, nrow(dots))
	nuclei$dataset = rep(dataset, nrow(nuclei))
	nuclei$cell_type = rep(cell_type, nrow(nuclei))

	# Prepare allele by channel table --------------------------------------
	aldata = as.numeric(dots$Allele)
	aldata = dots[0 < aldata & !is.na(aldata),]
	if( 0 != nrow(aldata) ) {
		aldata$universal = paste(
			aldata$File, aldata$Channel, aldata$cell_ID, sep = "~")
		alleles = do.call(rbind, by(aldata, aldata$universal,
			FUN = function(subt) {
			d = subt[1, c("File", "Channel", "cell_ID", "G1")]
			d_3d = subt[1, c("x", "y", "z")] - subt[2, c("x", "y", "z")]
			d$d_3d = sqrt(sum(((d_3d) * aspect)^2))
			d$d_lamin = abs(diff(subt$lamin_dist))
			d$d_lamin_norm = abs(diff(subt$lamin_dist_norm))
			d$d_centr = abs(diff(subt$centr_dist))
			d$d_centr_norm = abs(diff(subt$centr_dist_norm))
			d$angle = subt$angle[1]
			d$dataset = dataset
			d$label = label
			d$probe_label = subt$probe_label[1]
			d$cell_type = cell_type
			return(d)
		}))
		rownames(alleles) = c()
	} else {
		cat(paste0("Warning: no allele couples found in ", flag, ".\n"))
		alleles = NULL
	}

	# Output ---------------------------------------------------------------
	return(list(dots = dots, nuclei = nuclei, alleles = alleles))
})

# Remove skipped
l2 = l2[!is.null(l2)]

# Output
alleles = lapply(l2, FUN = function(x) x[[3]])
dots = do.call(rbind, lapply(l2, FUN = function(x) x[[1]]))
nuclei = do.call(rbind, lapply(l2, FUN = function(x) x[[2]]))
alleles = do.call(rbind, alleles[!is.null(alleles)])

# Write output
if( !dir.exists(outdir) ) dir.create(outdir)
write.table(dots, paste0(outdir, "/", Sys.Date(), "_dots.merged.tsv"),
	col.names = T, row.names = F, quote = F, sep = "\t")
write.table(alleles, paste0(outdir, "/", Sys.Date(), "_copies.merged.tsv"),
	col.names = T, row.names = F, quote = F, sep = "\t")
write.table(nuclei, paste0(outdir, "/", Sys.Date(), "_nuclei.merged.tsv"),
	col.names = T, row.names = F, quote = F, sep = "\t")

# END --------------------------------------------------------------------------

################################################################################
