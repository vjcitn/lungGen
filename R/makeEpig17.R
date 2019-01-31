#' collect epigenetic states lifted to hg38, chr17 only!
#' @export
makeEpig17 = function() {
    ermaset = GenomicFiles::GenomicFiles(files = dir(system.file("bed", 
        package = "lungGen"), full.names = TRUE, pattern = "gz$"))
    data(epigColdata, package="lungGen")
    rownames(epigColdata) = epigColdata[,2]
    ftags = gsub("_.*", "", basename(files(ermaset)))
    ermaset@colData = epigColdata[ftags, ]
    genome(ermaset) = "GRCh38"
    ermaset # return plain old GenomicFiles
}

