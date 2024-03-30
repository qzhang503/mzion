find_mzml_indexes <- function (spec, is_3d = TRUE)
{
  ## the first scan
  x <- spec[[1]]
  ids <- .Internal(strsplit(xml2::xml_attr(x, "id"), " ", fixed = TRUE, 
                            perl = FALSE, useBytes = FALSE))[[1]]
  ids <- .Internal(strsplit(ids, "=", fixed = TRUE, 
                            perl = FALSE, useBytes = FALSE))
  id_nms <- lapply(ids, `[[`, 1)
  idx_demux <- which(id_nms == "demux")
  is_demux <- if (length(idx_demux)) TRUE else FALSE
  
  if (is_3d) {
    idx_sc <- which(id_nms == "scan")
    if (!length(idx_osc <- which(id_nms == "originalScan"))) idx_osc <- idx_sc
    idx_sc_start <- idx_sc_end <- idx_osc
  } else {
    idx_sc_start <- which(id_nms == "scanStart") # 3
    idx_sc_end <- which(id_nms == "scanEnd") # 4
    idx_osc <- idx_sc <- idx_sc_start
    # c(idx_sc_start, idx_sc_end) # length(2)
  }
  
  xc <- xml2::xml_children(x)
  xcp_attrs <- xml2::xml_attrs(xc[which(xml2::xml_name(xc) == "cvParam")])
  xcp_names <- lapply(xcp_attrs, `[[`, "name")
  xcp_vals <- lapply(xcp_attrs, `[[`, "value")
  idx_mslev <- which(xcp_names == "ms level")
  idx_title <- which(xcp_names == "spectrum title")
  rm(list = c("x", "ids", "xc", "xcp_attrs", "xcp_names", "xcp_vals", "id_nms"))
  
  # MS2 indexes (no MS1 with DDA -> MS2 goes first)
  len <- length(spec)
  rng <- max(floor(len/2L), 1L):len
  allowance <- min(100L, len)
  count <- 0L
  
  if (!len) {
    stop("No spectrum data found.")
  }
  
  for (i in rng) {
    # i=46
    # i = 6000
    x <- spec[[i]]
    xc <- xml2::xml_children(x)
    
    if (xml2::xml_attr(xc[[idx_mslev]], "value") == "2") {
      if (is_3d) {
        idx_scan_lwr_2 <- grep("lowest observed m/z", xc)
        if (!length(idx_scan_lwr_2)) {
          warning("Fields of `lowest observed m/z` not found.")
          idx_scan_lwr_2 <- 8L
        }
      } else {
        idx_scan_lwr_2 <- grep("ion mobility lower limit", xc)
        if (!length(idx_scan_lwr_2)) {
          warning("Fields of `ion mobility lower limit` not found.")
          idx_scan_lwr_2 <- 8L
        }
      }
      
      if (is_3d) {
        idx_scan_upr_2 <- grep("highest observed m/z", xc)
        if (!length(idx_scan_upr_2)) {
          warning("Fields of `highest observed m/z` not found.")
          idx_scan_upr_2 <- 9L
        }
      } else {
        idx_scan_upr_2 <- grep("ion mobility upper limit", xc)
        if (!length(idx_scan_upr_2)) {
          warning("Fields of `ion mobility upper limit` not found.")
          idx_scan_upr_2 <- 9L
        }
      }
      
      idx_precursor_2 <- grep("precursorList", xc)
      if (!length(idx_precursor_2)) {
        warning("Fields of `precursorList` not found.")
        idx_precursor_2 <- if (is_3d) 12L else 11L
      }
      
      idx_scanList_2 <- grep("scanList", xc)
      if (!length(idx_scanList_2)) {
        warning("Fields of `scanList` not found.")
        idx_scanList_2 <- if (is_3d) 11L else 10L
      }
      
      idx_bin_2 <- grep("binaryDataArrayList", xc)
      if (!length(idx_bin_2)) {
        warning("Fields of `binaryDataArrayList` not found.")
        idx_bin_2 <- if (is_3d) 13L else 12L
      }
      
      scanList <- xml2::xml_children(xc[[idx_scanList_2]])
      
      # !!! a range if timsTOF
      idx_rt_2 <- which(xml2::xml_name(scanList) == "scan")
      if (length(idx_rt_2)) {
        idx_rt_2 <- idx_rt_2[[1]]
      } else {
        warning("Fields of `scan` not found.")
        idx_rt_2 <- 2L
      }
      
      scanList_ret <- xml2::xml_children(scanList[[idx_rt_2]])
      scanList_ret_attrs <- xml2::xml_attr(scanList_ret, "name")
      
      idx_scan_start_2 <- which(scanList_ret_attrs == "scan start time") # 1
      if (!length(idx_scan_start_2)) {
        warning("Fields of `scan start time` not found.")
        idx_scan_start_2 <- 1L
      }
      
      if (is_3d) {
        if (is_pwiz) {
          idx_ms2_reso <- which(scanList_ret_attrs == "mass resolving power")
          if (length(idx_ms2_reso)) {
            ms2_reso <- scanList_ret[[idx_ms2_reso]]
            ms2_reso <- as.integer(xml2::xml_attr(ms2_reso, "value"))
          } else {
            warning("Fields of `mass resolving power` not found.")
            ms2_reso <- 120000L
          }
        } else { # MSFragger timsTOF
          idx_ms2_reso <- which(scanList_ret_attrs == "inverse reduced ion mobility")
          if (length(idx_ms2_reso)) {
            ms2_reso <- scanList_ret[[idx_ms2_reso]]
            ms2_reso <- as.numeric(xml2::xml_attr(ms2_reso, "value"))
          } else {
            warning("Fields of `inverse reduced ion mobility` not found.")
            ms2_reso <- 1L
          }
        }
      } else {
        idx_ms2_reso <- which(scanList_ret_attrs == "inverse reduced ion mobility")
        if (length(idx_ms2_reso)) {
          ms2_reso <- scanList_ret[[idx_ms2_reso]]
          ms2_reso <- as.numeric(xml2::xml_attr(ms2_reso, "value"))
        } else {
          warning("Fields of `inverse reduced ion mobility` not found.")
          ms2_reso <- 1
        }
      }
      
      # entire MS2 is empty
      if (length(xc) < idx_precursor_2) {
        next
      }
      
      precursorList <- xml2::xml_children(xc[[idx_precursor_2]])
      precursor <- precursorList[[1]] # assume one precursor
      precursorc <- xml2::xml_children(precursor)
      
      idx_selectedIonList <- grep("selectedIonList", precursorc)
      if (!length(idx_selectedIonList)) {
        warning("Fields of `selectedIonList` not found.")
        idx_selectedIonList <- 2L
      }
      
      idx_isolationWindow <- grep("isolationWindow", precursorc)
      if (!length(idx_isolationWindow)) {
        warning("Fields of `isolationWindow` not found.")
        idx_isolationWindow <- 1L
      }
      
      isolationWindowc <- xml2::xml_children(precursorc[[idx_isolationWindow]])
      iso_nms <- lapply(isolationWindowc, function (x) xml2::xml_attr(x, "name"))
      
      idx_ctrmz <- which(iso_nms == "isolation window target m/z")
      if (!length(idx_ctrmz)) {
        warning("Fields of `isolation window target m/z` not found.")
        idx_ctrmz <- 1L
      }
      
      idx_lwrmz <- which(iso_nms == "isolation window lower offset")
      if (!length(idx_lwrmz)) {
        warning("Fields of `isolation window lower offset` not found.")
        idx_lwrmz <- 2L
      }
      
      idx_uprmz <- which(iso_nms == "isolation window upper offset")
      if (!length(idx_uprmz)) {
        warning("Fields of `isolation window upper offset` not found.")
        idx_uprmz <- 3L
      }
      
      lwr_offset <- 
        as.numeric(xml2::xml_attr(isolationWindowc[[idx_lwrmz]], "value"))
      upr_off_set <- 
        as.numeric(xml2::xml_attr(isolationWindowc[[idx_uprmz]], "value"))
      iso_width <- upr_off_set + lwr_offset
      
      selectedIon <- xml2::xml_child(precursorc[[idx_selectedIonList]], 1)
      selectedIonc <- xml2::xml_children(selectedIon)
      selion_nms <- lapply(selectedIonc, function (x) xml2::xml_attr(x, "name"))
      
      idx_moverz <- which(selion_nms == "selected ion m/z")
      if (!length(idx_moverz)) {
        warning("Fields of `selected ion m/z` not found.")
        idx_moverz <- 1L
      }
      
      # not with timsTOF?
      idx_ms1int <- which(selion_nms == "peak intensity") # DIA: zero intensity; timsTOF: length(0)
      if (!length(idx_ms1int)) {
        if (count <= allowance) {
          count <- count + 1L
          next
        }
        else {
          warning("The field of `peak intensity` not found.")
          idx_ms1int <- 3L
        }
      }
      
      idx_charge <- which(selion_nms == "charge state") # DIA: no "charge state"
      is_dia <- if (length(idx_charge)) FALSE else TRUE
      
      # timsTOF
      if (is_dia && !is_3d) is_dia <- FALSE
      
      rm(list = c("scanList", "scanList_ret", "scanList_ret_attrs", 
                  "precursorList", "precursor", "precursorc", 
                  "selectedIon", "selectedIonc", "selion_nms"))
      break
    }
  }
  
  # MS1 indexes
  for (i in rng) {
    x <- spec[[i]]
    xc <- xml2::xml_children(x)
    
    if (as.integer(xml2::xml_attr(xc[[idx_mslev]], "value")) == 1) {
      if (is_3d) {
        idx_scan_lwr_1 <- grep("lowest observed m/z", xc) # "ion mobility lower limit"
        if (!length(idx_scan_lwr_1)) {
          warning("Fields of `lowest observed m/z` not found.")
          idx_scan_lwr_1 <- 8L
        }
        
        idx_scan_upr_1 <- grep("highest observed m/z", xc) # "ion mobility upper limit"
        if (!length(idx_scan_upr_1)) {
          warning("Fields of `highest observed m/z` not found.")
          idx_scan_upr_1 <- 9L
        }
      } else {
        idx_scan_lwr_1 <- grep("ion mobility lower limit", xc) # "ion mobility lower limit"
        if (!length(idx_scan_lwr_1)) {
          warning("Fields of `ion mobility lower limit` not found.")
          idx_scan_lwr_1 <- 8L
        }
        
        idx_scan_upr_1 <- grep("ion mobility upper limit", xc) # "ion mobility upper limit"
        if (!length(idx_scan_upr_1)) {
          warning("Fields of `ion mobility upper limit` not found.")
          idx_scan_upr_1 <- 9L
        }
      }
      
      idx_scanList_1 <- grep("scanList", xc)
      if (!length(idx_scanList_1)) {
        warning("Fields of `scanList` not found.")
        idx_scanList_1 <- if (is_3d) 11L else 10L
      }
      
      idx_bin_1 <- grep("binaryDataArrayList", xc)
      if (!length(idx_bin_1)) {
        warning("Fields of `binaryDataArrayList` not found.")
        idx_bin_1 <- if (is_3d) 12L else 11L
      }
      
      local({
        # [[1]] name="m/z array"; [[2]] name="intensity array"; [[3]] name="mean inverse reduced ion mobility array"
        binData <- xml2::xml_children(xc[[idx_bin_1]])[[1]]
        binDataC <- xml2::xml_children(binData)
        binDataC <- binDataC[which(xml2::xml_name(binDataC) == "cvParam")]
        
        lapply(binDataC, function (v) {
          attrs <- xml2::xml_attrs(v)
          oks <- grepl("compression", attrs)
          
          if (any(oks))
            if (attrs[oks][[1]] == "zlib compression")
              stop("Please uncheck zlib compression in  mzML generations.")
        })
      })
      
      scanList <- xml2::xml_children(xc[[idx_scanList_1]])
      idx_rt_1 <- which(xml2::xml_name(scanList) == "scan")
      if (!length(idx_rt_1)) {
        warning("Fields of `scan` not found.")
        idx_rt_1 <- 2L
      }
      
      scanList_ret <- xml2::xml_children(scanList[[idx_rt_1]])
      
      scan_ret_attrs <- xml2::xml_attr(scanList_ret, "name")
      ## for later removal of the assumption of a fixed resolution with Astral
      # ms1_reso <- scanList_ret[[which(scan_ret_attrs == "mass resolving power")]]
      # ms1_reso <- as.integer(xml2::xml_attr(ms1_reso, "value"))
      idx_scan_start_1 <- which(scan_ret_attrs == "scan start time")
      if (!length(idx_scan_start_1)) {
        warning("Fields of `scan start time` not found.")
        idx_scan_start_1 <- 1L
      }
      
      # idx_scan_start_im <- which(scan_ret_attrs == "inverse reduced ion mobility")
      
      rm(list = c("scanList", "scanList_ret", "scan_ret_attrs"))
      break
    }
  }
  
  rm(list = c("x", "xc", "i"))
  
  if (is_dia) {
    filename <- extrDIA(
      spec = spec, raw_file = raw_file, 
      temp_dir = temp_dir, is_demux = is_demux, 
      idx_sc = idx_sc, idx_osc = idx_osc, idx_mslev = idx_mslev, 
      idx_title = idx_title, idx_scanList_1 = idx_scanList_1, 
      idx_scanList_2 = idx_scanList_2, idx_rt_1 = idx_rt_1, 
      idx_rt_2 = idx_rt_2, idx_scan_start_1 = idx_scan_start_1, 
      idx_scan_start_2 = idx_scan_start_2, 
      idx_precursor_2 = idx_precursor_2, 
      idx_isolationWindow = idx_isolationWindow, 
      idx_ctrmz = idx_ctrmz, idx_lwrmz = idx_lwrmz, 
      idx_uprmz = idx_uprmz, 
      idx_selectedIonList = idx_selectedIonList, 
      idx_demux = idx_demux, idx_bin_1 = idx_bin_1, 
      idx_bin_2 = idx_bin_2, 
      idx_scan_lwr_1 = idx_scan_lwr_1, 
      idx_scan_upr_1 = idx_scan_upr_1, 
      idx_scan_lwr_2 = idx_scan_lwr_2, 
      idx_scan_upr_2 = idx_scan_upr_2)
  }
  else {
    filename <- extrDDA(
      spec = spec, raw_file = raw_file, temp_dir = temp_dir, 
      idx_sc = idx_sc, idx_osc = idx_osc, 
      idx_mslev = idx_mslev, idx_title = idx_title, 
      idx_scanList_2 = idx_scanList_2, idx_rt_2 = idx_rt_2, 
      idx_scan_start_2 = idx_scan_start_2, 
      idx_precursor_2 = idx_precursor_2, 
      idx_isolationWindow = idx_isolationWindow, 
      idx_ctrmz = idx_ctrmz, idx_lwrmz = idx_lwrmz,
      idx_uprmz = idx_uprmz, 
      idx_selectedIonList = idx_selectedIonList, 
      idx_moverz = idx_moverz, idx_charge = idx_charge, 
      idx_ms1int = idx_ms1int, idx_bin_2 = idx_bin_2, 
      idx_scanList_1 = idx_scanList_1, 
      idx_rt_1 = idx_rt_1, idx_scan_start_1 = idx_scan_start_1, 
      idx_bin_1 = idx_bin_1, 
      idx_scan_lwr_1 = idx_scan_lwr_1, idx_scan_upr_1 = idx_scan_upr_1, 
      idx_scan_lwr_2 = idx_scan_lwr_2, idx_scan_upr_2 = idx_scan_upr_2)
  }
  
  attr(filename, "is_dia") <- is_dia
  attr(filename, "iso_width") <- iso_width
  attr(filename, "mzml_type") <- mzml_type
  
  invisible(filename)
}