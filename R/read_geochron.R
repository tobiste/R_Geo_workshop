clean_header <- function(x) {
  require(stringr)
  x <- stringr::str_replace(x, "[ ]+", "_")
  x <- stringr::str_replace(x, pattern = "[(]", replacement = "")
  x <- stringr::str_replace(x, pattern = "[)]", replacement = "")
  x <- stringr::str_replace(x, " ", "_")
  x <- stringr::str_replace(x, " ", "_")
  x <- stringr::str_squish(x)
  return(x)
}


#' @title Read data from Geochron database
#'
#' @description Converts Excel spreadsheets from a Geochon download into a meta
#' data and isotope data table. \code{meta} contains all the sample information.
#' \code{isotopes} contains all the information of the measurements (i.e. isotopes, ages)
#'
#' @param path Path to the xls file.
#' @param method Character string defining the input type of the dataset. One of
#' "Ar-Ar", "Fission Track", "(U-Th)/He", and "U-Pb"
#' @return list
#' @export
#' @importFrom dplyr "%>%" filter select
#' @importFrom stringr str_replace
#' @importFrom readxl read_excel
read_geochron <- function(path, method = "U-Pb") {
  require(dplyr)
  require(readxl)

  meta.header <- NULL


  if (method == "Ar-Ar") {
  }
  if (method == "(U-Th)/He") {
  }
  if (method == "Fission Track") {
    details0 <- readxl::read_excel(path = path, sheet = "Details", col_names = FALSE)

    details <- data.frame(t(details0$`...2`))
    names(details) <- clean_header(details0$`...1`)

    meta0 <- readxl::read_excel(path = path, sheet = "Metadata", col_names = FALSE)
    meta <- data.frame(t(meta0$`...2`))
    names(meta) <- clean_header(meta0$`...1`)

    meta <- cbind(meta, details)
    ages <- readxl::read_excel(path = path, sheet = "Ages")

    return(list(meta = meta, tracks = ages))
  }
  if (method == "U-Pb") {
    # meta
    meta0 <- readxl::read_excel(
      path = path,
      sheet = 1,
      skip = 7
    )
    meta.header <- colnames(meta0)
    colnames(meta0) <- clean_header(colnames(meta0))

    meta <- meta0 %>%
      filter(!is.na(Country)) %>%
      filter(Unique_ID != "Unique ID") %>%
      mutate(
        Longitude = as.numeric(Longitude),
        Latitude = as.numeric(Latitude),
        Max_Age_Ma = as.numeric(Max_Age_Ma),
        Min_Age_Ma = as.numeric(Min_Age_Ma),
        Oldest_Frac._Date_Ma = as.numeric(Oldest_Frac._Date_Ma),
        Youngest_Frac._Date_Ma = as.numeric(Youngest_Frac._Date_Ma)
      )

    # isotopes
    isotopes.header <- t(meta0[2, ])
    rownames(isotopes.header) <- NULL
    isotopes.header[1] <- "Sample ID"
    isotopes.header[2] <- "Unique ID"
    isotopes.header <-
      stringr::str_replace(isotopes.header, "[ ]", "_")
    isotopes.header <-
      stringr::str_replace(isotopes.header, "/", ".")
    isotopes.header <-
      stringr::str_replace(isotopes.header, "[ ]", "_")

    isotopes0 <-
      readxl::read_excel(
        path = path,
        sheet = 1,
        skip = 7
      )
    colnames(isotopes0) <- isotopes.header


    isotopes <- isotopes0[FALSE, ]
    for (i in seq_along(isotopes0$Sample_ID)) {
      if (!is.na(isotopes0$Sample_ID[i]) &
          !is.na(isotopes0$Unique_ID[i])) {
        sample_id <- isotopes0$Sample_ID[i]
        unique_id <- isotopes0$Unique_ID[i]
      }
      if (is.na(isotopes0$Sample_ID[i]) &
          is.na(isotopes0$Unique_ID[i])) {
        isotopes.i <- isotopes0[i, ]
        isotopes.i$Sample_ID <- sample_id
        isotopes.i$Unique_ID <- unique_id

        isotopes <- rbind(isotopes, isotopes.i)
      }
    }
    isotopes <- isotopes %>%
      dplyr::select(na.omit(isotopes.header)) %>%
      mutate(
        `206.238_Age` = as.numeric(`206.238_Age`),
        `206.238_Age_Error` = as.numeric(`206.238_Age_Error`),
        `207.235_Age` = as.numeric(`207.235_Age`),
        `207.235_Age_Error` = as.numeric(`207.235_Age_Error`),
        `207.206_Age` = as.numeric(`207.206_Age`),
        `207.206_Age_Error` = as.numeric(`207.206_Age_Error`),
        `Pb*.Pbc` = as.numeric(`Pb*.Pbc`),
        Rho = as.numeric(Rho),
        Rho_Error = as.numeric(Rho_Error),
        `206.238` = as.numeric(`206.238`),
        `206.238_Error` = as.numeric(`206.238_Error`),
        `206.204` = as.numeric(`206.204`),
        `208.206` = as.numeric(`208.206`),
        conc_U = as.numeric(conc_U),
        `Th.U_samp` = as.numeric(`Th.U_samp`),
        `207.235_Age` = as.numeric(`207.235_Age`),
        `207.235_Age_Error` = as.numeric(`207.235_Age_Error`),
        `207.206_Age` = as.numeric(`207.206_Age`),
        `207.206_Age_Error` = as.numeric(`207.206_Age_Error`),
        `206.238xTh_Age` = as.numeric(`206.238xTh_Age`),
        `206.238xTh_Age_Error` = as.numeric(`206.238xTh_Age_Error`),
        `207.235xPa_Age` = as.numeric(`207.235xPa_Age`),
        `207.235xPa_Age_Error` = as.numeric(`207.235xPa_Age_Error`),
        `207.206xTh_Age` = as.numeric(`207.206xTh_Age`),
        `207.206xTh_Age_Error` = as.numeric(`207.206xTh_Age_Error`),
        `207.206xPa_Age` = as.numeric(`207.206xPa_Age`),
        `207.206xPa_Age_Error` = as.numeric(`207.206xPa_Age_Error`)
      )
    colnames(isotopes) <- isotopes.header <- c(
      "Sample_ID", "Unique_ID", "Fraction_ID",
      "t.Pb206U238", "st.Pb206U238",
      "t.Pb207U235", "st.Pb207U235",
      "t.Pb207Pb206", "st.Pb207Pb206",
      "PbPb.cor", "rho", "s.rho",
      "Pb206U238", "errPb206U238",
      "Pb206Pb204", "Pb208Pb206", "U", "ThU",
      "Age_206.238xTh", "Age_Error_206.238xTh",
      "Age_207.235xPa", "Age_Error_207.235xPar",
      "Age_207.206xTh", "Age_Error_207.206xTh",
      "Age_207.206xPa_Age", "Age_Error_207.206xPa"
    )


    return(list(meta = meta, isotopes = isotopes))
  }
}

geochron_to_dz <- function(x) {
  lst <- list()
  for (i in x$Sample_ID) {
    x.i <- subset(x, x$Sample_ID == i)
    ages <- c(x.i$t)
    lst[[i]] <- ages
  }

  class(lst) <- "detritals"
  return(lst)
}
