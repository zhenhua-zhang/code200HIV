# Create date: 2019-Nov-5
# Last update: 2019-Dec-5
# Version    : 0.1.0
# Author     : Zhenhua Zhang
# Email      : zhenhua.zhang217@gmail.com

# NOTE:
#     1. In the future, this script will be implemented into a packages called
#     `zzhmisc`, which is the short for  ZhangZhenHuaMisc

# Do not load the loaded again
if (! "data.table" %in% (.packages())) {
    library(data.table)
}

# A smarter read function
smtread <- function(file_path, ..., idxc = "id", kpr = NULL, rmr = NULL, kpc = NULL, rmc = NULL,
                    trps = FALSE) {
    # A function to read your file smarter
    # 1. Remove / choose columns and rows after loading the file.
    # 2. Transform the data.frame after loading the file.
    # 3. It does remove and choose first, and then do the transform, if all
    #    of them are requested.
    # 4. You also would like to give the id columns, which will help make life
    #    easier

    dtfm <- fread(
        file_path,
        data.table = FALSE, stringsAsFactors = FALSE, verbose = FALSE, ...
    )

    col_names <- colnames(dtfm)
    if (idxc %in% col_names) {
        rownames(dtfm) <- dtfm[, idxc]
        dtfm <- dtfm[, col_names[!col_names %in% c(idxc)]]
    } else {
        stop("The given `idxc = ", idxc, "` is not in the column names")
    }

    if (!is.null(kpr)) {
        kpt_rows <- kpr
    } else {
        kpt_rows <- rownames(dtfm)
    }

    if (!is.null(rmr)) {
        kpt_rows <- kpt_rows[!kpt_rows %in% rmr]
    }

    if (!is.null(kpc)) {
        kpt_cols <- kpc
    } else {
        kpt_cols <- colnames(dtfm)
    }

    if (!is.null(rmc)) {
        kpt_cols <- kpt_cols[!kpt_cols %in% rmc]
    }

    if (!all(rownames(dtfm) %in% kpt_rows)) {
        if (length(kpt_rows) > 1) {
            dtfm <- dtfm[kpt_rows, ]
        } else {
            singlten <- as.list(dtfm[kpt_rows, ])
            names(singlten) <- colnames(dtfm)
            dtfm <- data.frame(singlten, row.names = kpt_rows)
        }
    }

    if (!all(colnames(dtfm) %in% kpt_cols)) {
        if (length(kpt_cols) > 1) {
            dtfm <- dtfm[, kpt_cols]
        } else {
            singlten <- dtfm[, kpt_cols]
            dtfm <- data.frame(singlten, row.names = rownames(dtfm))
            colnames(dtfm) <- kpt_cols
        }
    }

    if (trps) {
        # Any is character
        anic <- any(sapply(dtfm, function(e) {
            return(typeof(e) == "character")
        }))

        if (anic) {
            stop("There's character in the data.frame, plase remove them then try transform again")
        }
        dtfm <- as.data.frame(t(dtfm))
    }

    return(dtfm)
}

lambda <- function(stvc, sttp="PVAL", rnd=3) {
    if(sttp == "Z") {
        z = stvc
    } else if (sttp == "CHISQ") {
        z = sqrt(stvc)
    } else if (sttp == "PVAL") {
        z = qnorm(stvc / 2)
    }

    return(round(median(z^2) / qchisq(0.5, 1), rnd))
}
