#!/bin/awk
BEGIN{
    print "Job started"
    THRE = 5e-2
}
{
    if (NR == 1) {
        COLNAMES = $1"\t"$5
    } else {
        if ($5 <= THRE) {
            OUTPUT_FILE = $2"_qtls_"THRE".txt"
            print $1"\t"$5 >> OUTPUT_FILE
        }
        close(OUTPUT_FILE)
    }
}
END{print "Job was done!"}
