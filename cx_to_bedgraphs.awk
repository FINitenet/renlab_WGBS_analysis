BEGIN {
    OFS = "\t"
    while ((getline line < sizes_file) > 0) {
        split(line, f, "\t")
        chrom_size[tolower(f[1])] = f[2]
    }
    close(sizes_file)
    contexts[1] = "CG"; contexts[2] = "CHG"; contexts[3] = "CHH"; contexts[4] = "mC"
    print "sample", "group", "chrom", "start", "end", "context", \
          "methylated", "unmethylated", "covered_sites", "beta" > window_file
}

function flush_window(    i, c, e, cov) {
    if (current_chrom == "") return
    e = current_start + window_size
    if (e > chrom_size[current_chrom]) e = chrom_size[current_chrom]
    for (i = 1; i <= 4; i++) {
        c = contexts[i]
        cov = methylated[c] + unmethylated[c]
        if (cov > 0)
            print sample, group, current_chrom, current_start, e, c, \
                  methylated[c], unmethylated[c], covered[c], methylated[c] / cov > window_file
        methylated[c] = 0; unmethylated[c] = 0; covered[c] = 0
    }
}

{
    chrom = tolower($1)
    pos0 = $2 - 1
    meth = $4 + 0
    unmeth = $5 + 0
    context = $6
    cov = meth + unmeth
    if (!(context == "CG" || context == "CHG" || context == "CHH") || cov == 0) next

    window_start = int(pos0 / window_size) * window_size
    if (chrom != current_chrom || window_start != current_start) {
        flush_window()
        current_chrom = chrom
        current_start = window_start
    }
    methylated[context] += meth; unmethylated[context] += unmeth; covered[context]++
    methylated["mC"] += meth; unmethylated["mC"] += unmeth; covered["mC"]++

    if (cov >= min_site_coverage)
        print chrom, pos0, pos0 + 1, meth / cov > (prefix "." context ".bedGraph")
}

END {
    flush_window()
    close(window_file)
    for (i = 1; i <= 3; i++) close(prefix "." contexts[i] ".bedGraph")
}
