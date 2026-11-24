// A library handed straight to an entry point (run_reclassify.nf without --dedup,
// run_ensemble.nf) skips CAT_TE_LIBRARIES and CD_HIT_EST, which keep IDs unique; yet
// evidence, model votes and provenance are keyed by the ID before '#'. Stop at launch
// when two records share one, instead of letting them share evidence and votes.
// Repeated sequences under distinct IDs are fine.
def checkLibraryIds(path) {
    def seen = [] as Set
    def dup  = [] as LinkedHashSet
    file(path).withReader { r ->
        r.eachLine { l ->
            if (l.startsWith('>')) {
                def b = l.substring(1).trim().split(/\s+/)[0].split('#')[0]
                if (!seen.add(b)) dup << b
            }
        }
    }
    if (dup) {
        error "Library ${path}: ${dup.size()} ID(s) occur more than once before '#' " +
              "(e.g. ${dup.take(3).join(', ')}). Give every record a unique ID, " +
              "or run run_reclassify.nf with --dedup."
    }
}
