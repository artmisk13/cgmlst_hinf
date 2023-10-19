#!/usr/bin/env python3

import sys
import dendropy
import pandas as pd

outfile = sys.argv[2]

if __name__ == "__main__":
    t = sys.argv[1]
    tree = dendropy.Tree.get(path=t,
                             schema='newick',
                             preserve_underscores=True)

    d = {}
    pdm = tree.phylogenetic_distance_matrix()
    for taxon1 in tree.taxon_namespace:
        d[taxon1.label] = d.get(taxon1.label, {})
        for taxon2 in tree.taxon_namespace:
            if taxon2.label not in d[taxon1.label].keys():
                d[taxon1.label][taxon2.label] = pdm.patristic_distance(taxon1, taxon2)
                
    m = pd.DataFrame(d)
    m = m.reindex(m.columns)
    m.to_csv(outfile,
             sep='\t')
