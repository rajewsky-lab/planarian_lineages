# Cell Type Atlas and Lineage Tree Reconstruction of Whole Adult Animals by Single Cell Sequencing


### Graph Abstraction

The following notebooks provide the analyses based on *partition-based graph abstraction (PAGA)* [(Wolf *et al.*, 2017)](https://doi.org/10.1101/208819):

- [*planaria*](graph_abstraction/planaria.ipynb): inferring the lineage tree of Planaria
- [*epidermal-lineage*](graph_abstraction/epidermal-lineage.ipynb): zooming into the pseudotime-reconstruction of the epidermal lineage
- [*subsampling*](graph_abstraction/subsampling.ipynb): robustness of the tree under subsampling
- [*wildtype*](graph_abstraction/wildtype.ipynb): robustness when using only wild-type samples
- [*vary-neighbors*](graph_abstraction/vary-neighbors.ipynb): robustness under variation of the number of neighbors of the single-cell graph

### Velocyto 

[*PlanariaVelocytoAnalysis*](PlanariaVelocytoAnalysis.ipynb) provides the analysis using the *velocyto.py* package [(Gioele La Manno *et. al.*, 2017)](https://doi.org/10.1101/206052).

### Monocle

[*run_monocle2*](run_monocle2.R) reproduces the analysis using Monocle 2 [(Qiu *et al.*, 2017)](https://doi.org/10.1038/nmeth.4402).
