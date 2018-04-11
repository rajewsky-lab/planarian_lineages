# Cell Type Atlas and Lineage Tree of A Whole Complex Animal by Single-Cell Transcriptomics

### PAGA

The following notebooks provide the analyses based on *partition-based graph abstraction (PAGA)* [(Wolf *et al.*, 2017)](https://doi.org/10.1101/208819):

- [*planaria*](paga/planaria.ipynb): inferring the lineage tree of Planaria
- [*epidermal-lineage*](paga/epidermal-lineage.ipynb): zooming into the pseudotime-reconstruction of the epidermal lineage
- [*subsampling*](paga/subsampling.ipynb): robustness of the tree under subsampling
- [*wildtype*](paga/wildtype.ipynb): robustness when using only wild-type samples
- [*vary-neighbors*](paga/vary-neighbors.ipynb): robustness under variation of the number of neighbors of the single-cell graph

Note: In some of these notebooks, the layout of the graphs differ from the one in the supplement of the paper. See the versions of the notebooks in this [commit](https://github.com/rajewsky-lab/planarian_lineages/tree/00a3d752869613e4facedb4c722af442dc52a0f5) to reproduce the exact same layout.

### Velocyto 

[*PlanariaVelocytoAnalysis*](PlanariaVelocytoAnalysis.ipynb) provides the analysis using the *velocyto.py* package [(Gioele La Manno *et. al.*, 2017)](https://doi.org/10.1101/206052).

### Monocle

[*run_monocle2*](run_monocle2.R) reproduces the analysis using Monocle 2 [(Qiu *et al.*, 2017)](https://doi.org/10.1038/nmeth.4402).
