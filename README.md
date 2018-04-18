# Cell Type Atlas and Lineage Tree of A Whole Complex Animal by Single-Cell Transcriptomics

### PAGA

The following notebooks provide the analyses based on *partition-based graph abstraction (PAGA)* [(Wolf *et al.*, 2017)](https://doi.org/10.1101/208819):

- [*planaria*](https://nbviewer.jupyter.org/github/rajewsky-lab/planarian_lineages/blob/master/paga/planaria.ipynb): inferring the lineage tree from the topology of the Planarian cell atlas
- [*preprocessing*](https://nbviewer.jupyter.org/github/rajewsky-lab/planarian_lineages/blob/master/paga/preprocessing.ipynb): same as *planaria*, but including all preprocessing steps
- [*epidermal-lineage*](https://nbviewer.jupyter.org/github/rajewsky-lab/planarian_lineages/blob/master/paga/epidermal-lineage.ipynb): zooming into the pseudotime-reconstruction of the epidermal lineage
- [*subsampling*](https://nbviewer.jupyter.org/github/rajewsky-lab/planarian_lineages/blob/master/paga/subsampling.ipynb): robustness of the tree under subsampling
- [*wildtype*](https://nbviewer.jupyter.org/github/rajewsky-lab/planarian_lineages/blob/master/paga/wildtype.ipynb): robustness when using only wild-type samples
- [*vary-neighbors*](https://nbviewer.jupyter.org/github/rajewsky-lab/planarian_lineages/blob/master/paga/vary-neighbors.ipynb): robustness under variation of the number of neighbors in the single-cell graph

### Velocyto 

[*PlanariaVelocytoAnalysis*](PlanariaVelocytoAnalysis.ipynb) provides the analysis using the *velocyto.py* package [(Gioele La Manno *et. al.*, 2017)](https://doi.org/10.1101/206052).

### Monocle

[*run_monocle2*](run_monocle2.R) reproduces the analysis using Monocle 2 [(Qiu *et al.*, 2017)](https://doi.org/10.1038/nmeth.4402).
