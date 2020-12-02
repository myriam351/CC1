R Notebook
================

  - [Workflow for Microbiome Data Analysis: from raw reads to community
    analyses -
    PhyloSeq](#workflow-for-microbiome-data-analysis-from-raw-reads-to-community-analyses---phyloseq)
      - [Bonus: Handoff to phyloseq](#bonus-handoff-to-phyloseq)
      - [Ordination](#ordination)
      - [Bar plot](#bar-plot)
      - [Using phyloseq](#using-phyloseq)
      - [Prevalence Filtering](#prevalence-filtering)
      - [Agglomerate taxa](#agglomerate-taxa)
      - [Abundance value
        transformation](#abundance-value-transformation)
      - [Subset by taxonomy](#subset-by-taxonomy)
      - [Preprocessing](#preprocessing)
      - [Different Ordination
        Projections](#different-ordination-projections)
      - [Why are the ordination plots so far from
        square?](#why-are-the-ordination-plots-so-far-from-square)
      - [Canonical correspondence](#canonical-correspondence)
      - [Supervised learning](#supervised-learning)
      - [Graph-based analyses](#graph-based-analyses)
      - [Graph-based two-sample tests](#graph-based-two-sample-tests)
      - [Minimum Spanning Tree (MST)](#minimum-spanning-tree-mst)
      - [Nearest neighbors](#nearest-neighbors)
      - [Linear modeling](#linear-modeling)
      - [Hierarchical multiple testing](#hierarchical-multiple-testing)
      - [Conclusion :](#conclusion)

# Workflow for Microbiome Data Analysis: from raw reads to community analyses - PhyloSeq

## Bonus: Handoff to phyloseq

``` r
load("02_stat-analysis-with-DADA2_FinalENV")
```

### Commentaire : ici on a lié les données de Dada2 dans le fichier 02\_stat-analysis à ce fichier pour effectuer la suite des analyses.

``` r
library(phyloseq); packageVersion("phyloseq")
```

    ## [1] '1.34.0'

``` r
library(Biostrings); packageVersion("Biostrings")
```

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     distance

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## [1] '2.58.0'

``` r
library(ggplot2); packageVersion("ggplot2")
```

    ## [1] '3.3.2'

### Commentaire : Ici on a importer directement dans le Phyloseq, les tableaux produits par le pipeline Dada2.

``` r
theme_set(theme_bw())
```

``` r
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out
```

### Commentaire : Construction d’un échantillon data.frame à partir de données d’échantillons dans des fichiers.

``` r
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
```

### Commentaire : Construction d’un objet phyloseq à partir des sorties de Dada2.

``` r
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 232 taxa and 19 samples ]
    ## sample_data() Sample Data:       [ 19 samples by 4 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 232 taxa by 7 taxonomic ranks ]
    ## refseq()      DNAStringSet:      [ 232 reference sequences ]

### Commentaire : Utilisation de noms plus courts pour faire au plus simple. Par exemple: on a renommer nos taxons en chaîne plus courte. Ainsi, on aura des nouveaux noms plus courts qui apparaît dans les tableaux.

### Commentaire : Maintenant on peut utiliser PhyloSeq.

``` r
plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")
```

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->
\#\#\# Commentaire : Par exemple : on peut visualiser la diversité
alpha. \#\#\# Interprétation: On observe qu’il n’y a pas de différence
marquante en terme de diversité alpha entre les échantillons précoces et
tardis.

## Ordination

``` r
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
```

    ## Run 0 stress 0.08043117 
    ## Run 1 stress 0.08076346 
    ## ... Procrustes: rmse 0.01062291  max resid 0.03271618 
    ## Run 2 stress 0.09477284 
    ## Run 3 stress 0.08616061 
    ## Run 4 stress 0.089892 
    ## Run 5 stress 0.0898921 
    ## Run 6 stress 0.08076345 
    ## ... Procrustes: rmse 0.01061518  max resid 0.03269139 
    ## Run 7 stress 0.08616093 
    ## Run 8 stress 0.08988948 
    ## Run 9 stress 0.08616061 
    ## Run 10 stress 0.08616061 
    ## Run 11 stress 0.08043117 
    ## ... Procrustes: rmse 1.720233e-06  max resid 3.941785e-06 
    ## ... Similar to previous best
    ## Run 12 stress 0.1010632 
    ## Run 13 stress 0.08989187 
    ## Run 14 stress 0.08616063 
    ## Run 15 stress 0.08989265 
    ## Run 16 stress 0.08043116 
    ## ... New best solution
    ## ... Procrustes: rmse 4.169631e-06  max resid 1.056527e-05 
    ## ... Similar to previous best
    ## Run 17 stress 0.08043117 
    ## ... Procrustes: rmse 3.658694e-06  max resid 8.865144e-06 
    ## ... Similar to previous best
    ## Run 18 stress 0.1212047 
    ## Run 19 stress 0.1010633 
    ## Run 20 stress 0.08988953 
    ## *** Solution reached

``` r
plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")
```

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->
\#\#\# Commentaire : Transformation des données en proportions
appropriées pour les distances Bray-Curtis. On a ensuite fait de
l’ordination (NMDS). Comme on peut le voir, l’ordination a effectué
une séparation nette entre les échantillons précoces et tardifs.

## Bar plot

``` r
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")
```

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->
\#\#\# Commentaire : Rien de marquant ne ressort de la distribution
taxonomique des 20 premières séquences pour expliquer la différenciation
précoce-tardive.

## Using phyloseq

### Taxonomic Filtering

``` r
ps_connect<-url("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/ps.rds")
ps=readRDS(ps_connect)
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 389 taxa and 360 samples ]
    ## sample_data() Sample Data:       [ 360 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 389 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 389 tips and 387 internal nodes ]

### Commentaire : On récupéré les données des auteurs pour effectué la suite sur Phyloseq, étant donné qu’’il y a des différences entre les deux données (dont l’un qu’on a utilisé sur Dada2 et l’autre sur Phyloseq).

``` r
rank_names(ps)
```

    ## [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"

### Commentaire : Ici on a affiché les rangs présents dans l’ensemble des données.

``` r
table(tax_table(ps)[, "Phylum"], exclude = NULL)
```

    ## 
    ##              Actinobacteria               Bacteroidetes 
    ##                          13                          23 
    ## Candidatus_Saccharibacteria   Cyanobacteria/Chloroplast 
    ##                           1                           4 
    ##         Deinococcus-Thermus                  Firmicutes 
    ##                           1                         327 
    ##                Fusobacteria              Proteobacteria 
    ##                           1                          11 
    ##                 Tenericutes             Verrucomicrobia 
    ##                           1                           1 
    ##                        <NA> 
    ##                           6

### Commentaire : Ensuite, on a créer un tableau avec le nombre de caractéristiques pour chaque Phyla. Ainsi, on peut voir quelques Phyla qui présentent des caractéristiques.

``` r
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
```

### Commentaire: Ici on a fait en sorte que les caractéristiques avec une annotation ambiguë soient supprimées.

### Commentaire : Ensuite, on va explorer la prévalence des caractéristiques dans l’ensemble des données.

``` r
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
```

### Commentaire : Ici on a calculé la prévalence de chaque caractéristique, stocker sous forme de data.frame

``` r
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))
```

### Commentaire : Ensuite, on a ajouté la taxonomie et le nombre de lectures dans nos données.

``` r
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
```

    ##                         Phylum         1     2
    ## 1               Actinobacteria 120.15385  1562
    ## 2                Bacteroidetes 265.52174  6107
    ## 3  Candidatus_Saccharibacteria 280.00000   280
    ## 4    Cyanobacteria/Chloroplast  64.25000   257
    ## 5          Deinococcus-Thermus  52.00000    52
    ## 6                   Firmicutes 179.24771 58614
    ## 7                 Fusobacteria   2.00000     2
    ## 8               Proteobacteria  59.09091   650
    ## 9                  Tenericutes 234.00000   234
    ## 10             Verrucomicrobia 104.00000   104

### Commentaire : On a calculé ici la prévalence totale et moyenne des caractéristiques de chaque Phylum.

``` r
filterPhyla = c("Fusobacteria", "Deinococcus-Thermus")
```

### Commentaire : On a définit les phyla pour les filtrer par la suite.

``` r
ps1 = subset_taxa(ps, !Phylum %in% filterPhyla)
ps1
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 381 taxa and 360 samples ]
    ## sample_data() Sample Data:       [ 360 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 381 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 381 tips and 379 internal nodes ]

### Commentaire : On a Filtré les entrées avec les Phylum non identifiés. On voit bien que les données on bien été filtrées en se référant aux données non filtrées :

### otu\_table() OTU Table: \[ 232 taxa and 19 samples \]

### sample\_data() Sample Data: \[ 19 samples by 4 sample variables \]

### tax\_table() Taxonomy Table: \[ 232 taxa by 6 taxonomic ranks \]

### refseq() DNAStringSet: \[ 232 reference sequences \]

## Prevalence Filtering

``` r
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) + geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
```

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->
\#\#\# Commentaire : Ici on a fait un sous-ensemble du reste du Phylum
en incluant une estimation pour le paramètre. Comme on voit sur les
graphiques, on a en abscisse le total d’abondance (par comptage) et en
ordonnée on a la prévalence. Par exemple : on peut voir que pour les
Firmicutes, on une prévalence importante. Et au contraire, on peut voir
que pour les Tenericutes, on a une abondance très faible, quand on les
comparent.

``` r
prevalenceThreshold = 0.05 * nsamples(ps)
prevalenceThreshold
```

    ## [1] 18

### Commentaire : On a Définit le seuil de prévalence de 5 % du total des échantillons.

``` r
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps)
```

### Commentaire : On a exécuté le filtre de prévalence (soit de 5 %), en utilisant la fonction "prune\_taxa()\`.

## Agglomerate taxa

``` r
length(get_taxa_unique(ps2, taxonomic.rank = "Genus"))
```

    ## [1] 49

### Commentaire : Ici c’est un exemple de code, afin de voir combien de genres seraient présents après le filtrage. Le nombre de genres après filtrage serait de 49.

``` r
ps3 = tax_glom(ps2, "Genus", NArm = TRUE)
```

``` r
h1 = 0.4
ps4 = tip_glom(ps2, h = h1)
```

``` r
h1 = 0.4
ps4 = tip_glom(ps2, h = h1)
```

``` r
multiPlotTitleTextSize = 15
p2tree = plot_tree(ps2, method = "treeonly",
                   ladderize = "left",
                   title = "Before Agglomeration") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p3tree = plot_tree(ps3, method = "treeonly",
                   ladderize = "left", title = "By Genus") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p4tree = plot_tree(ps4, method = "treeonly",
                   ladderize = "left", title = "By Height") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
```

### Commentaire : On a comparé les données orginiales non filtrées avec la fonction plot\_tree().

``` r
library(gridExtra)
```

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     combine

``` r
grid.arrange(nrow = 1, p2tree, p3tree, p4tree)
```

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->
\#\#\# Commentaire : Après comparaison des données, elles sont stockées
sous forme d’objets de parcelle séparés (dont p2tree, p3tree, p4tree ).
Ces derniers sont mis ensemble sous forme de graphique comme on peut le
voir ci-dessus, avec différent types d’agglomération.

## Abundance value transformation

``` r
plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(physeq, Phylum %in% c("Firmicutes"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "sex",y = "Abundance",
                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}
```

### Commentaire : Ici nous transformons les données de comptage pour tenir compte de différents paramètres (tels que la taille de la bibliothèque, la variance, l’echelle et autres). Par exemple, on peut voir ici que la fonction plot\_abundance() est utilisée, afin de définir un graphique d’abondance relative.

``` r
ps3ra = transform_sample_counts(ps3, function(x){x / sum(x)})
```

### Commentaire : En utilisant la fonction transform\_sample\_counts(), on a convertit les comptages de chaque échantillon en abondance relative (ou en proportion).

``` r
plotBefore = plot_abundance(ps3,"")
plotAfter = plot_abundance(ps3ra,"")
```

``` r
grid.arrange(nrow = 2,  plotBefore, plotAfter)
```

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->
\#\#\# Commentaire: Ici nous avons tracés les valeurs obtenues
précédemment avant et après la transformation. De ce fait, nous avons
une comparaison des abondances initiales avec les données transformées.
On a également combiné chaque parcelle en un seul graphique. Ainsi, la
figure obtenue montre bien la comparaison des abondances initiales et
des abondances relatives.

## Subset by taxonomy

``` r
psOrd = subset_taxa(ps3ra, Order == "Lactobacillales")
plot_abundance(psOrd, Facet = "Genus", Color = NULL)
```

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->
\#\#\# Commentaire : Dans le graphique ci-dessus on observe les
abondances relatives de l’ordre taxonomique des “Lactobacillales”, qui
sont regroupés par les genres mâle/femelle de l’hôte.

## Preprocessing

``` r
qplot(sample_data(ps)$age, geom = "histogram",binwidth=20) + xlab("age")
```

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->
\#\#\# Commentaire : on a pu effectuer des analyses de données
complémentaires. On peut voir que l’histogramme ci-dessus démontre que
l’âge des souris se répartit en plusieurs groupes.

``` r
qplot(log10(rowSums(otu_table(ps))),binwidth=0.2) +
  xlab("Logged counts-per-sample")
```

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->
\#\#\# Commentaire : La figure suivante ci-dessus, présente un
histogramme comparant les profondeurs de lectures brutes et transformées
par log.

``` r
sample_data(ps)$age_binned <- cut(sample_data(ps)$age,
                          breaks = c(0, 100, 200, 400))
levels(sample_data(ps)$age_binned) <- list(Young100="(0,100]", Mid100to200="(100,200]", Old200="(200,400]")
sample_data(ps)$family_relationship=gsub(" ","",sample_data(ps)$family_relationship)
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
out.wuf.log <- ordinate(pslog, method = "MDS", distance = "wunifrac")
```

    ## Warning in UniFrac(physeq, weighted = TRUE, ...): Randomly assigning root as --
    ## GCAAGCGTTATCCGGAATTACTGGGTGTAAAGGGTGAGTAGGCGGCATGGTAAGCCAGATGTGAAAGCCTTGGGCTTAACCCAAGGATTGCATTTGGAACTATCAAGCTAGAGTACAGGAGAGGAAAGCGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAAGAACACCAGTGGCGAAGGCGGCTTTCTGGACTGAAACTGACGCTGAGGCACGAAAGCGTGGGGAG
    ## -- in the phylogenetic tree in the data you provided.

``` r
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned") +
  labs(col = "Binned Age") +
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-39-1.png)<!-- -->
\#\#\# Commentaire : Par la suite, nous examinons l’analyse des
coordonnées principales (PCoA) avec la dissimilitude de Bray-Curtis sur
la distance Unifrac pondérée. On peut voir que l’ordination des données
révèle quelques valeurs abbérantes comme on peut le voir sur le
graphique.

``` r
rel_abund <- t(apply(otu_table(ps), 1, function(x) x / sum(x)))
qplot(rel_abund[, 12], geom = "histogram",binwidth=0.05) +
  xlab("Relative abundance")
```

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-40-1.png)<!-- -->
\#\#\# Commentaire : Dans le graphique ci-dessus, nous pouvons voir que
les échantillons abbérants ont été majoritaire pour une seule séquence
d’amplicon variable (ASV), comme on peut le voir à gauche.

## Different Ordination Projections

``` r
outliers <- c("F5D165", "F6D165", "M3D175", "M4D175", "M5D175", "M6D175")
ps <- prune_samples(!(sample_names(ps) %in% outliers), ps)
```

### Commentaire : Nous retirons en premier toutes les valeurs abbérantes

``` r
which(!rowSums(otu_table(ps)) > 1000)
```

    ## F5D145 M1D149   M1D9 M2D125  M2D19 M3D148 M3D149   M3D3   M3D5   M3D8 
    ##     69    185    200    204    218    243    244    252    256    260

### Commentaire : Nous retirons également les échantillons de moins de 1000 lectures.

``` r
ps <- prune_samples(rowSums(otu_table(ps)) > 1000, ps)
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
```

### Commentaire : Nous retirons également les échantillons de moins de 1000 lectures.

``` r
out.pcoa.log <- ordinate(pslog,  method = "MDS", distance = "bray")
evals <- out.pcoa.log$values[,1]
plot_ordination(pslog, out.pcoa.log, color = "age_binned",
                  shape = "family_relationship") +
  labs(col = "Binned Age", shape = "Litter")+
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->
\#\#\# Commentaire : Nous faisons une analyse des coordonnées
principales (PCoA) en utilisant l’indice de dissimilitude de
Bray-Curtis, entre les échantillons. Nous pouvons observer sur le
graphique ci-dessus,on voit qu’il y a un effet âge (=observation via les
couleurs: du plus jeune en orange, d’âge moyen en vert et du plus vieux
en bleu) chez les souris mâles/femelles et qui proviennent également de
différentes portées (observation sous forme de rond et de triangle).

``` r
out.dpcoa.log <- ordinate(pslog, method = "DPCoA")
evals <- out.dpcoa.log$eig
plot_ordination(pslog, out.dpcoa.log, color = "age_binned", label= "SampleID",
                  shape = "family_relationship") +
  labs(col = "Binned Age", shape = "Litter")+
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-45-1.png)<!-- -->
\#\#\# Commentaire : Par la suite nous effectuons une analyse de
coordonnées principales double (DPCoA) cette fois-ci, cela fournit une
représentation en bioplots des échantillons et des catégories
taxonomiques. Comme on peut le voir dans le graphique on observe des
résultats très condensés. Ici le premier axe explique 75% de la
variabilité.

``` r
plot_ordination(pslog, out.dpcoa.log, type = "species", color = "Phylum") +
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-46-1.png)<!-- -->
\#\#\# Commentaire : Nous avons pu identifier via la graphique
ci-dessus, les taxons responsabes des axes 1 et 2, comme on l’a pu voir
précédemment.

``` r
out.wuf.log <- ordinate(pslog, method = "PCoA", distance ="wunifrac")
```

    ## Warning in UniFrac(physeq, weighted = TRUE, ...): Randomly assigning root as --
    ## GCGAGCGTTGTCCGGAATTACTGGGTGTAAAGGGAGCGTAGGCGGGAGAGCAAGTCGACTGTGAAATCTATGGGCTTAACCCATAGCTGCGATCGAAACTGTTCATCTTGAGTGAAGTAGAGGCAGGCGGAATTCCTAGTGTAGCGGTGAAATGCGTAAATATTAGGAGGAACACCAGTGGCGAAGGCGGCCTGCTGGGCTTTAACTGACGCTGAGGCTCGAAAGCGTGGGTAGC
    ## -- in the phylogenetic tree in the data you provided.

``` r
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned",
                  shape = "family_relationship") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  labs(col = "Binned Age", shape = "Litter")
```

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-47-1.png)<!-- -->
\#\#\# Commentaire : Ici on a effectué une analyse des coordonnées
principales (PCoA) en utilisant Unifrac pondéré entre les échantillons.
Comme précédemment, nous voyons que le deuxième axe est associé à un
effet de l’âge, qui est assez similaire à celui du PCoA. Cela est tout à
fait normal, car dans les deux méthodes d’ordination pennent en compte
l’abondance. Mais on a pu voir que la méthode DPCoA a donné une
représentation et une interprétation des résultats beaucoup plus nette
pour la partie de l’axe 2.

## Why are the ordination plots so far from square?

### Aspect ratio of ordination plots

### PCA on ranks

``` r
abund <- otu_table(pslog)
abund_ranks <- t(apply(abund, 1, rank))
```

### Commentaire : Ici on a créer une nouvelle matrice qui représente les abondances par leur rang. Les microbes sont classés selon les rangs. Par exemple, dans un échantillon le microbe le plus petit sera classé dans le rang 1, le deuxième plus petit dans le rang 2, le troisième le plus petit dans le rang 3 etc…

``` r
abund_ranks <- abund_ranks - 329
abund_ranks[abund_ranks < 1] <- 1
```

### Commentaire : On a mis l’abondance des microbes (dont le rang est inférieur à un certain seuil) égal à 1. Et les autres sont décalés vers le bas. Cela a été effectué afin d’éviter des écartsimportant entre les rangs.

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following objects are masked from 'package:Biostrings':
    ## 
    ##     collapse, intersect, setdiff, setequal, union

    ## The following object is masked from 'package:XVector':
    ## 
    ##     slice

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(reshape2)
abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
```

    ## Joining, by = c("Var1", "Var2")

``` r
colnames(abund_df) <- c("sample", "seq", "abund", "rank")

abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
```

    ## Joining, by = c("Var1", "Var2")

``` r
colnames(abund_df) <- c("sample", "seq", "abund", "rank")

sample_ix <- sample(1:nrow(abund_df), 8)
ggplot(abund_df %>%
         filter(sample %in% abund_df$sample[sample_ix])) +
  geom_point(aes(x = abund, y = rank, col = sample),
             position = position_jitter(width = 0.2), size = 1.5) +
  labs(x = "Abundance", y = "Thresholded rank") +
  scale_color_brewer(palette = "Set2")
```

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-50-1.png)<!-- -->
\#\#\# Commentaire : Ici on a transformé le seuil de classement des
échantillons. Au niveau du graphique, on peut voir que pour certains
échantillons on a une association entre l’abondance et le rang.

``` r
library(ade4)
```

    ## 
    ## Attaching package: 'ade4'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     score

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     score

``` r
ranks_pca <- dudi.pca(abund_ranks, scannf = F, nf = 3)
row_scores <- data.frame(li = ranks_pca$li,
                         SampleID = rownames(abund_ranks))
col_scores <- data.frame(co = ranks_pca$co,
                         seq = colnames(abund_ranks))
tax <- tax_table(ps) %>%
  data.frame(stringsAsFactors = FALSE)
tax$seq <- rownames(tax)
main_orders <- c("Clostridiales", "Bacteroidales", "Lactobacillales",
                 "Coriobacteriales")
tax$Order[!(tax$Order %in% main_orders)] <- "Other"
tax$Order <- factor(tax$Order, levels = c(main_orders, "Other"))
tax$otu_id <- seq_len(ncol(otu_table(ps)))
row_scores <- row_scores %>%
  left_join(sample_data(pslog))
```

    ## Joining, by = "SampleID"

    ## Warning in class(x) <- c(setdiff(subclass, tibble_class), tibble_class): Setting
    ## class(x) to multiple strings ("tbl_df", "tbl", ...); result will no longer be an
    ## S4 object

``` r
col_scores <- col_scores %>%
  left_join(tax) 
```

    ## Joining, by = "seq"

``` r
evals_prop <- 100 * (ranks_pca$eig / sum(ranks_pca$eig))
ggplot() + 
  geom_point(data = row_scores, aes(x = li.Axis1, y = li.Axis2), shape = 2) +
  geom_point(data = col_scores, aes(x = 25 * co.Comp1, y = 25 * co.Comp2, col = Order),
             size = .3, alpha = 0.6) +
  scale_color_brewer(palette = "Set2") +
  facet_grid(~ age_binned) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
       y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  coord_fixed(sqrt(ranks_pca$eig[2] / ranks_pca$eig[1])) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-51-1.png)<!-- -->
\#\#\# Commentaire : On voit que le biplot ci-dessus montre des
résultats similaires après la transformation des échantillons effectuée
précédemment. Les auteurs ont même suggérés que “cela renforce leur
confiance dans l’analyse des données originales”.

## Canonical correspondence

``` r
ps_ccpna <- ordinate(pslog, "CCA", formula = pslog ~ age_binned + family_relationship)
```

### Commentaire : Ici on donne des éléments nécessaires avant d’effectuer une autre méthode d’ordination, qui est l’analyse des correspondances canonique (CCpnA). Nous effectuons cela pour comparer les résultats avec ceux obtenus précédemment avec PCoA et DPCoA.

``` r
library(ggrepel)
ps_scores <- vegan::scores(ps_ccpna)
sites <- data.frame(ps_scores$sites)
sites$SampleID <- rownames(sites)
sites <- sites %>%
  left_join(sample_data(ps))
```

    ## Joining, by = "SampleID"

    ## Warning in class(x) <- c(setdiff(subclass, tibble_class), tibble_class): Setting
    ## class(x) to multiple strings ("tbl_df", "tbl", ...); result will no longer be an
    ## S4 object

``` r
species <- data.frame(ps_scores$species)
species$otu_id <- seq_along(colnames(otu_table(ps)))
species <- species %>%
  left_join(tax)
```

    ## Joining, by = "otu_id"

``` r
evals_prop <- 100 * ps_ccpna$CCA$eig[1:2] / sum(ps_ccpna$CA$eig)
ggplot() +
  geom_point(data = sites, aes(x = CCA1, y = CCA2), shape = 2, alpha = 0.5) +
  geom_point(data = species, aes(x = CCA1, y = CCA2, col = Order), size = 0.5) +
  geom_text_repel(data = species %>% filter(CCA2 < -2),
                    aes(x = CCA1, y = CCA2, label = otu_id),
            size = 1.5, segment.size = 0.1) +
  facet_grid(. ~ family_relationship) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
        y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  scale_color_brewer(palette = "Set2") +
  coord_fixed(sqrt(ps_ccpna$CCA$eig[2] / ps_ccpna$CCA$eig[1])*0.45   ) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-53-1.png)<!-- -->
\#\#\# Commentaire : Ici sur le biplot on peut observer qu’on a un score
de souris et de bactéries générés via CCpna. Nous pouvons observer que
certains ordres ressort plus majoritairement par rapport à d’autres. Par
exemple, on peut voir que l’ordre Clostridiales est en majorité, mais
également regroupé avec certaines espèces bactériennes isolées.

## Supervised learning

``` r
library(lattice)
```

``` r
library(caret)
sample_data(pslog)$age2 <- cut(sample_data(pslog)$age, c(0, 100, 400))
dataMatrix <- data.frame(age = sample_data(pslog)$age2, otu_table(pslog))
# take 8 mice at random to be the training set, and the remaining 4 the test set
trainingMice <- sample(unique(sample_data(pslog)$host_subject_id), size = 8)
inTrain <- which(sample_data(pslog)$host_subject_id %in% trainingMice)
training <- dataMatrix[inTrain,]
testing <- dataMatrix[-inTrain,]
plsFit <- train(age ~ ., data = training,
                method = "pls", preProc = "center")
```

### Commentaire : Ici on a divisé les données avec différentes affectations.

``` r
plsClasses <- predict(plsFit, newdata = testing)
table(plsClasses, testing$age)
```

    ##            
    ## plsClasses  (0,100] (100,400]
    ##   (0,100]        63         1
    ##   (100,400]       1        44

### Commentaire : En utilisant la fonction “predict” on a prédit les différents étiquettes de la classe en comparaison avec les données réelles. Et comme on peut le voir sur le tableau croisé ci-dessus, on a une bonne méthode de prédiction de l’âge.

``` r
library(ggplot2)
```

``` r
library(randomForest)
```

    ## randomForest 4.6-14

    ## Type rfNews() to see new features/changes/bug fixes.

    ## 
    ## Attaching package: 'randomForest'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     margin

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     combine

``` r
rfFit <- train(age ~ ., data = training, method = "rf",
               preProc = "center", proximity = TRUE)
rfClasses <- predict(rfFit, newdata = testing)
table(rfClasses, testing$age)
```

    ##            
    ## rfClasses   (0,100] (100,400]
    ##   (0,100]        64         5
    ##   (100,400]       0        40

### Commentaire : En utilisant la méthode précédente, voici un exemple avec des données aléatoires de forêts.

``` r
library(permute)
```

``` r
library(vegan)
```

    ## This is vegan 2.5-6

    ## 
    ## Attaching package: 'vegan'

    ## The following object is masked from 'package:caret':
    ## 
    ##     tolerance

``` r
pls_biplot <- list("loadings" = loadings(plsFit$finalModel),
                   "scores" = scores(plsFit$finalModel))
class(pls_biplot$scores) <- "matrix"

pls_biplot$scores <- data.frame(sample_data(pslog)[inTrain, ],
                                pls_biplot$scores)

tax <- tax_table(ps)@.Data %>%
  data.frame(stringsAsFactors = FALSE)
main_orders <- c("Clostridiales", "Bacteroidales", "Lactobacillales",
                 "Coriobacteriales")
tax$Order[!(tax$Order %in% main_orders)] <- "Other"
tax$Order <- factor(tax$Order, levels = c(main_orders, "Other"))
class(pls_biplot$loadings) <- "matrix"
pls_biplot$loadings <- data.frame(tax, pls_biplot$loadings)
ggplot() +
  geom_point(data = pls_biplot$scores,
             aes(x = Comp.1, y = Comp.2), shape = 2) +
  geom_point(data = pls_biplot$loadings,
             aes(x = 25 * Comp.1, y = 25 * Comp.2, col = Order),
             size = 0.3, alpha = 0.6) +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Axis1", y = "Axis2", col = "Binned Age") +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  facet_grid( ~ age2) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-61-1.png)<!-- -->
\#\#\# Commentaire : Pour continuer l’exemple, on a fait un biplot pour
séparer les échantillons par une variable de réponse.La différence ici
par rapport au graphique d’ordination c’est qu’on a identifié un
sous-espace pour maximiser la descrimination entre les classes.

``` r
rf_prox <- cmdscale(1 - rfFit$finalModel$proximity) %>%
  data.frame(sample_data(pslog)[inTrain, ])

ggplot(rf_prox) +
  geom_point(aes(x = X1, y = X2, col = age_binned),
             size = 1, alpha = 0.7) +
  scale_color_manual(values = c("#A66EB8", "#238DB5", "#748B4F")) +
  guides(col = guide_legend(override.aes = list(size = 4))) +
  labs(col = "Binned Age", x = "Axis1", y = "Axis2")
```

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-62-1.png)<!-- -->
\#\#\# Commentaire : On continue toujours avec l’exemple. Afin de
produire le graphique ci-dessous, une distance est calculée entre les
échantillons en prenant en compte la fréquence de répartition des
échantillons dans la même partition d’échantillon d’arbre. En utilisant
PCoA, on a un aperçu du mécanisme de classification des forêts
aléatoires.On peut alors savoir quel échantillon est facile à classé.

``` r
as.vector(tax_table(ps)[which.max(importance(rfFit$finalModel)), c("Family", "Genus")])
```

    ## [1] "Lachnospiraceae" NA

``` r
impOtu <- as.vector(otu_table(pslog)[,which.max(importance(rfFit$finalModel))])
maxImpDF <- data.frame(sample_data(pslog), abund = impOtu)
ggplot(maxImpDF) +   geom_histogram(aes(x = abund)) +
  facet_grid(age2 ~ .) +
  labs(x = "Abundance of discriminative bacteria", y = "Number of samples")
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-64-1.png)<!-- -->
\#\#\# Commentaire : Ici on a pu identifier le microbe qui a le plus
d’influence dans la prédiction de la forêt aléatoire. Ce microbe fait
partie de la famille des Lachnospiraceae et du genre Roseburia. Dans le
graphique on a en abscisse l’abondance des bactéries et en ordonné le
nombre d’échantillon de 0 à 100 jours (ce qui correspond au cadre en
haut) et de 100 à 400 jours (ce qui correspond au cadre en bas). Par
exemple : On peut voir que ce microbe est présent dans les échantillons.
Il est plus important de 100 à 400 jours.

## Graph-based analyses

### Creating and plotting graphs

``` r
library("phyloseqGraphTest")
library("igraph")
```

    ## 
    ## Attaching package: 'igraph'

    ## The following object is masked from 'package:vegan':
    ## 
    ##     diversity

    ## The following object is masked from 'package:permute':
    ## 
    ##     permute

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     as_data_frame, groups, union

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     union

    ## The following object is masked from 'package:XVector':
    ## 
    ##     path

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     union

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     union

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     normalize, path, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

``` r
library("ggnetwork")
net <- make_network(ps, max.dist=0.35)
sampledata <- data.frame(sample_data(ps))
V(net)$id <- sampledata[names(V(net)), "host_subject_id"]
V(net)$litter <- sampledata[names(V(net)), "family_relationship"]
```

``` r
library(igraph)
```

``` r
net_graph <-ggnetwork(net)
```

``` r
ggplot(net_graph, aes(x = x, y = y, xend = xend, yend = yend), layout = "fruchtermanreingold") +
  geom_edges(color = "darkgray") +
  geom_nodes(aes(color = id, shape = litter),  size = 3 ) +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        legend.key.height = unit(0.5,"line")) +
  guides(col = guide_legend(override.aes = list(size = .5)))
```

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-68-1.png)<!-- -->
\#\#\# Commentaire : Ici on a créer un réseau avec un seuil de matrice
de dissimilarité Jaccard. Les différentes formes correspondent aux
litières (1 et 2) et les couleurs correspondent aux différentes souris.
On peut observer dans ce réseau qu’il y a un regroupement
d’échantillons.

## Graph-based two-sample tests

## Minimum Spanning Tree (MST)

``` r
gt <- graph_perm_test(ps, "family_relationship", grouping = "host_subject_id",
                      distance = "jaccard", type = "mst")
gt$pval
```

    ## [1] 0.002

``` r
plotNet1=plot_test_network(gt) + theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))
plotPerm1=plot_permutations(gt)
grid.arrange(ncol = 2,  plotNet1, plotPerm1)
```

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-70-1.png)<!-- -->
\#\#\# Commentaire : On veut savoir si les deux portées proviennent de
la même distribution. Pour se faire, on effectue un test en utilisant
MST et l’indice de dissimilarité Jaccard. En gardant la structure des
données intacte on doit permuter les étiquettes. On peut observer le
graphique et l’histogramme de pertumation ci-dessus.

## Nearest neighbors

``` r
gt <- graph_perm_test(ps, "family_relationship", grouping = "host_subject_id",
                      distance = "jaccard", type = "knn", knn = 1)
```

``` r
plotNet2=plot_test_network(gt) + theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))
plotPerm2=plot_permutations(gt)
grid.arrange(ncol = 2,  plotNet2, plotPerm2)
```

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-72-1.png)<!-- -->
\#\#\# Commentaire : Ici on a effectué un réseau de plus proche voisin,
k-nearest, et l’histogramme de permutation. Le réseau montre que si une
paire d’échantillon est relié entre eux, cela veut dire qu’il son des
voisins et donc qu’il serait probable qu’ils soient dans la même portée.

## Linear modeling

``` r
library("nlme")
```

    ## 
    ## Attaching package: 'nlme'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     collapse

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     collapse

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     collapse

``` r
library("reshape2")
library(dplyr)
library(Biostrings)
library(IRanges)
ps_alpha_div <- estimate_richness(ps, split = TRUE, measure = "Shannon")
ps_alpha_div$SampleID <- rownames(ps_alpha_div) %>%
  as.factor()
ps_samp <- sample_data(ps) %>%
  unclass() %>%
  data.frame() %>%
  left_join(ps_alpha_div, by = "SampleID") %>%
  melt(measure.vars = "Shannon",
       variable.name = "diversity_measure",
       value.name = "alpha_diversity")
```

``` r
diversity_means <- ps_samp %>%
  group_by(host_subject_id) %>%
  summarise(mean_div = mean(alpha_diversity)) %>%
  arrange(mean_div)
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

``` r
ps_samp$host_subject_id <- factor(ps_samp$host_subject_id)
```

``` r
alpha_div_model <- lme(fixed = alpha_diversity ~ age_binned, data = ps_samp,
                       random = ~ 1 | host_subject_id)
```

``` r
new_data <- expand.grid(host_subject_id = levels(ps_samp$host_subject_id),
                        age_binned = levels(ps_samp$age_binned))
new_data$pred <- predict(alpha_div_model, newdata = new_data)
X <- model.matrix(eval(eval(alpha_div_model$call$fixed)[-2]),
                  new_data[-ncol(new_data)])
pred_var_fixed <- diag(X %*% alpha_div_model$varFix %*% t(X))
new_data$pred_var <- pred_var_fixed + alpha_div_model$sigma ^ 2
```

``` r
ggplot(ps_samp %>% left_join(new_data)) +
  geom_errorbar(aes(x = age_binned, ymin = pred - 2 * sqrt(pred_var),
                    ymax = pred + 2 * sqrt(pred_var)),
                col = "#858585", size = .1) +
  geom_point(aes(x = age_binned, y = alpha_diversity,
                 col = family_relationship), size = 0.8) +
  facet_wrap(~host_subject_id) +
  scale_y_continuous(limits = c(2.4, 4.6), breaks = seq(0, 5, .5)) +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Binned Age", y = "Shannon Diversity", color = "Litter") +
  guides(col = guide_legend(override.aes = list(size = 4))) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)),
        axis.text.x = element_text(angle = -90, size = 6),
        axis.text.y = element_text(size = 6))
```

    ## Joining, by = c("host_subject_id", "age_binned")

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-77-1.png)<!-- -->
\#\#\# Commentaire : Nous avons effectué ici une modélisation linéaire.
On a pris comme exemple, un modèle à effet mixte, afin d’étudier la
relation entre la diversité de la communauté microbienne de souris et
les variables d’âges, mais également la portée.

### Commentaire : Pourquoi on a utilisé un modèle a effet mixte ? Pour pouvoir formaliser nos observations.

### Commentaire : En premier lieu, on a utilisé l’indice de shannon pour chaque échantillon. Ensuite, on a assemblé le tout avec l’annotation de l’échantillon.

### Commentaire : Ensuite, on a réorganiser la facette de la diversité de la plus basse à la plus haute (pour les échantillons). Par la suite, on ajusté les valeurs, avec les barres d’erreur.

### Commentaire : On peut voir sur les graphiques ci-dessus les différentes portées selon les couleurs, avec sur les graphiques en oordonnée on a la diversité de Shannon et en abscisse on a l’âge limite.

## Hierarchical multiple testing

``` r
library("reshape2")
library("DESeq2")
```

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     sampleNames

``` r
library(matrixStats)
```

### Commentaire : nous effectuons ici des tests multiples hiérarchiques. Cela va permettre d’identifier les microbes individuels dont l’abondance est liée aux variables d’échantillonnage qui nous intéresse. Par la suite, nous faisons des calculs de transformation stabilisatrice de variance sur les objets.

``` r
sample_data(ps)$age_binned <- cut(sample_data(ps)$age,
                          breaks = c(0, 100, 200, 400))
levels(sample_data(ps)$age_binned) <- list(Young100="(0,100]", Mid100to200="(100,200]", Old200="(200,400]")
sample_data(ps)$family_relationship = gsub(" ", "", sample_data(ps)$family_relationship)
ps_dds <- phyloseq_to_deseq2(ps, design = ~ age_binned + family_relationship)
```

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

``` r
library(DESeq2)
```

``` r
geo_mean_protected <- function(x) {
  if (all(x == 0)) {
    return (0)
  }
  exp(mean(log(x[x != 0])))
}

geoMeans <- apply(counts(ps_dds), 1, geo_mean_protected)
ps_dds <- estimateSizeFactors(ps_dds, geoMeans = geoMeans)
ps_dds <- estimateDispersions(ps_dds)
```

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

``` r
abund <- getVarianceStabilizedData(ps_dds)
```

### Commentaire : ici on a fait une moyenne géométrique. Par exemple :if (all(x == 0)) {return (0)} –\> signifie qu’on fait une mise à zéro lorsque tous les coordonnées sont égales à zéro.

``` r
short_names <- substr(rownames(abund), 1, 5)%>%
  make.names(unique = TRUE)
rownames(abund) <- short_names
```

### Commentaire : On va faire les tests hiérarchiques. Mais avant on raccourcis les noms de chaque taxon/ASV.

``` r
abund_sums <- rbind(data.frame(sum = colSums(abund),
                               sample = colnames(abund),
                               type = "DESeq2"),
                    data.frame(sum = rowSums(otu_table(pslog)),
                               sample = rownames(otu_table(pslog)),
                               type = "log(1 + x)"))

ggplot(abund_sums) +
  geom_histogram(aes(x = sum), binwidth = 20) +
  facet_grid(type ~ .) +
  xlab("Total abundance within sample")
```

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-83-1.png)<!-- -->
\#\#\# Commentaire : Ici on a construit un histogramme qui donne
l’abondance totale dans chaque échantillon. On a deux histogrammes :
celui d’en haut c’est avec DESeq2 et l’autre en-dessous c’est avec le
log.

### Commentaire : En effet, on peut voir qu’il y a une variation de l’abondance totale entre chaque échantillon.

``` r
library("structSSI")
el <- phy_tree(pslog)$edge
el0 <- el
el0 <- el0[nrow(el):1, ]
el_names <- c(short_names, seq_len(phy_tree(pslog)$Nnode))
el[, 1] <- el_names[el0[, 1]]
el[, 2] <- el_names[as.numeric(el0[, 2])]
unadj_p <- treePValues(el, abund, sample_data(pslog)$age_binned)
```

### Commentaire : Pour effectuer des tests hiérarchique ce la nécessite des tests unidimensionnels pour chaque groupe taxonomique. Pour ce faire on a utilisé comme vous pouvez le voir, la fonction treePValues().

``` bash
wget https://cran.r-project.org/src/contrib/Archive/structSSI/structSSI_1.1.1.tar.gz
```

    ## --2020-12-02 15:09:39--  https://cran.r-project.org/src/contrib/Archive/structSSI/structSSI_1.1.1.tar.gz
    ## Resolving cran.r-project.org (cran.r-project.org)... 137.208.57.37
    ## Connecting to cran.r-project.org (cran.r-project.org)|137.208.57.37|:443... connected.
    ## HTTP request sent, awaiting response... 200 OK
    ## Length: 25591 (25K) [application/x-gzip]
    ## Saving to: ‘structSSI_1.1.1.tar.gz.1’
    ## 
    ##      0K .......... .......... ....                            100% 1.12M=0.02s
    ## 
    ## 2020-12-02 15:09:39 (1.12 MB/s) - ‘structSSI_1.1.1.tar.gz.1’ saved [25591/25591]

### Commentaire : pour pouvoir effectuer la suite des analyses on a récupérer des données venant d’un site internet, cran.r-project.org. Si on ne change pas de données les résultats ne seront pas similaires à ceux des auteurs.

``` r
library("structSSI")
```

``` r
library(devtools)
```

    ## Loading required package: usethis

    ## 
    ## Attaching package: 'devtools'

    ## The following object is masked from 'package:permute':
    ## 
    ##     check

``` r
install_local("./structSSI_1.1.1.tar.gz")
```

    ## Skipping 1 packages not available: multtest

    ##      checking for file ‘/tmp/Rtmpr2Ntsi/remotes6907271df8bb/structSSI/DESCRIPTION’ ...  ✓  checking for file ‘/tmp/Rtmpr2Ntsi/remotes6907271df8bb/structSSI/DESCRIPTION’ (363ms)
    ##   ─  preparing ‘structSSI’:
    ##    checking DESCRIPTION meta-information ...  ✓  checking DESCRIPTION meta-information
    ##   ─  checking for LF line-endings in source and make files and shell scripts
    ##   ─  checking for empty or unneeded directories
    ## ─  looking to see if a ‘data/datalist’ file should be added
    ##   ─  building ‘structSSI_1.1.1.tar.gz’
    ##      
    ## 

    ## Installing package into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)

``` r
library("structSSI")
el <- phy_tree(pslog)$edge
el0 <- el
el0 <- el0[nrow(el):1, ]
el_names <- c(short_names, seq_len(phy_tree(pslog)$Nnode))
el[, 1] <- el_names[el0[, 1]]
el[, 2] <- el_names[as.numeric(el0[, 2])]
unadj_p <- treePValues(el, abund, sample_data(pslog)$age_binned)
```

``` r
hfdr_res <- hFDR.adjust(unadj_p, el, .75)
summary(hfdr_res)
```

    ## Number of hypotheses: 764 
    ## Number of tree discoveries: 579 
    ## Estimated tree FDR: 1 
    ## Number of tip discoveries: 280 
    ## Estimated tips FDR: 1 
    ## 
    ##  hFDR adjusted p-values: 
    ##                 unadjp         adjp adj.significance
    ## GCAAG.95  1.861873e-82 3.723745e-82              ***
    ## GCAAG.70  1.131975e-75 2.263950e-75              ***
    ## GCAAG.187 5.148758e-59 1.029752e-58              ***
    ## GCAAG.251 3.519276e-50 7.038553e-50              ***
    ## GCAAG.148 1.274481e-49 2.548962e-49              ***
    ## GCAAG.30  9.925218e-49 1.985044e-48              ***
    ## GCGAG.76  1.722591e-46 3.445183e-46              ***
    ## GCAAG.167 6.249050e-43 1.249810e-42              ***
    ## 255       8.785479e-40 1.757096e-39              ***
    ## GCAAG.64  2.727610e-36 5.455219e-36              ***
    ## [only 10 most significant hypotheses shown] 
    ## --- 
    ## Signif. codes:  0 '***' 0.015 '**' 0.15 '*' 0.75 '.' 1.5 '-' 1

``` r
#interactive part: not run
plot(hfdr_res, height = 5000) # opens in a browser
```

### Commentaire : Ici on devrait avoir un lien qui nous emmène vers un sous- arbre contenant des bactéries avec des abondances différentielles. Cela a été obtenue à partir des tests hiérarchiques. Normalement on aurait dû voir des noeuds aux bactéries. Via l’arbre, on pouvait voir certaines associations entre le groupe d’âge et l’abondance bactérienne selon des groupes taxonomiques.

### Commentaire : Malheureusement, après l’ouverture de la page, celle-ci reste blanche, sans aucun arbre ni écriture. Ce qui est bien dommage.

``` r
tax <- tax_table(pslog)[, c("Family", "Genus")] %>%
  data.frame()
tax$seq <- short_names
```

### Commentaire : Les auteurs ont expliqués qu’il serait intéressant de donner un contexte aux résultats obtenus précédemment (malgrès que vous ne pouvez le voir dans mon fichier github). Pour ce faire on a récupéré l’identité taxonomique des hypothèses rejetées.

``` r
options(digits=3)
hfdr_res@p.vals$seq <- rownames(hfdr_res@p.vals)
tax %>%
  left_join(hfdr_res@p.vals) %>%
  arrange(adjp) %>% head(10)
```

    ## Joining, by = "seq"

    ##             Family            Genus       seq   unadjp     adjp
    ## 1  Lachnospiraceae             <NA>  GCAAG.95 1.86e-82 3.72e-82
    ## 2  Lachnospiraceae        Roseburia  GCAAG.70 1.13e-75 2.26e-75
    ## 3  Lachnospiraceae Clostridium_XlVa GCAAG.187 5.15e-59 1.03e-58
    ## 4  Lachnospiraceae             <NA> GCAAG.251 3.52e-50 7.04e-50
    ## 5  Lachnospiraceae Clostridium_XlVa GCAAG.148 1.27e-49 2.55e-49
    ## 6  Lachnospiraceae             <NA>  GCAAG.30 9.93e-49 1.99e-48
    ## 7  Ruminococcaceae     Ruminococcus  GCGAG.76 1.72e-46 3.45e-46
    ## 8  Lachnospiraceae Clostridium_XlVa GCAAG.167 6.25e-43 1.25e-42
    ## 9  Lachnospiraceae        Roseburia  GCAAG.64 2.73e-36 5.46e-36
    ## 10            <NA>             <NA>   GCAAG.1 5.22e-35 1.04e-34
    ##    adj.significance
    ## 1               ***
    ## 2               ***
    ## 3               ***
    ## 4               ***
    ## 5               ***
    ## 6               ***
    ## 7               ***
    ## 8               ***
    ## 9               ***
    ## 10              ***

### Commentaire : Après utilisation des données précédentes, on peut voir le résultats ci-dessus. En effet, on peut voir qu’il semblerait que les bactéries le plus associées, appartiennent toutes à la famille des Lachnospiraceae. On peut voir cela en observant la première colonne. Il faut se rappeler que ces résultats séraient en parfaite accord avec les résultats qu’on avait obtenus précédemment sur les fôrets. On avait également identifié la famille des Lachnospiraceae.

### Multitable techniques

``` r
metab <- read.csv("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/metabolites.csv",row.names = 1)
microbe_connect <-url("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/microbe.rda")
load(microbe_connect)
microbe
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 20609 taxa and 12 samples ]
    ## tax_table()   Taxonomy Table:    [ 20609 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 20609 tips and 20607 internal nodes ]

### Commentaire : Comme les données utilisées précédemment ne présentaient qu’un seule table, nous utilson de nouvelles données, qui ont été obtenus à partir des liens ci-dessus.Comme vous pouvez le voir on a deux types de données, l’un pour les métabolites (en premier) et l’autres pour les bactéries (en second). On peut voir qu’on a un total de 12 échantillons, avec 6 rangs taxonomiques et 20 607 noeuds internes.

``` r
library("genefilter")
```

    ## 
    ## Attaching package: 'genefilter'

    ## The following objects are masked from 'package:MatrixGenerics':
    ## 
    ##     rowSds, rowVars

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     rowSds, rowVars

``` r
keep_ix <- rowSums(metab == 0) <= 3
metab <- metab[keep_ix, ]
microbe <- prune_taxa(taxa_sums(microbe) > 4, microbe)
microbe <- filter_taxa(microbe, filterfun(kOverA(3, 2)), TRUE)
metab <- log(1 + metab, base = 10)
X <- otu_table(microbe)
X[X > 50] <- 50
dim(X)
```

    ## [1] 174  12

### Commentaire : Ici on a récupéré les données précédentes et on les a filtrés.

``` r
dim(metab)
```

    ## [1] 405  12

### Commentaire : En observant le résultats précédemment et celui-ci, on peut voir qu’on a 12 colonnes pour metab et X.Après avoir effectué cela on va comparer les caractéristiques dans des tableaux de données à haute dimension.Pour faire cela, on va appliquer une DPA clairsemée. Cette dernière correspond à une procédure de sélection à différence des méthodes d’ordination effectuées au cours des scripts précédemment.

``` r
library(PMA)
cca_res <- CCA(t(X),  t(metab), penaltyx = .15, penaltyz = .15)
```

    ## 123456789101112131415

### Commentaire : On a utilisé les paramètres “penaltyx”et “penaltyz” qui sont des pénalités de rateté. Pourquoi on utiliserait ces paramètres ? Par exemple, si on a des données de “penaltyx” faible cela va induire une réduction du nombre de microbes sélectionnés. Et nous avons cette même tendance pour la sélection des métabolites en utilisant “penaltyz”.

``` r
cca_res
```

    ## Call: CCA(x = t(X), z = t(metab), penaltyx = 0.15, penaltyz = 0.15)
    ## 
    ## 
    ## Num non-zeros u's:  5 
    ## Num non-zeros v's:  15 
    ## Type of x:  standard 
    ## Type of z:  standard 
    ## Penalty for x: L1 bound is  0.15 
    ## Penalty for z: L1 bound is  0.15 
    ## Cor(Xu,Zv):  0.974

### Commentaire : Comme on peut le voir, 5 microbes et 15 métabolites ont été sélectionnés. On peut voir qu’on a une corrélation de 0,974 entre les deux tableaux. On peut dire que les données bactériennes et métaboliques présentent des similaritudes en fonction des 20 caractéristiques sélectionnées.

``` r
combined <- cbind(t(X[cca_res$u != 0, ]),
                  t(metab[cca_res$v != 0, ]))
pca_res <- dudi.pca(combined, scannf = F, nf = 3)
```

### Commentaire : Les auteurs ont voulu relier les métabolistes et les OTU récupérés avec les caractéristiques des échantillons. Pour ce faire on a utilisé une PCA ordinaire.

``` r
genotype <- substr(rownames(pca_res$li), 1, 2)
sample_type <- substr(rownames(pca_res$l1), 3, 4)
feature_type <- grepl("\\.", colnames(combined))
feature_type <- ifelse(feature_type, "Metabolite", "OTU")
sample_info <- data.frame(pca_res$li, genotype, sample_type)
feature_info <- data.frame(pca_res$c1,
                           feature = substr(colnames(combined), 1, 6))
```

``` r
ggplot() +  geom_point(data = sample_info,
            aes(x = Axis1, y = Axis2, col = sample_type, shape = genotype), size = 3) + 
  geom_label_repel(data = feature_info,
                   aes(x = 5.5 * CS1, y = 5.5 * CS2, label = feature, fill = feature_type),
                   size = 2, segment.size = 0.3,
                   label.padding = unit(0.1, "lines"), label.size = 0) +
  geom_point(data = feature_info,
             aes(x = 5.5 * CS1, y = 5.5 * CS2, fill = feature_type),
             size = 1, shape = 23, col = "#383838") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("#a6d854", "#e78ac3")) +
  guides(fill = guide_legend(override.aes = list(shape = 32, size = 0))) +
  coord_fixed(sqrt(pca_res$eig[2] / pca_res$eig[2])) +
  labs(x = sprintf("Axis1 [%s%% Variance]",
                   100 * round(pca_res$eig[1] / sum(pca_res$eig), 2)),
       y = sprintf("Axis2 [%s%% Variance]",
                   100 * round(pca_res$eig[2] / sum(pca_res$eig), 2)),
       fill = "Feature Type", col = "Sample Type")
```

![](03_stat-analysisCC1_files/figure-gfm/unnamed-chunk-100-1.png)<!-- -->
\#\#\# Commentaire : En effet, on a fait une analyse en composantes
principales ordinaires (ou PCA ordinaire). On mis cela sous forme de
triplot avec des données multiples, comme on peut le voir. On a
différents types d’échantillons (PD et ST), différentes
caractéristiques (métabolite et OTU) et différents génotype (souris WT
et KO ou Knockout).

### Commentaire : En observant le triplot, on peut comparer les différents échantillons en fonctions des différents paramètres choisis. Par exemple, on peut voir que les échantillons sont séparés quand on se réfère aux types d’échantillons soit PD et ST, qui correspond au différents régimes alimentaires. Un autre exemple : Pour d’autres caractéristiques tel que le métabolite. En effet, on peut voir que les métabolites des souris WT sont assemblés en ligne.

## Conclusion :

### Au cours de l’application de ces scripts, on a pu avoir qu’à partir de données on pouvait effectuer différents types d’analyses et que cela dépendaient de ce qu’on cherche. On a pu voir qu’on pouvait projeter ces différents analyses de données via des arbres phylogénétiques, des biplots, des triplots, des analyses en réseaux et autres. Il est important de souligné qu’il est important de savoir ce que nous cherchons. Si nous ne le savons pas il sera difficile de décider par quoi nous allons commencer et surtout dans quel but.

### Via R, j’ai pu comprendre que c’est un outil statistique très important pour l’analyse de nos données.
