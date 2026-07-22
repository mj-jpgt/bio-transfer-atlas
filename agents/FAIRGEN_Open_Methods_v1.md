# FAIRGEN-Open: Mechanistic Decomposition of Polygenic Score Portability Risk
## Full Technical Methods — Architecture, Training, Evaluation, Stack

---

## 0. Conceptual Reframe

### Why the original framing fails

Standard PGS portability work measures *accuracy decay* as a function of ancestry distance and then attempts to correct it — either statistically (ancestry-residualized PRS) or representationally (adversarial disentanglement). Two results published in 2025–2026 undercut both approaches at the foundation:

1. **Harpak et al. (Nat Commun, 2026)**: Individual-level PRS accuracy is only weakly predicted by genome-wide genetic distance (R² < 0.5% for height), and is explained *comparably well* by socioeconomic measures. Portability trends are trait-specific in ways that reflect evolutionary history, not just ancestry distance. Precision and recall can move in opposite directions across the same ancestry gradient.

2. **Hu et al. (Nat Genet, 2025)**: LD differences and allele frequency divergence of causal variants — not heterogeneity in causal effect sizes — are the primary driver of low portability. Functional annotations help but cannot solve the problem, because the mechanism is structural.

These results collapse two independent things the field conflated: (a) ancestry-driven *distributional shift* of PRS values (measurable in 1000G, produced by AF/LD divergence) and (b) ancestry-driven *effect-size heterogeneity* (produced by differential selection, measurable in multi-ancestry GWAS comparison). Most fairness/portability tools target (a). Harpak et al. show that (b) matters differently by trait — and it's (b) that produces the catastrophic immune-trait failures, not (a).

### The core thesis

> *Portability risk is mechanistically heterogeneous: AF-divergence and LD-divergence produce recoverable distribution shift; selection-driven effect-size turnover produces unrecoverable structural failure. These mechanisms operate at different genomic scales (SNP → LD block → pathway) and can be jointly predicted from population-genetic features alone — without phenotypes — using cross-ancestry GWAS effect concordance as ground truth.*

This thesis makes three contributions simultaneously: it provides a **mechanistic explanation** for the trait-specific patterns Harpak et al. observe empirically; a **phenotype-free validation protocol** usable in 1000G without biobank access; and a **triage tool** that predicts which score-trait-pathway-population combinations will fail structurally vs. recover under standard corrections.

---

## 1. Data Sources

### 1.1 Individual-level genotypes

**1000 Genomes Project Phase 3** (GRCh38 alignment via IGSR)

- 2,504 individuals, 26 populations across 5 superpopulations (AFR, AMR, EAS, EUR, SAS)
- Includes admixed populations (ASW, ACB, ACB, CLM, MXL, PEL, PUR) — critical for secondary local-ancestry analysis
- Access: VCF per chromosome via `ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/`
- Phased haplotypes required for iHS and local ancestry inference

### 1.2 Population allele frequencies

**gnomAD v3.1.2** (GRCh38)

- 76,156 whole genomes across 8 genetic ancestry groups: AFR, AMR, ASJ, EAS, FIN, NFE, SAS, OTH
- Key fields: `AF_popmax`, per-ancestry `AF_*`, `nhomalt_*`, `AN_*` for missing-rate estimation
- Download: `gs://gcp-public-data--gnomad/release/3.1.2/vcf/genomes/` or via Hail on local cluster

### 1.3 Multi-ancestry GWAS summary statistics (concordance ground truth)

| Source | Ancestry | Access |
|---|---|---|
| Pan-UK Biobank | AFR, AMR, CSA, EAS, EUR, MID | `pan.ukbb.broadinstitute.org` |
| BioBank Japan | EAS (Japanese) | `jenger.riken.jp/pheweb` |
| FinnGen R10 | EUR (Finnish-enriched) | `r10.finngen.fi` |
| OpenGWAS | Meta-analyzed, variable | `gwas.mrcieu.ac.uk` |

Per-trait, per-ancestry effect estimates for the same variants enable computation of **cross-ancestry GWAS concordance labels** — the training target.

### 1.4 PGS Catalog scoring files

- Traits: T2D (∼40 published scores), CAD (∼30), BMI (∼25), LDL (∼20)
- Download via PGS Catalog REST API or bulk: `ftp.ebi.ac.uk/pub/databases/spot/pgs/`
- Retain scores with: genome build annotated, n\_variants > 50, development ancestry recorded
- These define the *index SNP sets* over which portability risk is computed

### 1.5 Pathway / gene-set resources

- Reactome hierarchy (DAG): `reactome.org/download/`
- MSigDB gene sets (C2: curated, C5: GO): `gsea-msigdb.org`
- Ensembl VEP offline cache for GRCh38
- gnomAD constraint scores (pLI, LOEUF) for gene-level evolutionary constraint

### 1.6 Selection statistics

- Pre-computed: 1000G PopHuman browser (`pophumanbrowser.upf.edu`) for iHS, PBS, Tajima's D per population
- Self-computed: iHS via selscan on 1000G phased haplotypes for populations not covered
- LINSIGHT and fitCons conservation scores for regulatory context

---

## 2. Data Processing Pipeline

### 2.1 Orchestration

All processing is managed by **Snakemake 8.x** with SLURM executor profile. DAG checkpointing used throughout. Docker/Apptainer containers pinned per rule.

```
workflow/
├── Snakefile
├── config.yaml
├── rules/
│   ├── 00_harmonize_genotypes.smk
│   ├── 01_compute_ld.smk
│   ├── 02_compute_afs.smk
│   ├── 03_selection_stats.smk
│   ├── 04_gwas_concordance.smk
│   ├── 05_feature_matrix.smk
│   ├── 06_train_model.smk
│   └── 07_evaluate.smk
└── envs/
    ├── genomics.yaml
    └── ml.yaml
```

### 2.2 Genotype preprocessing

```bash
# Lift PGS Catalog scoring files from GRCh37 to GRCh38 where necessary
CrossMap.py vcf hg19ToHg38.over.chain.gz {input.vcf} {ref.fa} {output.vcf}

# LD pruning for PCA (not for score computation — scores use all index SNPs)
plink2 --bfile 1kg_all \
  --indep-pairwise 1000kb 1 0.1 \
  --out pruned_for_pca

# Compute ancestry PCs (top 40)
plink2 --bfile 1kg_all \
  --extract pruned_for_pca.prune.in \
  --pca 40 approx \
  --out 1kg_pca
```

### 2.3 Population-specific LD reference panels

LD must be computed *within* each superpopulation reference panel, because LD-divergence across populations is one of the three mechanism components. Computed for each of the 5 superpopulations + all 26 individual populations.

```bash
# Per-superpopulation LD score computation (for LD-divergence features)
for POP in AFR AMR EAS EUR SAS; do
  plink2 --bfile 1kg_${POP} \
    --r2 --ld-window 1000kb --ld-window-r2 0.01 \
    --out ld_ref/${POP}
  
  # LD scores per SNP (sum of r² to all SNPs within 1cM)
  ldsc.py --bfile 1kg_${POP} \
    --l2 --ld-wind-cm 1 \
    --out ld_scores/${POP}
done
```

### 2.4 Variant-level feature extraction

For each index SNP `v` in each PGS Catalog scoring file:

**AF-divergence features (dim = 34)**
```python
import hail as hl

gnomad = hl.read_matrix_table("gs://gcp-public-data--gnomad/release/3.1.2/...")

# Per-population AFs: AFR, AMR, ASJ, EAS, FIN, NFE, SAS (7)
# Pairwise Fst relative to discovery ancestry (n_pop - 1 = 6 for EUR discovery)
# Max pairwise AF difference (1)
# Ancestry informativeness score In (1) = sum_p [p_i * log(p_i / p_bar)]
# Derived allele frequency (1)
# gnomAD popmax AF (1)
# gnomAD AF in admixed populations (AMR breakdown: 3)
# Missing rate per population (7)
# Total dim: 7 + 6 + 1 + 1 + 1 + 1 + 3 + 7 = 27 core + 7 flag = 34
```

**LD-divergence features (dim = 30)**
```python
# Per-superpopulation LD score (5)
# Per-superpopulation r² to nearest genome-wide-sig SNP in locus (5)
# Ratio of LD scores (EUR / target pop) — captures LD deflation (5)
# LD-block membership consistency across populations:
#   Bloque entropy = H(block assignments across populations) (1)
# r² to discovery GWAS index SNP, per superpopulation panel (5)
# Pairwise LD-score correlation between superpopulations at this locus (10)
# Total: 5 + 5 + 5 + 1 + 5 + 10 = 31 → rounded to 30 after dedup
```

**Selection-turnover features (dim = 45)**
```python
# iHS per population (26 populations from selscan): captures recent positive selection
# PBS (Population Branch Statistic) per triplet:
#   (AFR, EAS, EUR), (AFR, SAS, EUR), (EAS, SAS, EUR) → 3 values
# XP-CLR score vs EUR reference, per non-EUR superpop (4)
# Tajima's D per superpopulation (5)
# Fay & Wu's H per superpopulation (5)
# Per-SNP Fst (pairwise between discovery ancestry and each target pop) (5)
# gnomAD constraint: pLI, LOEUF, syn_z, mis_z for gene (4)
# phastCons100way, phyloP100way conservation (2)
# LINSIGHT score (1)
# Total: 26 + 3 + 4 + 5 + 5 + 5 + 4 + 2 + 1 = 55 → subset to 45 non-colinear
```

**Total per-SNP feature vector: 34 + 30 + 45 = 109 dimensions**
These are stored as an HDF5 feature store indexed by (rsid, GRCh38 position).

### 2.5 Cross-ancestry GWAS concordance labels

This is the **training target** — the phenotype-free ground truth for effect-size heterogeneity.

For each index SNP `v` in PGS Catalog with matching GWAS available in ≥2 ancestry-stratified sources:

```python
import pandas as pd
import numpy as np
from scipy import stats

def compute_concordance_labels(snp_id, trait, gwas_sources):
    """
    Returns per-SNP concordance metrics across all available ancestry-GWAS pairs.
    
    Labels:
        sign_concordance_matrix : (n_pop x n_pop) bool — same effect direction
        cochran_Q               : float — heterogeneity statistic
        I_squared               : float — proportion of variance due to heterogeneity
        beta_correlation        : float — Pearson r of betas across ancestries
        portability_risk_class  : {0: low, 1: medium, 2: high} — 3-class label
    """
    betas, ses, pops = [], [], []
    for source, pop in gwas_sources:
        row = source.query(f"SNP == '{snp_id}' and trait == '{trait}'")
        if len(row) > 0 and row['pval'].values[0] < 0.05:
            betas.append(row['beta'].values[0])
            ses.append(row['se'].values[0])
            pops.append(pop)
    
    if len(betas) < 2:
        return None  # insufficient data
    
    betas = np.array(betas)
    ses = np.array(ses)
    weights = 1 / ses**2
    
    # Cochran's Q
    beta_mean = np.average(betas, weights=weights)
    Q = np.sum(weights * (betas - beta_mean)**2)
    df = len(betas) - 1
    I2 = max(0, (Q - df) / Q)
    
    # Sign concordance: fraction of population pairs with same sign
    signs = np.sign(betas)
    n = len(signs)
    concordant_pairs = sum(
        signs[i] == signs[j]
        for i in range(n) for j in range(i+1, n)
    )
    total_pairs = n * (n - 1) / 2
    sign_concordance = concordant_pairs / total_pairs
    
    # Risk class: high risk = high I², low sign concordance
    risk_class = 0
    if I2 > 0.25 or sign_concordance < 0.80:
        risk_class = 1
    if I2 > 0.50 or sign_concordance < 0.60:
        risk_class = 2
    
    return {
        'I2': I2,
        'Q': Q,
        'sign_concordance': sign_concordance,
        'beta_corr': np.corrcoef(betas, np.arange(len(betas)))[0,1],
        'risk_class': risk_class,
        'n_ancestries': len(betas)
    }
```

**Expected labeled dataset size**: Across T2D, CAD, BMI, LDL with PGS Catalog scores and Pan-UKB multi-ancestry GWAS, we project ~15,000–45,000 labeled SNPs with ≥2 ancestry-stratified effect estimates after intersection. This is sufficient for supervised training with the architecture described below.

---

## 3. Architecture: Hierarchical Portability Risk Network (HPRN)

### 3.1 Overview

HPRN operates at three hierarchical levels: variant → LD-block/locus → gene → pathway. The key architectural novelties are:

1. **PopSpec Encoder**: treats the multi-population AF profile of each variant as a "population spectrogram" and encodes it via transformer self-attention over population tokens — learning population relationships from variant evolution rather than from fixed phylogenetic structure.

2. **LD-Graph Attention Network**: encodes the *divergence* of tagging structure across population-specific LD panels as graph edge features, capturing LD-divergence mechanism at the locus level.

3. **Structured Mechanism Bottleneck**: forces factored attribution of portability risk to the three mechanism components (AF-divergence, LD-divergence, selection-turnover) by routing each feature group through separate encoders with adversarial isolation, producing interpretable per-mechanism risk scores.

4. **Reactome-Hierarchical Pathway Pooling**: pathway-level risk aggregation following the Reactome DAG hierarchy, so that a portability-risk score is produced not just for terminal pathways but for all levels of the biological hierarchy.

### 3.2 PopSpec Encoder

```python
import torch
import torch.nn as nn
from einops import rearrange

class PopSpecEncoder(nn.Module):
    """
    Encode the AF profile of a variant across populations using
    transformer self-attention over population tokens.
    
    Each population is a 'token' with features:
        - AF of this variant in that population
        - iHS (if available, else 0 + missing flag)
        - Tajima's D
        - LD score in that population's LD panel
        - Learned population embedding (from population metadata)
    
    Self-attention over population tokens learns which populations
    cluster together in their treatment of this variant — effectively
    inferring local population structure from the variant's allele
    frequency landscape, without imposing a fixed phylogeny.
    """
    
    def __init__(
        self,
        n_populations=26,         # 1000G populations
        pop_feat_dim=6,           # per-population scalar features
        pop_embed_dim=16,         # learned population embedding dimension
        d_model=128,
        n_heads=8,
        n_layers=3,
        dropout=0.1
    ):
        super().__init__()
        
        # Learned population embeddings (encodes geographic/demographic context)
        self.pop_embedding = nn.Embedding(n_populations, pop_embed_dim)
        
        # Project per-population features to model dimension
        self.input_proj = nn.Linear(pop_feat_dim + pop_embed_dim, d_model)
        
        # Transformer encoder over population tokens
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=n_heads,
            dim_feedforward=d_model * 4,
            dropout=dropout,
            batch_first=True
        )
        self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=n_layers)
        
        # CLS-style aggregation token
        self.cls_token = nn.Parameter(torch.randn(1, 1, d_model))
        
        # Final projection
        self.out_proj = nn.Linear(d_model, d_model)
    
    def forward(self, pop_features, pop_ids, missing_mask=None):
        """
        pop_features : (batch, n_pop, pop_feat_dim)
        pop_ids      : (batch, n_pop) — population integer IDs
        missing_mask : (batch, n_pop) — True where AF is missing
        """
        B, N, _ = pop_features.shape
        
        pop_embeds = self.pop_embedding(pop_ids)  # (B, N, pop_embed_dim)
        x = torch.cat([pop_features, pop_embeds], dim=-1)  # (B, N, pop_feat_dim + pop_embed_dim)
        x = self.input_proj(x)  # (B, N, d_model)
        
        # Prepend CLS token
        cls = self.cls_token.expand(B, -1, -1)
        x = torch.cat([cls, x], dim=1)  # (B, N+1, d_model)
        
        # Extend mask for CLS token (always attend)
        if missing_mask is not None:
            cls_mask = torch.zeros(B, 1, dtype=torch.bool, device=x.device)
            key_padding_mask = torch.cat([cls_mask, missing_mask], dim=1)
        else:
            key_padding_mask = None
        
        x = self.transformer(x, src_key_padding_mask=key_padding_mask)
        
        # Return CLS representation as variant-level encoding
        return self.out_proj(x[:, 0, :])  # (B, d_model)
```

**Why this is novel**: No prior PGS portability tool encodes the multi-population AF profile as a structured sequence with learned inter-population attention. Existing tools use fixed Fst, fixed PCA distance, or concatenated AF vectors with no structural inductive bias. PopSpec learns *which populations covary* in their treatment of this variant — effectively learning local phylogeny from the data — which allows it to generalize to population pairs not seen during training.

### 3.3 LD-Graph Attention Network

```python
import torch_geometric.nn as pyg_nn
from torch_geometric.data import Data

class LDGraphEncoder(nn.Module):
    """
    Graph attention network over SNPs within an LD block.
    
    Nodes: index SNPs within ±500kb window
    Node features: PopSpec encoding + selection features
    
    Edges: pairwise LD r² between SNPs, per superpopulation panel
    Edge features: [r2_AFR, r2_AMR, r2_EAS, r2_EUR, r2_SAS]
                   — captures LD-divergence across populations
    
    Key: population-specific edge weights allow the GNN to learn
    that a tagging SNP may be in strong LD in EUR but weak LD in EAS,
    capturing the structural source of LD-driven portability failure.
    """
    
    def __init__(self, node_dim=128, edge_dim=5, hidden_dim=256, n_layers=3):
        super().__init__()
        
        # Edge feature MLP: population-specific LD r² → edge embedding
        self.edge_encoder = nn.Sequential(
            nn.Linear(edge_dim, 32),
            nn.GELU(),
            nn.Linear(32, 32)
        )
        
        # GATv2 layers with edge features
        self.convs = nn.ModuleList([
            pyg_nn.GATv2Conv(
                in_channels=node_dim if i == 0 else hidden_dim,
                out_channels=hidden_dim // 8,
                heads=8,
                edge_dim=32,
                concat=True,
                dropout=0.1
            )
            for i in range(n_layers)
        ])
        
        self.norms = nn.ModuleList([
            nn.LayerNorm(hidden_dim) for _ in range(n_layers)
        ])
        
        # Global mean pooling → locus embedding
        self.pool = pyg_nn.global_mean_pool
        
        # LD-divergence score head: explicit population-pair LD divergence
        self.ld_divergence_head = nn.Sequential(
            nn.Linear(hidden_dim, 64),
            nn.GELU(),
            nn.Linear(64, 10)  # pairwise LD divergence score for 5C2=10 superpop pairs
        )
    
    def forward(self, data):
        x, edge_index, edge_attr, batch = (
            data.x, data.edge_index, data.edge_attr, data.batch
        )
        
        edge_emb = self.edge_encoder(edge_attr)
        
        for conv, norm in zip(self.convs, self.norms):
            x = conv(x, edge_index, edge_attr=edge_emb)
            x = norm(x)
            x = torch.nn.functional.gelu(x)
        
        # Locus-level pooling
        locus_emb = self.pool(x, batch)  # (n_loci, hidden_dim)
        ld_divergence_scores = self.ld_divergence_head(locus_emb)
        
        return locus_emb, ld_divergence_scores
```

### 3.4 Structured Mechanism Bottleneck

This is the core architectural innovation for interpretability. Each mechanism component has its own encoder, and adversarial isolation ensures no feature leakage between components.

```python
class MechanismBottleneck(nn.Module):
    """
    Three factored encoders, each consuming only its mechanism's features.
    
    Adversarial isolation: a gradient-reversal discriminator attempts to
    predict which mechanism encoder a representation came from — pushing
    each encoder to be *complementary* to the others rather than redundant.
    
    Outputs:
        z_AF  : embedding of AF-divergence mechanism (dim=64)
        z_LD  : embedding of LD-divergence mechanism  (dim=64)
        z_SEL : embedding of selection-turnover mechanism (dim=64)
        z_combined : concatenation (dim=192) for prediction heads
    """
    
    def __init__(self, af_dim=34, ld_dim=30, sel_dim=45, bottleneck_dim=64):
        super().__init__()
        
        # Separate encoders per mechanism
        self.af_encoder = nn.Sequential(
            nn.Linear(af_dim, 128), nn.LayerNorm(128), nn.GELU(),
            nn.Linear(128, 128), nn.GELU(),
            nn.Linear(128, bottleneck_dim)
        )
        self.ld_encoder = nn.Sequential(
            nn.Linear(ld_dim, 128), nn.LayerNorm(128), nn.GELU(),
            nn.Linear(128, 128), nn.GELU(),
            nn.Linear(128, bottleneck_dim)
        )
        self.sel_encoder = nn.Sequential(
            nn.Linear(sel_dim, 128), nn.LayerNorm(128), nn.GELU(),
            nn.Linear(128, 128), nn.GELU(),
            nn.Linear(128, bottleneck_dim)
        )
        
        # Adversarial discriminator: which encoder did this come from?
        # Gradient reversal is applied before this during backprop
        self.mechanism_discriminator = nn.Sequential(
            nn.Linear(bottleneck_dim, 64), nn.GELU(),
            nn.Linear(64, 3)  # 3-class: AF / LD / SEL
        )
        
    def forward(self, af_features, ld_features, sel_features, grl_lambda=1.0):
        z_AF  = self.af_encoder(af_features)
        z_LD  = self.ld_encoder(ld_features)
        z_SEL = self.sel_encoder(sel_features)
        
        # Adversarial mechanism isolation
        # (gradient reversal applied in training loop)
        adv_logits = {
            'AF':  self.mechanism_discriminator(z_AF),
            'LD':  self.mechanism_discriminator(z_LD),
            'SEL': self.mechanism_discriminator(z_SEL)
        }
        
        z_combined = torch.cat([z_AF, z_LD, z_SEL], dim=-1)  # (B, 192)
        
        return z_AF, z_LD, z_SEL, z_combined, adv_logits
```

### 3.5 Reactome-Hierarchical Pathway Pooling

```python
from torch_geometric.nn import global_add_pool, global_mean_pool
import networkx as nx

class ReactomeHierarchicalPooling(nn.Module):
    """
    Aggregate variant/gene embeddings up the Reactome DAG hierarchy.
    
    Reactome is a DAG where pathways have parent-child relationships.
    We pool bottom-up: terminal pathways first, then parent pathways
    inherit weighted contributions from children.
    
    This produces a portability risk score at every level of biological
    organization — from specific molecular functions up to broad
    biological process categories.
    """
    
    def __init__(self, embedding_dim=192, n_pathways=2000):
        super().__init__()
        
        # Attention-weighted gene→pathway aggregation
        self.gene_to_pathway_attn = nn.MultiheadAttention(
            embed_dim=embedding_dim,
            num_heads=8,
            batch_first=True,
            dropout=0.1
        )
        
        # Pathway→parent aggregation MLP
        self.pathway_to_parent = nn.Sequential(
            nn.Linear(embedding_dim, embedding_dim),
            nn.LayerNorm(embedding_dim),
            nn.GELU()
        )
        
        # Per-pathway portability risk head
        self.pathway_risk_head = nn.Linear(embedding_dim, 3)  # 3-class risk
        
        # Per-mechanism attribution head at pathway level
        self.mechanism_attribution_head = nn.Sequential(
            nn.Linear(embedding_dim, 32),
            nn.GELU(),
            nn.Linear(32, 3),
            nn.Softmax(dim=-1)  # (p_AF, p_LD, p_SEL) attribution
        )
    
    def forward(self, gene_embeddings, reactome_dag, pathway_gene_membership):
        """
        Bottom-up pass through Reactome DAG.
        Returns pathway_risk_logits and mechanism_attribution per pathway.
        """
        pathway_embeddings = {}
        
        # Topological sort of DAG (leaves first)
        for pathway_id in nx.topological_sort(reactome_dag.reverse()):
            member_genes = pathway_gene_membership.get(pathway_id, [])
            
            if member_genes:
                gene_embs = torch.stack([
                    gene_embeddings[g] for g in member_genes
                    if g in gene_embeddings
                ])  # (n_genes, embedding_dim)
                
                # Attention-pool gene embeddings → pathway embedding
                query = gene_embs.mean(dim=0, keepdim=True).unsqueeze(0)
                key_value = gene_embs.unsqueeze(0)
                pathway_emb, _ = self.gene_to_pathway_attn(query, key_value, key_value)
                pathway_emb = pathway_emb.squeeze()
            else:
                # Internal node: aggregate from child pathways
                child_embs = torch.stack([
                    pathway_embeddings[child]
                    for child in reactome_dag.successors(pathway_id)
                    if child in pathway_embeddings
                ])
                pathway_emb = self.pathway_to_parent(child_embs.mean(dim=0))
            
            pathway_embeddings[pathway_id] = pathway_emb
        
        # Compute per-pathway outputs
        all_pathway_embs = torch.stack(list(pathway_embeddings.values()))
        risk_logits = self.pathway_risk_head(all_pathway_embs)
        attributions = self.mechanism_attribution_head(all_pathway_embs)
        
        return risk_logits, attributions, pathway_embeddings
```

### 3.6 Full HPRN Model

```python
class HPRN(nn.Module):
    """
    Hierarchical Portability Risk Network.
    
    Input:  per-SNP feature vectors (AF, LD, selection) + LD graph structure
    Output: per-pathway portability risk class + mechanism attribution
    
    Training targets:
        1. Cross-ancestry GWAS sign concordance (binary, per pop-pair)
        2. I² heterogeneity (continuous regression)
        3. 3-class portability risk (low / medium / high)
        4. Mechanism attribution (implicit, via bottleneck isolation)
    """
    
    def __init__(self, config):
        super().__init__()
        
        self.popspec = PopSpecEncoder(
            n_populations=config.n_populations,
            d_model=config.d_model
        )
        self.ld_graph = LDGraphEncoder(
            node_dim=config.d_model,
            edge_dim=5,
            hidden_dim=config.hidden_dim
        )
        self.mechanism_bottleneck = MechanismBottleneck(
            af_dim=34, ld_dim=30, sel_dim=45,
            bottleneck_dim=64
        )
        self.pathway_pooling = ReactomeHierarchicalPooling(
            embedding_dim=192  # 3 * 64
        )
        
        # Variant-level prediction heads
        self.sign_concordance_head = nn.Sequential(
            nn.Linear(192 + config.d_model + config.hidden_dim, 128),
            nn.GELU(),
            nn.Linear(128, 1),
            nn.Sigmoid()
        )
        self.heterogeneity_head = nn.Sequential(
            nn.Linear(192 + config.d_model + config.hidden_dim, 128),
            nn.GELU(),
            nn.Linear(128, 1)  # I² regression
        )
        self.risk_class_head = nn.Sequential(
            nn.Linear(192 + config.d_model + config.hidden_dim, 128),
            nn.GELU(),
            nn.Linear(128, 3)  # 3-class
        )
    
    def forward(self, batch):
        # 1. PopSpec encoding of AF landscape
        popspec_emb = self.popspec(
            batch['pop_features'],
            batch['pop_ids'],
            batch['missing_mask']
        )
        
        # 2. LD-graph locus encoding
        locus_emb, ld_divergence_scores = self.ld_graph(batch['locus_graph'])
        
        # 3. Mechanism factorization
        z_AF, z_LD, z_SEL, z_combined, adv_logits = self.mechanism_bottleneck(
            batch['af_features'],
            batch['ld_features'],
            batch['sel_features']
        )
        
        # 4. Fuse representations
        fused = torch.cat([z_combined, popspec_emb, locus_emb], dim=-1)
        
        # 5. Variant-level predictions
        sign_concordance = self.sign_concordance_head(fused)
        I2_pred = self.heterogeneity_head(fused)
        risk_logits = self.risk_class_head(fused)
        
        # 6. Pathway-level aggregation (done externally over batched loci)
        
        return {
            'sign_concordance': sign_concordance,
            'I2': I2_pred,
            'risk_logits': risk_logits,
            'z_AF': z_AF,
            'z_LD': z_LD,
            'z_SEL': z_SEL,
            'ld_divergence': ld_divergence_scores,
            'adv_logits': adv_logits
        }
```

---

## 4. Pre-Training: Self-Supervised Contrastive Learning on Population AF Profiles

Before supervised training on concordance labels (which cover only ~15–45k SNPs), the PopSpec encoder is pre-trained self-supervisedly over all ~84 million gnomAD variants — giving it a rich prior over population-genetic structure.

### 4.1 Contrastive objective

**Positive pairs**: two views of the same variant's AF profile, one with stochastic population masking (mask 3 random populations, impute from phylogenetic neighbors) + AF jitter (±5% multiplicative noise). These should have similar representations.

**Negative pairs**: variants from different genomic regions with different evolutionary histories (sampled to match AF distribution so the model can't cheat by memorizing AF level).

**Loss**: InfoNCE / NT-Xent over a queue of 4096 negatives (following MoCo v2 design, suitable for large genome-scale pre-training without memory issues):

```python
import torch.nn.functional as F

def nt_xent_loss(z1, z2, temperature=0.1, queue=None):
    """
    Normalized temperature-scaled cross-entropy.
    z1, z2: (batch, dim) L2-normalized embeddings of two views.
    queue:  (K, dim) L2-normalized negative queue (optional).
    """
    z1 = F.normalize(z1, dim=-1)
    z2 = F.normalize(z2, dim=-1)
    
    if queue is not None:
        negatives = queue.clone().detach()
        logits_pos = torch.sum(z1 * z2, dim=-1, keepdim=True) / temperature
        logits_neg = (z1 @ negatives.T) / temperature
        logits = torch.cat([logits_pos, logits_neg], dim=1)
        labels = torch.zeros(len(z1), dtype=torch.long, device=z1.device)
    else:
        B = z1.shape[0]
        z = torch.cat([z1, z2], dim=0)
        sim = (z @ z.T) / temperature
        sim.fill_diagonal_(-1e9)
        labels = torch.arange(B, device=z1.device)
        labels = torch.cat([labels + B, labels])
        logits = sim
    
    return F.cross_entropy(logits, labels)
```

**Pre-training data**: ~84M gnomAD v3.1.2 variants with AF in ≥3 populations. Compute variant-level features on-the-fly. Pre-train for 50 epochs, batch size 2048, on 4x A100 or equivalent, or 8x RTX 3090 for budget training. Expected wall time: ~72 hours.

Alternatively (no GPU cluster): use gnomAD variant subset restricted to chromosomes 1–3 (~25M variants) for a lighter pre-training run.

---

## 5. Supervised Training

### 5.1 Dataset construction

```python
from torch.utils.data import Dataset

class ConcordanceDataset(Dataset):
    """
    Labeled dataset of index SNPs with cross-ancestry GWAS concordance labels.
    
    Stratified by trait × score × discovery ancestry to avoid leakage:
        Train : T2D (EUR discovery) + CAD (EUR discovery)
        Val   : BMI (multi-ancestry) 
        Test  : LDL (EUR) + T2D (multi-ancestry discovery)
    """
    
    def __init__(self, snp_list, feature_store, label_store, reactome_dag,
                 ld_graph_store, pop_af_store):
        self.snps = snp_list
        self.features = feature_store   # HDF5
        self.labels = label_store       # pandas DataFrame
        self.dag = reactome_dag
        self.ld_graphs = ld_graph_store
        self.pop_afs = pop_af_store
    
    def __getitem__(self, idx):
        snp = self.snps[idx]
        return {
            'af_features':   self.features[snp]['af'],     # (34,)
            'ld_features':   self.features[snp]['ld'],     # (30,)
            'sel_features':  self.features[snp]['sel'],    # (45,)
            'pop_features':  self.pop_afs[snp],            # (26, 6)
            'pop_ids':       torch.arange(26),
            'locus_graph':   self.ld_graphs[snp],          # PyG Data object
            'sign_conc':     self.labels.loc[snp, 'sign_concordance'],
            'I2':            self.labels.loc[snp, 'I2'],
            'risk_class':    self.labels.loc[snp, 'risk_class']
        }
```

### 5.2 Loss function

```python
class HPRNLoss(nn.Module):
    
    def __init__(self, lambda_concordance=1.0, lambda_I2=0.5,
                 lambda_risk=1.0, lambda_adv=0.3, lambda_ld_div=0.2):
        super().__init__()
        self.lambdas = {
            'concordance': lambda_concordance,
            'I2':          lambda_I2,
            'risk':        lambda_risk,
            'adv':         lambda_adv,
            'ld_div':      lambda_ld_div
        }
    
    def forward(self, outputs, targets):
        # Primary: sign concordance prediction (BCE)
        L_concordance = F.binary_cross_entropy(
            outputs['sign_concordance'].squeeze(),
            targets['sign_conc'].float()
        )
        
        # Primary: I² heterogeneity regression (Huber loss for outlier robustness)
        L_I2 = F.huber_loss(outputs['I2'].squeeze(), targets['I2'].float())
        
        # Primary: 3-class portability risk (focal loss for class imbalance)
        L_risk = sigmoid_focal_loss(
            outputs['risk_logits'],
            targets['risk_class'],
            alpha=0.25, gamma=2.0
        )
        
        # Adversarial mechanism isolation:
        # Each mechanism encoder should NOT be able to predict which mechanism it is
        # (gradient reversal already applied; here we just compute cross-entropy)
        adv_labels = {
            'AF':  torch.zeros(len(outputs['adv_logits']['AF']), dtype=torch.long),
            'LD':  torch.ones(len(outputs['adv_logits']['LD']), dtype=torch.long),
            'SEL': torch.full((len(outputs['adv_logits']['SEL']),), 2, dtype=torch.long)
        }
        L_adv = sum(
            F.cross_entropy(outputs['adv_logits'][m], adv_labels[m].to(device))
            for m in ['AF', 'LD', 'SEL']
        )
        
        total = (
            self.lambdas['concordance'] * L_concordance +
            self.lambdas['I2']          * L_I2 +
            self.lambdas['risk']        * L_risk +
            self.lambdas['adv']         * L_adv
        )
        
        return total, {
            'L_concordance': L_concordance.item(),
            'L_I2': L_I2.item(),
            'L_risk': L_risk.item(),
            'L_adv': L_adv.item()
        }
```

### 5.3 Optimization

```python
optimizer = torch.optim.AdamW(
    model.parameters(),
    lr=1e-4,
    weight_decay=1e-2,
    betas=(0.9, 0.999)
)

# Warmup + cosine decay
scheduler = torch.optim.lr_scheduler.OneCycleLR(
    optimizer,
    max_lr=3e-4,
    total_steps=n_epochs * steps_per_epoch,
    pct_start=0.05,
    anneal_strategy='cos'
)

# Gradient clipping
torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
```

**Training regime**:
- Phase 1 (frozen PopSpec): train mechanism bottleneck + graph encoder + heads for 10 epochs with PopSpec frozen (loads pre-trained weights)
- Phase 2 (full fine-tuning): unfreeze all, train for 40 more epochs with 10x lower LR for PopSpec

**Batch size**: 512 variants per batch (each with their locus graph). Gradient accumulation ×4 for effective batch of 2048.

**Compute**: ~8 hours on 2x A100 80GB for full training. Single RTX 3090: ~36 hours with mixed precision.

---

## 6. Evaluation Protocol

### 6.1 Primary metrics

| Task | Metric | Baseline |
|---|---|---|
| Sign concordance (binary) | AUROC, AUPRC | Fst-only logistic regression |
| I² regression | Pearson r, Spearman ρ, RMSE | LD-score regression prediction |
| 3-class risk | Macro F1, calibration ECE | Hu et al. AF+LD feature rule |
| Pathway risk ranking | NDCG@10, Spearman ρ (pathway level) | Pathway Fst |
| Mechanism attribution | Silhouette score (ground truth: immune vs metabolic pathways) | n/a |

### 6.2 Baselines (mandatory)

Every contribution claim is tested against all of these:

1. **Global Fst only** — replicates Ding et al. approach
2. **AF-divergence only** — replicates Hu et al. primary finding
3. **LD-divergence only** — ablation of selection-turnover
4. **Selection statistics only** — ablation of AF/LD
5. **Concatenated features + random forest** — non-hierarchical flat baseline
6. **PopSpec encoder without pathway aggregation** — ablation of hierarchy
7. **HPRN without mechanism bottleneck** — ablation of attribution structure
8. **FairPRS** — existing fairness tool (reimplemented for concordance prediction)

### 6.3 Hold-out splits

**Trait hold-out**: Train on T2D + CAD; validate on BMI; test on LDL.
This tests generalization across biological domains.

**Population hold-out**: Train using Pan-UKB AFR/CSA/EAS as concordance sources; test concordance predictions against BioBank Japan (independent EAS GWAS) and FinnGen (independent EUR-Finnish GWAS). This tests whether the model generalizes to GWAS sources not in training.

**Score generation strategy hold-out**: Train on C+T scores; test on LDPred2 and PRS-CS scores. Tests whether portability risk predictions hold across PGS construction methodology.

### 6.4 Calibration

```python
from sklearn.calibration import calibration_curve
import matplotlib.pyplot as plt

def evaluate_calibration(risk_probs, labels, n_bins=10):
    """
    Reliability diagram + ECE for the 3-class portability risk classifier.
    Well-calibrated model: predicted probabilities match observed frequencies.
    """
    ece_per_class = []
    for cls in range(3):
        prob_true, prob_pred = calibration_curve(
            (labels == cls).astype(int),
            risk_probs[:, cls],
            n_bins=n_bins
        )
        ece = np.mean(np.abs(prob_true - prob_pred))
        ece_per_class.append(ece)
    return np.mean(ece_per_class)
```

### 6.5 Interpretability validation

**Biological ground truth test**: The mechanism attribution head should assign high *selection-turnover* weight to SNPs in immune/pathogen-response pathways and low weight to SNPs in core housekeeping pathways — consistent with Harpak et al.'s hypothesis. We test this by computing mean mechanism attribution per Reactome pathway category and measuring Spearman correlation with independent measures of evolutionary constraint (LOEUF, pLI) and selection signals (pathway-level PBS aggregated from population genetics databases).

**Replication of Harpak et al. finding**: For lymphocyte-count PGS index SNPs (the "smoking gun" trait from their paper), the model should assign high *selection-turnover* risk class and produce high I² predictions. This is a direct named prediction from the theory that the model should satisfy.

---

## 7. Secondary Analysis: Local Ancestry in Admixed 1000G Individuals

This is a direct empirical test of Harpak et al.'s stated open question: *"whether refined measures of genetic distance that capture local ancestry better explain portability."*

### 7.1 Local ancestry inference

```bash
# Phase haplotypes in admixed samples (already phased in 1000G Phase 3)
# Reference panels: unadmixed 1000G superpopulations

# Run RFMix2 on admixed populations
for POP in ASW ACB CLM MXL PEL PUR; do
  rfmix \
    --query-file 1kg_${POP}_phased.vcf.gz \
    --reference-file 1kg_reference_panels.vcf.gz \
    --sample-map 1kg_reference_superpops.map \
    --genetic-map 1kg_genetic_map_GRCh38.txt \
    --output-basename local_ancestry/${POP} \
    --n-generations 8 \
    --chromosome {chrom}
done
```

### 7.2 Local-ancestry-weighted PRS

For each admixed individual `i` and each index SNP `v`:

```python
def local_ancestry_weighted_prs(individual, pgs_scoring_file, local_ancestry_calls):
    """
    Compute local-ancestry-weighted PRS deviation.
    
    For each SNP, weight the dosage by the probability that the local
    ancestry matches the GWAS discovery population.
    
    Compare to: global-ancestry-weighted PRS (standard approach).
    
    Test: does knowing local ancestry at each SNP better predict
    per-SNP PRS deviation than knowing global ancestry PCs?
    This directly answers Harpak et al.'s open question.
    """
    global_prs = 0
    local_prs = 0
    
    for snp, weight in pgs_scoring_file.items():
        dosage = individual.genotype(snp)
        global_ancestry_correction = individual.global_pc_dist
        
        # Local ancestry: probability of EUR ancestry at this specific locus
        p_discovery_ancestry = local_ancestry_calls[individual][snp]['p_EUR']
        
        global_prs += dosage * weight
        local_prs  += dosage * weight * p_discovery_ancestry
    
    return global_prs, local_prs
```

**Test statistic**: Variance explained in PRS deviation by local ancestry at index SNP loci vs. global ancestry PCs, measured within each admixed population. If local ancestry explains significantly more variance, this validates HPRN's LD-divergence component (which is implicitly local) over global ancestry metrics.

---

## 8. Technology Stack

### 8.1 Core genomics

| Tool | Version | Purpose |
|---|---|---|
| PLINK 2.0 | ≥2.00a6 | LD computation, score calculation, Fst, PCA |
| PLINK 1.9 | ≥1.90b7 | LD clumping (legacy compatibility) |
| bcftools | ≥1.18 | VCF manipulation, normalization, intersection |
| tabix / bgzip | htslib ≥1.18 | VCF indexing |
| selscan | ≥2.0 | iHS, XP-CLR, XP-EHH per population |
| CrossMap | ≥0.6.6 | GRCh37 → GRCh38 liftover |
| RFMix2 | ≥2.03 | Local ancestry inference in admixed samples |
| Ensembl VEP | ≥112 | Variant → gene mapping (offline GRCh38 cache) |
| LDSC | ≥1.0.1 | LD score regression, LD score computation |

### 8.2 Python environment (ML / analysis)

```yaml
# environment.yaml
name: fairgen_open
channels: [conda-forge, bioconda, pytorch]
dependencies:
  - python=3.11
  - pytorch=2.3
  - torchvision=0.18
  - pytorch-cuda=12.1
  - torch-geometric=2.5
  - torch-scatter=2.1
  - torch-sparse=0.6
  - einops=0.7
  - hail=0.2.130           # gnomAD processing (optional; can use pandas)
  - polars=0.20            # fast DataFrame for large VCF tables
  - pandas=2.2
  - numpy=1.26
  - scipy=1.13
  - scikit-learn=1.5
  - statsmodels=0.14
  - networkx=3.3           # Reactome DAG manipulation
  - pyensembl=2.3          # Ensembl gene model
  - pysam=0.22             # VCF reading in Python
  - cyvcf2=0.30            # Fast VCF parsing (preferred over pysam for speed)
  - h5py=3.11              # HDF5 feature store
  - wandb=0.17             # Experiment tracking
  - hydra-core=1.3         # Config management
  - snakemake=8.16         # Pipeline orchestration
  - matplotlib=3.9
  - seaborn=0.13
  - plotly=5.22
  - rich=13.7              # CLI output formatting
```

### 8.3 R environment (statistical analysis)

```r
# renv.lock key packages
renv::install(c(
  "tidyverse",     # Data manipulation
  "data.table",    # Fast I/O for large summary stat files
  "ggplot2",       # Visualization
  "patchwork",     # Figure composition
  "broom",         # Tidy statistical outputs
  "metafor",       # Meta-analysis (Cochran's Q, I², heterogeneity)
  "GenomicRanges", # Genomic interval operations
  "rtracklayer",   # BigWig / BED import for conservation scores
  "corrplot",      # Correlation matrix visualization
  "ggrepel"        # Non-overlapping labels for Manhattan-style plots
))
```

### 8.4 Infrastructure

```
Compute:        Local GPU workstation (2x RTX 3090 24GB) or cloud (Lambda/Vast.ai)
Storage:        ~2TB for full 1000G + gnomAD features + LD graphs
Pipeline:       Snakemake 8.x with --executor slurm or --executor local
Containers:     Apptainer (Singularity) for genomics tools; Docker for ML
Tracking:       Weights & Biases (wandb) for all training runs
Versioning:     DVC for large data assets; git for code
Config:         Hydra for hierarchical experiment configuration
```

---

## 9. Expected Contributions

### 9.1 Primary theoretical contribution

**A structured mechanistic decomposition of PGS portability risk into three mutually-informative, empirically-separable components**: AF-divergence (recoverable via AF-informed weighting), LD-divergence (recoverable via local LD re-weighting), and selection-turnover (structurally unrecoverable — requires new GWAS in the target population). This is the first formal framework that *distinguishes* these three sources of portability failure from each other rather than treating "ancestry divergence" as a monolithic cause.

Concrete prediction: AF-divergence and LD-divergence portability risk should be *predictive of distributional PRS shift* in 1000G (measurable) but *not predictive of effect-size concordance discordance*. Selection-turnover risk should predict *both* — and should be disproportionately elevated in immune/pathogen-response pathways. These are falsifiable, quantitative predictions.

### 9.2 Primary architectural contribution

**PopSpec Encoder**: a population-spectrogram transformer that encodes the multi-population AF landscape of a variant via self-attention over population tokens. Unlike any existing model in the PGS literature, this learns the local population relationship structure from variant evolutionary history rather than imposing a fixed phylogeny or PC-distance prior. Pre-trained on 84M gnomAD variants, it functions as a general-purpose evolutionary context encoder for any downstream population-genetic task.

### 9.3 Primary methodological contribution

**Phenotype-free portability risk validation protocol**: cross-ancestry GWAS sign concordance and I² heterogeneity as the training and evaluation target, replacing phenotypic prediction accuracy (unavailable in 1000G and confounded by socioeconomic factors per Harpak et al.). This enables rigorous training and evaluation of portability risk models on open data without biobank access — and targets the *mechanism* that actually produces structural portability failure rather than its downstream accuracy symptom.

### 9.4 Secondary empirical contribution

**First direct test of local-ancestry-specificity in PRS portability**: using RFMix2-derived local ancestry in 1000G admixed populations to determine whether locus-specific ancestry is more predictive of per-SNP PRS deviation than global ancestry PCs. Direct answer to the first open question in Harpak et al. (2026), using the exact 1000G populations they gesture toward.

### 9.5 Resource contribution

**PRSM-PORT**: a pathway-resolved portability risk atlas for T2D, CAD, BMI, LDL across all 26 1000G populations, covering all PGS Catalog scores available for these traits. For each (score, trait, pathway, population) combination, PRSM-PORT reports:
- Predicted portability risk class (low/medium/high)
- Dominant mechanism (AF-divergence / LD-divergence / selection-turnover)
- Evidence from cross-ancestry GWAS concordance
- Actionable implication (can this pathway's score be corrected statistically, or does it require new EAS/AFR GWAS?)

This is the first trait-agnostic, pathway-resolved, mechanistically-explained portability risk resource built from entirely open data.

---

## 10. Claimed Advances Over SOTA

| Claim | Prior SOTA | This work |
|---|---|---|
| Portability mechanism decomposition | Hu et al. (AF/LD vs causal effects, binary) | Three-way (AF / LD / selection), quantified per SNP/pathway |
| Portability prediction target | R², AUC against phenotypic outcomes | Cross-ancestry effect concordance (I², sign) — phenotype-free |
| Population encoding | Fixed PCs or Fst scalar | PopSpec transformer over population tokens |
| Biological resolution | Score-level or SNP-level | Pathway-level, Reactome-hierarchical |
| Attribution | None (or post-hoc SHAP) | Structured mechanism bottleneck with adversarial isolation |
| Validation in admixed individuals | Not done empirically | RFMix2 local-ancestry experiment in 1000G admixed populations |
| Data requirement | UK Biobank / proprietary phenotype data | 100% open access |

---

*Document version: 1.0 | Framework: FAIRGEN-Open | Architecture: HPRN + PopSpec*
