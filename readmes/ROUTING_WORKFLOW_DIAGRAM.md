# Routing Workflow Diagram

[← Back to Main README](../README.md)

---

This diagram shows how the wrapper handles the two entry modes (`--taxon` create vs.
`--input` query routing), and how each branches depending on assembly-set scale
(roughly "1k" vs "25k" genomes) and the `--megatree` / `--megatree-lazy` / `--placement`
flags. See `readmes/README_ADVANCED_FEATURES.md` §10 and `readmes/README_PIPELINE_PHASES.md`
for full prose detail on each step.

```mermaid
flowchart TD
    Start(["User provides input"]) --> Mode{"Entry mode?"}

    Mode -->|"--taxon NAME<br/>(create new DB)"| Download["Download raw genomes from NCBI<br/>(gather_filter_asms.sh --download-only)"]
    Mode -->|"--input assemblies.tsv<br/>(route query assemblies)"| Route["assembly_router.py:<br/>compare query taxonomy<br/>vs all DB taxonomies"]

    Route -->|"No match:<br/>NOVEL taxon"| Download
    Route -->|"Match found:<br/>DATABASE taxon"| TieBreak{"Tied DBs at top<br/>specificity? (megatree<br/>backbone + subclades<br/>share one taxonomy)"}

    %% ---------------- NOVEL / TAXON-CREATE PATH ----------------
    Download --> RawCheck{"Raw genome count vs<br/>--max-tree-genomes<br/>(default 2000)"}

    RawCheck -->|"under ceiling<br/>('1k' case)"| SingleTree["QC (CheckM2) then OrthoPhyl.sh<br/>ONE tree, ONE database<br/>(is_subclade=false)"]

    RawCheck -->|"over ceiling<br/>('25k' case)"| BigMode{"--megatree flag?"}

    BigMode -->|"No (default)"| Subsample["subsample_genomes.py:<br/>greedy MASH max-min down to<br/>--subsample-size (default 500)"]
    Subsample --> SingleTree

    BigMode -->|"Yes"| CeilingCheck{"Raw count vs<br/>--max-total-genomes<br/>(default 25000)"}
    CeilingCheck -->|"over ceiling"| Refuse["REFUSED:<br/>partition cost too high<br/>(~8hr / ~22GB at 75k)"]
    CeilingCheck -->|"under ceiling<br/>('25k' case)"| Partition["subclade_partition.py:<br/>MASH triangle + UPGMA<br/>into subclades of<br/>at most --subclade-size (150)"]

    Partition --> BuildEager{"--megatree-lazy?"}
    BuildEager -->|"No (default):<br/>build ALL subclades"| BuildAll["QC + OrthoPhyl per subclade<br/>(is_subclade=true)"]
    BuildEager -->|"Yes: build only subclades<br/>holding a query genome"| BuildSome["Build subclades WITH query;<br/>REGISTER (built=false)<br/>subclades with NO query"]

    BuildAll --> Backbone["Pick backbone-reps (5) diverse<br/>reps per subclade, build<br/>BACKBONE tree"]
    BuildSome --> Backbone
    Backbone --> Graft["megatree_graft.py:<br/>graft subclade trees onto<br/>backbone reps -> ONE merged tree"]
    Graft --> MegaDB["MULTIPLE databases created:<br/>backbone (is_backbone=true)<br/>+ one per subclade (is_subclade=true),<br/>all sharing one taxonomy string"]

    %% ---------------- DATABASE / EXISTING-DB PATH ----------------
    TieBreak -->|"No: single DB match<br/>('1k' DB - plain or<br/>subsampled single tree)"| ReLeafDirect["ReLeaf.sh directly<br/>onto matched database"]

    TieBreak -->|"Yes: megatree parent<br/>('25k' DB - backbone +<br/>subclades tied)"| Placement{"--placement flag?"}
    Placement -->|"backbone"| ToBackbone["Route to is_backbone=true DB<br/>(sparse overview, no MASH)"]
    Placement -->|"subclade (default)"| MashPick["MASH-sketch query, compare vs<br/>each tied subclade's sketch,<br/>pick nearest"]

    MashPick --> BuiltCheck{"Chosen subclade<br/>built=true?"}
    BuiltCheck -->|"Yes"| ReLeafDirect
    BuiltCheck -->|"No (lazy,<br/>never built)"| OnDemand["OrthoPhyl_subclade_build:<br/>QC + build that subclade from<br/>ITS OWN raw members,<br/>promote to built=true"]
    OnDemand --> ReLeafDirect
    ToBackbone --> ReLeafDirect
```

## Key defaults referenced

| Flag | Default | Meaning |
|---|---|---|
| `--max-tree-genomes` | 2000 | Ceiling for a single tree before capping strategy kicks in |
| `--subsample-size` | 500 | Diverse-subsample target size (default capping strategy) |
| `--max-total-genomes` | 25000 | Ceiling for the opt-in `--megatree` partition step (bounds the O(n²) MASH distance array) |
| `--subclade-size` | 150 | Max genomes per partitioned subclade |
| `--backbone-reps` | 5 | Diverse reps picked per subclade to build the megatree backbone |

## Reading the diagram

- **Novel-taxon branch is shared.** A `--taxon` create run and an unmatched query both
  fall through the same `Download → RawCheck` logic, so "1k vs 25k" plays out
  identically whether it was user-initiated or query-triggered.
- **Database-taxon branch only diverges once a query hits an existing DB.** If that DB
  was built small (~1k, no partitioning), it's a direct ReLeaf. If it was built as a
  megatree (~25k, partitioned), the tie-break/placement/MASH machinery runs first.
- **`--megatree-lazy` defers building, not partitioning.** The MASH triangle + UPGMA
  clustering step always runs up front for `--megatree`; lazy only decides whether an
  individual subclade's tree gets built immediately or registered as `built=false` for
  on-demand construction later.
