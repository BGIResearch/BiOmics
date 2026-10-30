# BiOmics-KG Data Source Selection Guide

This guide helps agents choose appropriate data sources for KG/RAG questions. The goal is not to query every database, but to first identify the required entity types, relation types, and evidence strength, then select 2–6 of the most relevant sources.

> BiOmics-KG coverage (from SupplementaryTable1): `11M` nodes, `356M` edges, and `16` node types.

## 0. Data model before querying with BRICK

BiOmics-KG stores database provenance on edges, not on nodes. Entity types and evidence sources must be handled separately.

- Nodes have no field for filtering by database source. Node types only describe entity categories, such as `Gene`, `Cell`, `Disease`, `Chemical`, `Pathway`, and `Protein`.
- Database and information sources are stored in edge metadata. Common output columns are `path.1.info_source` and `path.1.info_source_length`.
- `brick qr` does not currently send an `info_source` query parameter to `/v1/relation`. Filter sources after the query returns by inspecting the CSV or by running `brick reason normalize/score`.
- Multi-hop queries are optional. Use `brick qr --hops 2` or `brick qr --hops 3`, with a maximum of 3 hops. Inspect every intermediate node and every edge's source and relation in multi-hop results.
- `co_mentioned` is weak evidence only. Before using it, inspect `path.1.original_relation`, `path.1.info_source`, `path.1.condition`, and `path.1.relation_confidence`.

Example command templates:

```bash
brick qr -i genes.txt --st Gene --tt Cell -r marker_of --task celltype -n 100 -o marker_hits.csv
brick qr "TP53" --st Gene --tt Disease -r causally_related_to -n 50 -o tp53_disease.csv
brick qr "EGFR" --st Protein --tt Protein -r interact_with -n 100 -o egfr_ppi.csv
brick qr "TP53" --st Gene --tt Disease --hops 3 --directed -n 50 -o tp53_disease_paths.csv
```

Recommended decision order: identify the question and entity types, choose relation types, and then review evidence by `path.1.info_source` in the results. Do not query the entire graph by default.

## 1. Principles for agents

1. **Identify the task first**: entity normalization, path explanation, candidate discovery, drug mechanisms, disease phenotypes, and single-cell annotation require different sources.
2. **Resolve IDs in the primary database, then add evidence edges**: for example, map diseases to DOID before adding MONDO; map phenotypes to HPO; map proteins to UniProt before querying ClinVar/GOA/Reactome/DrugBank evidence.
3. **Ontologies provide semantics and hierarchies**: for example `is_a`, `has_part`, `equals_to`, canonical names, synonyms, and cross-database mappings.
4. **Public databases provide factual evidence**: such as gene-disease, drug-target, protein interaction, pathway participation, cell marker, and ligand-receptor evidence.
5. **Broad-coverage networks are best for candidate expansion**: PrimeKG, BioSNAP, STRING, and STITCH have broad coverage, but key conclusions should be verified in curated or source-specific databases.
6. **Relation names have different meanings across sources**: `interact_with` may mean different things in BioGRID, STITCH, and CellTalkDB. Always state the source context.
7. **Preserve raw spellings but interpret them correctly**: `develope_into` should be interpreted as `develops_into`; `ClinicalTrails` usually refers to a ClinicalTrials-like source, so queries should support both spellings.

## 2. Primary databases by node type

| Node type | Primary database | How agents should use it |
| --- | --- | --- |
| Species | NCBItaxon | Entry point for species, hosts, model organisms, and taxon lineages. |
| Tissue | UBERON | Entry point for tissues, organs, and anatomical structures. |
| Cell | CL | Entry point for cell types; use CellMarker/PanglaoDB for marker evidence. |
| Cell Component | GO | Entry point for cellular components. |
| Stage | HSAPDV, MMUSDV | Entry point for developmental stages: HSAPDV for humans and MMUSDV for mice. |
| Protein | UNIPROT | Entry point for protein entities; use GOA/Reactome_anno for functional annotation. |
| Gene | NCBI GENE | Entry point for gene entities; use Gene2Protein/UniProt for protein mappings. |
| Chemical | MESH | Entry point for chemicals; add ChEBI/PubChem for detailed chemical classes. |
| Disease | DOID | Entry point for diseases in KG tables; add MONDO for cross-database disease integration and mappings. |
| Phenotype | HPO | Entry point for human clinical phenotypes. |
| Pathway | KEGG | Entry point for pathways; use Reactome_anno for member annotation. |
| Function | GO | Entry point for molecular functions. |
| Process | GO | Entry point for biological processes. |

## 3. Dictionary/NER coverage summary

These counts come from SupplementaryTable4 and indicate dictionary coverage by entity type.

| Entity type | Dictionary name count |
| --- | --- |
| Cell | 5,845 |
| Disease | 32,140 |
| Phenotype | 61,082 |
| Stage | 931 |
| Pathway | 2,276 |
| Gene | 317,049 |
| Tissue | 68,383 |
| Function | 54,868 |
| Process | 156,990 |
| CellComponent | 14,157 |

## 4. Recommended sources by question type

| Question/task type | Preferred sources | Usage notes |
| --- | --- | --- |
| Entity normalization / synonyms / cross-database mappings | Bioportal, Mesh, DOID, MONDO, NCIT, SNOMED_CT, MEDGEN, CUI_MEDGEN_merged, UBERON, CL, GO, CHEBI, NCBItaxon, UNIPROT_mapping | Map to primary IDs first; do not infer mechanisms from equals_to/mapping alone. |
| Cell type annotation / scRNA markers | CL, CellMarker, CellMarker_mm, PanglaoDB, CellTaxonomy, CLO, Cellosaurus | Use CL for canonical cell types and CellMarker/PanglaoDB for marker evidence; use CLO/Cellosaurus for cell lines. |
| Cell communication / ligand-receptor | CellChat, CellTalkDB, CellTalkDB_mm, STITCH, IntAct, Reactome, SIGNOR | Use CellChat/CellTalkDB for ligand-receptor evidence and Reactome/SIGNOR for mechanisms. |
| Tissues, organs, anatomy, and spatial location | UBERON, BTO, EMAPA, MA, EFO, GO, NCIT | Use UBERON across species, EMAPA for mouse development, MA for adult mice, and BTO for sample provenance. |
| Developmental stage / time | HSAPDV, MMUSDV, EMAPA, EFO, UBERON | Use HSAPDV for humans, MMUSDV for mice, and EMAPA/UBERON for developing anatomical structures. |
| Basic gene information / gene IDs | NCBI_GENE, ENSEMBL, Gene2Protein, Gene2Protein_mm, UNIPROT_mapping | Use NCBI_GENE for genes; move to UniProt/GOA for proteins or functions. |
| Protein function / GO annotation | UNIPROT_human, UNIPROT_rodents, GO, GOA, GOA_mm, GO-PLUS | GO is the ontology; GOA provides gene/protein-to-GO annotation facts. |
| Protein interactions / complexes / network expansion | BioGRID_hs, BioGRID_mm, IntAct, StringDB, CORUM, SIGNOR, STITCH, BIOSNAP_GG-EnhancedTissue | BioGRID/IntAct are curated; STRING is broad; CORUM complexes are not equivalent to direct interactions. |
| Signaling / metabolic pathways / reactions | Reactome, Reactome_anno, KEGG, KEGG_anno, KEGG_anno_mm, RHEA, METACYC, SMPDB, PANTHER_PATHWAY, PW | KEGG is the main KG pathway entry; use Reactome for mechanisms and RHEA/MetaCyc for reactions; note KEGG licensing. |
| Disease normalization / disease hierarchy | MONDO, DOID, NCIT, MEDGEN, SNOMED_CT, ORPHANET, IDO | Use DOID in KG, MONDO for integration, NCIT for cancer, Orphanet for rare diseases, and IDO for infectious diseases. |
| Phenotypes / symptoms / disease-to-phenotype | HPO, HPO_anno, HPO_d2p_anno, HPO_g2p_anno, PATO, MPHENO, UPHENO, NBO | Use HPO for human phenotypes, MPHENO for animal models, and UPHENO for cross-species mappings. |
| Mutations / variants / genetic disease | ClinVar, ClinVar_gene, OMIM, OMIM_anno, GENO, HPO_g2p_anno, MONDO | Use ClinVar for clinical significance, OMIM for genetic mechanisms, DOID in KG, and MONDO for disease IDs. |
| Drug targets / mechanisms / drug-disease | DrugBank, DrugCentral, DRUGCENTRAL, DrugCentral_hs, DrugCentral_mm, IUPHAR, PHARMGKB, STITCH, BIOSNAP_ChG-Miner, PrimeKG | Prefer DrugBank/IUPHAR/DrugCentral for mechanisms, PharmGKB for PGx, and STITCH/PrimeKG for candidate networks. |
| Side effects / adverse reactions / combination therapy | SIDER, MedDRA, MEDRA, BIOSNAP_ChSe-Decagon, BIOSNAP_ChChSe-Decagon, BIOSNAP_ChCh-Miner, DrugBank | Use MedDRA for side-effect terminology, Decagon for combination effects, and DrugBank for drug interactions. |
| Drug repurposing / clinical evidence | REPODB, ClinicalTrails, DrugCentral, DrugBank, PrimeKG, BIOSNAP_DCh-Miner, BIOSNAP_DG-Miner | Use repoDB for successful and failed repositioning and ClinicalTrails/ClinicalTrials for trial evidence. |
| Transcriptional regulation / TF-target / PTM | TRRUST, SCENIC, SCENIC_mm, PhosphoSitePlus, SIGNOR, Reactome | TRRUST/PhosphoSitePlus are more curated; SCENIC provides inferred regulons for single-cell analysis. |
| Food / nutritional chemicals / metabolites | FOODB_FOOD, FOODB_COMPOUND, FOODON, HMDB, CHEBI, PUBCHEMCID, RHEA, METACYC | Use FOODON/FooDB for foods, ChEBI/PubChem for chemicals, and HMDB for human metabolites. |
| Species / hosts / model organisms | NCBItaxon, AGRKB_GENE, AGRKB_EXPRESSION, AGRKB_Disease, GOA_mm, Gene2Protein_mm, BioGRID_mm, UNIPROT_rodents, MA, EMAPA, MMUSDV | Use NCBItaxon for species and combine AGRKB with species-specific evidence. |

## 5. Relation type quick reference

| Edge type | Agent interpretation |
| --- | --- |
| is_a | Class/subclass or inheritance: A is a type of B. |
| has_part | Whole-part: A has B as a component. |
| equals_to | Cross-database equivalence or synonym mapping; not a mechanism. |
| species_relation | Relation to a species or taxonomic unit. |
| spatially_related_to | Spatial or anatomical relation, such as occurring in, located in, or adjacent to. |
| temporally_related_to | Relation to time or developmental stage. |
| sequentially_related_to | Sequence, ordering, mapping chains, or entity transitions; common in gene-to-protein and protein-family mappings. |
| processly_related_to | Relation to a biological process; retain the raw spelling `processly`. |
| functionally_related_to | Functional association; not necessarily physical interaction. |
| interact_with | Interaction, including protein-protein, chemical-protein, drug-target, and others. |
| participate_in | Participation in a process, pathway, reaction, or phenotype annotation. |
| produce | Production, expression, secretion, or generation by an entity; interpret from source context. |
| regulated | Regulatory relation, including activation, inhibition, transcriptional regulation, and PTM regulation. |
| pharmaco_effect_on | Pharmacological effect of a drug on a target, disease, phenotype, or process. |
| causally_related_to | Causal association, including variant pathogenicity, drug side effects, and gene-to-phenotype effects; inspect source evidence. |
| capable_of | Capability to perform a function or process. |
| transport_from_to | Transport direction or source-to-destination relation. |
| converse_from_to | Conversion or transformation direction, such as from A to B. |
| develope_into | Develops into; retain the raw spelling and interpret it as `develops_into`. |
| chemicals_relation / chemical_relation | Relations between chemicals, such as similarity, interaction, class membership, combination, or chemical association; meaning depends on the source. |
| disease_relation | Specialized relation between diseases or related to diseases. |
| compare | Comparison, such as similarity, difference, or relative properties. |
| no_relation | Retained no-relation, negative, or placeholder edge; do not use as positive biological evidence. |

## 6. Ontology source catalog

| Source | What it is | Suitable questions | KG node types | KG edge types | Usage notes |
| --- | --- | --- | --- | --- | --- |
| BTO | BRENDA Tissue Ontology, a vocabulary for tissues, cells, cell lines, and cultures from biological materials. | Normalizing tissue/cell provenance, sample provenance, tissue hierarchies, and tissue-produced entities. | Cell, Tissue | spatially_related_to, produce, causally_related_to, is_a, processly_related_to, has_part, develope_into | Combine with UBERON, CL, and CLO; BTO emphasizes experimental samples and enzyme-related context. |
| CHEBI | Chemical Entities of Biological Interest, an ontology for small molecules, metabolites, ions, and chemical classes. | Chemical normalization, chemical classes, chemical participation in reactions, and links from drugs/metabolites to processes. | Chemical | interact_with, produce, spatially_related_to, causally_related_to, is_a, functionally_related_to, chemicals_relation, equals_to, has_part, participate_in | Add PubChem, HMDB, and DrugBank for structure/ID queries and RHEA/Reactome for mechanisms. |
| CL | Cell Ontology, the canonical ontology for animal cell types across species. | Cell type identification, single-cell annotation, cell hierarchy, and developmental/functional relations. | Cell | spatially_related_to, produce, no_relation, causally_related_to, is_a, functionally_related_to, participate_in, processly_related_to, equals_to, temporally_related_to, transport_from_to, capable_of, compare, converse_from_to, sequentially_related_to, has_part, regulated, develope_into | Primary cell-type database; add CellMarker and PanglaoDB for marker evidence. |
| CLO | Cell Line Ontology, an ontology for cell lines and their provenance, species, diseases, and tissue background. | Cell-line normalization, provenance tracking, and relations among cell lines, diseases, species, and tissues. | Cell | spatially_related_to, produce, causally_related_to, is_a, species_relation, functionally_related_to, participate_in, processly_related_to, temporally_related_to, transport_from_to, capable_of, compare, converse_from_to, sequentially_related_to, has_part, regulated, develope_into | Combine with Cellosaurus: CLO for hierarchy and Cellosaurus for entity records and cross-references. |
| Cellosaurus | A knowledge base of cell-line names, synonyms, provenance, species, diseases, and external cross-references. | Cell-line name normalization, synonym resolution, identity checks, and provenance review. | Cell | Not declared in the table | Typically used for Cell entity dictionaries and mappings rather than declared edge types. |
| DOID | Disease Ontology, an open ontology of human diseases. | Disease hierarchies, disease normalization, and broad disease-class/subclass queries. | Disease | is_a, participate_in | DOID is the KG primary disease database; add MONDO/ORPHANET for integration or rare diseases. |
| EFO | Experimental Factor Ontology, an ontology for experimental factors, phenotypes, diseases, tissues, cell lines, species, and other metadata. | Normalizing multi-omics metadata, GWAS/experimental conditions, and sample attributes. | Other, Process, Protein, Species, Tissue | spatially_related_to, produce, is_a, functionally_related_to, transport_from_to, regulated, capable_of, has_part, participate_in, develope_into | Useful for linking dataset conditions to biological entities; combine with MONDO, UBERON, CL, and GO. |
| EMAPA | Mouse Developmental Anatomy Ontology, covering embryonic and developmental mouse anatomy. | Developing mouse tissue structures, embryonic expression sites, and developmental anatomical hierarchies. | Tissue | is_a, temporally_related_to, has_part, develope_into | Use MA for adult mouse anatomy and UBERON for cross-species tissues. |
| GO | Gene Ontology, a standardized description of molecular functions, biological processes, and cellular components. | Gene/protein function interpretation, GO enrichment, and function/process/component hierarchies. | Cell_Component, Function, Process | spatially_related_to, is_a, temporally_related_to, has_part, regulated | Core ontology for Function, Process, and Cell_Component; add GOA/GOA_mm for annotation facts. |
| GO-PLUS | GO Plus, an enhanced GO that links GO to external ontologies such as ChEBI, CL, and UBERON through logical axioms. | Cross-ontology reasoning, such as where a process occurs or which chemicals it involves. | Not declared in the table | spatially_related_to, produce, no_relation, causally_related_to, is_a, species_relation, functionally_related_to, processly_related_to, temporally_related_to, transport_from_to, regulated, capable_of, converse_from_to, sequentially_related_to, has_part, participate_in, develope_into | Best for inference/explanation; use GO for simple GO-term lookup. |
| HPO | Human Phenotype Ontology, a canonical vocabulary of human disease phenotypic abnormalities. | Symptom/sign normalization, disease-phenotype mapping, and gene-phenotype interpretation. | Phenotype | is_a | Primary phenotype database; add HPO_anno, HPO_d2p_anno, and HPO_g2p_anno for annotation facts. |
| HSAPDV | Human Developmental Stages, an ontology of human developmental stages. | Human age and developmental-stage normalization, including embryo, fetus, and adult stages. | Stage | is_a | One of the main Stage databases; use MMUSDV for mice. |
| IDO | Infectious Disease Ontology, covering infectious diseases, hosts, pathogens, transmission, and pathogenesis. | Infectious disease normalization, pathogen-related diseases, and host-pathogen semantics. | Disease | is_a, participate_in | Appears as both Ontology and Public Database; use IDO + MONDO for infectious diseases. |
| MA | Adult Mouse Anatomy, an ontology of adult mouse anatomical structures. | Mouse tissue/organ normalization, adult mouse expression sites, and mouse tissue hierarchy. | Tissue | is_a, has_part | Use EMAPA for embryonic/developmental anatomy and UBERON for cross-species integration. |
| MMUSDV | Mouse Developmental Stages, a vocabulary of mouse developmental stages. | Mouse embryonic days, postnatal stages, and developmental timing normalization. | Stage | is_a | One of the main Stage databases; use HSAPDV for humans. |
| MONDO | Mondo Disease Ontology, a unified disease ontology integrating OMIM, Orphanet, DOID, NCIt, and other vocabularies. | Disease normalization, cross-database mapping, and rare/genetic/cancer disease hierarchies. | Disease | causally_related_to, is_a, equals_to, disease_relation, capable_of | Excellent for disease harmonization; DOID is the KG primary disease database, with ClinVar, OMIM, DrugBank, and PrimeKG for mechanisms or drug relations. |
| MPATH | Mouse Pathology Ontology, covering mouse pathological structures, lesions, and pathological processes. | Animal-model pathology, histopathology descriptions, and model phenotype interpretation. | Disease | is_a | Combine with MPHENO, MA, and EMAPA. |
| MPHENO | Mammalian Phenotype Ontology, covering mammalian/mouse phenotypes. | Knockout phenotypes, animal-model phenotypes, and cross-species phenotype comparisons. | Phenotype | is_a | Use HPO for human phenotypes and UPHENO for cross-species integration. |
| NBO | Neuro Behavior Ontology, covering neurobehavior, behavioral experiments, and behavioral phenotypes. | Neurobehavioral phenotypes, behavioral processes, and behavioral experiment interpretation. | Process | spatially_related_to, causally_related_to, is_a, species_relation, functionally_related_to, transport_from_to, temporally_related_to, capable_of, converse_from_to, has_part, regulated | Best for neuroscience and behavior questions; combine with HPO/MPHENO. |
| NCIT | NCI Thesaurus, a comprehensive terminology for cancer, diseases, drugs, genes, mutations, and clinical/biomedical entities. | Cancer terminology, cancer diseases, drug/gene/mutation normalization, and clinical semantic queries. | Cell, Chemical, Disease, Gene, Mutation, Other, Pathway, Process, Species, Tissue | spatially_related_to, produce, causally_related_to, is_a, pharmaco_effect_on, equals_to, temporally_related_to, disease_relation, capable_of, sequentially_related_to, has_part, participate_in | Very useful for cancer; align disease IDs to DOID, add MONDO for integration, and SNOMED_CT/MedDRA for clinical terminology. |
| PATO | Phenotype And Trait Ontology, an ontology of qualities such as size, morphology, color, and abnormality level. | Decomposing phenotypes into entity + quality and handling traits such as enlargement, reduction, or abnormality. | Phenotype | spatially_related_to, causally_related_to, is_a, functionally_related_to, transport_from_to, equals_to, compare, has_part, regulated | Mainly a supplement to phenotype composition; do not use it alone to infer disease mechanisms. |
| UBERON | A cross-species anatomy ontology unifying tissues, organs, anatomical structures, and their hierarchies/part relations. | Tissue/organ normalization, cross-species tissue mapping, anatomical hierarchies, and spatial relations. | Tissue | spatially_related_to, produce, no_relation, causally_related_to, is_a, species_relation, functionally_related_to, processly_related_to, temporally_related_to, transport_from_to, regulated, capable_of, converse_from_to, sequentially_related_to, has_part, participate_in, develope_into | Primary tissue database; add MA/EMAPA for mouse-specific anatomy and BTO for experimental sample provenance. |
| UPHENO | Unified Phenotype Ontology, integrating phenotype semantics across species. | Comparing human, mouse, and other phenotypes and mapping model-organism phenotypes to humans. | Phenotype | interact_with, species_relation, is_a, has_part, regulated, develope_into, no_relation, chemicals_relation, capable_of, compare, sequentially_related_to, spatially_related_to, functionally_related_to, processly_related_to, temporally_related_to, participate_in, produce, causally_related_to, transport_from_to, converse_from_to | Preferred for cross-species phenotype reasoning; use HPO for human clinical phenotypes. |

## 7. Public database source catalog

| Source | What it is | Suitable questions | KG node types | KG edge types | Usage notes |
| --- | --- | --- | --- | --- | --- |
| AGRKB_Disease | Alliance of Genome Resources disease knowledge, integrating model-organism genes/alleles and disease associations. | Model-organism disease models, gene-disease evidence, and cross-species disease inference. | Not declared in the table | sequentially_related_to, produce, functionally_related_to, no_relation | Combine with OMIM, MONDO, ClinVar, and MPHENO for model-to-human disease bridging. |
| AGRKB_EXPRESSION | Alliance of Genome Resources expression knowledge, recording gene expression evidence across tissues, cells, and stages. | “In which tissues or stages is this gene expressed?” and model-organism expression sites. | Not declared in the table | produce | Combine expression facts with normalized entities in UBERON, MA, EMAPA, and CL. |
| AGRKB_GENE | Alliance of Genome Resources gene knowledge, integrating model-organism genes, orthology, and functional information. | Model-organism gene information, interaction/orthology clues, and gene function background. | Not declared in the table | interact_with | Use NCBI_GENE as the primary gene ID source and UniProt for protein information. |
| BIOSNAP_ChCh-Miner | BioSNAP drug-drug or chemical-chemical relation network. | Drug interactions, compound similarity/co-action networks, and initial drug-combination risk screening. | Not declared in the table | chemicals_relation | `ChCh` usually means chemical-chemical; add Decagon/SIDER for side effects. |
| BIOSNAP_ChChSe-Decagon | BioSNAP Decagon drug-combination and side-effect network. | Multi-drug adverse effects and drug-drug-side-effect prediction/explanation. | Not declared in the table | causally_related_to | Best for combination risks, not for single-drug indication evidence. |
| BIOSNAP_ChG-Miner | BioSNAP chemical-gene/drug-target network. | Drug targets, compound-gene interactions, and candidate target retrieval. | Not declared in the table | interact_with | Cross-check with DrugBank, DrugCentral, IUPHAR, and STITCH. |
| BIOSNAP_ChSe-Decagon | BioSNAP Decagon drug side-effect network. | Single-drug side effects and adverse-reaction signals. | Not declared in the table | causally_related_to | Add MedDRA/SIDER for clinical side-effect terminology. |
| BIOSNAP_DCh-Miner | BioSNAP disease-chemical/drug-disease network. | Disease-drug treatment/associations and repurposing candidates. | Not declared in the table | pharmaco_effect_on | A supplement to drug-disease edges; evaluate evidence alongside DrugCentral/ClinicalTrials/REPODB. |
| BIOSNAP_DG-Miner | BioSNAP disease-gene network. | Disease-associated genes, candidate pathogenic genes, and disease-module expansion. | Not declared in the table | pharmaco_effect_on | Combine with ClinVar, OMIM, HPO_g2p_anno, and PrimeKG. |
| BIOSNAP_GG-EnhancedTissue | BioSNAP tissue-enhanced gene-gene network, emphasizing tissue-context gene relations. | Tissue-specific gene networks and gene interactions/functional associations in tissue context. | Not declared in the table | interact_with | Combine with UBERON/MA/EMAPA for tissue normalization and STRING/BioGRID for PPI. |
| BioGRID_hs | BioGRID human data, with curated human protein, genetic, chemical, and selected PTM interactions. | Human PPI, gene interactions, chemical-protein interactions, and interaction evidence expansion. | Chemical, Protein | interact_with, pharmaco_effect_on | High-confidence curated interaction source; add STRING for broad functional associations. |
| BioGRID_mm | BioGRID mouse data, covering mouse protein, genetic, and chemical interactions. | Mouse PPI, model-organism gene interactions, and mouse gene networks. | Chemical, Protein | Not declared in the table | KG edge types are not declared, but it is useful as mouse interaction evidence. |
| Bioportal | NCBO BioPortal's ontology collection and cross-ontology mappings. | Cross-ontology equals_to/mapping, term lookup, and ontology ID alignment. | Not declared in the table | equals_to | Best for “what is this concept called in other databases?” rather than mechanistic evidence. |
| CORUM | A curated database of mammalian protein complexes. | Complex membership, complex function, and whether proteins belong to the same complex. | Protein | functionally_related_to | Combine with BioGRID/IntAct/STRING to distinguish complex co-membership from direct interaction. |
| CUI_MEDGEN_merged | A merged set of diseases, phenotypes, and other clinical concepts centered on MedGen/UMLS CUIs. | Clinical concept normalization, CUI alignment, and disease/phenotype synonym merging. | Disease, Other, Phenotype | Not declared in the table | Commonly used for entity merging and mapping; KG edge types are not declared. |
| CellChat | CellChatDB, a literature-supported ligand-receptor cell-communication resource. | Cell-cell communication inference, ligand-receptor signaling, and communication interpretation in single-cell/spatial transcriptomics. | Not declared in the table | regulated | Fits the `regulated` relation; combine with CellTalkDB, STITCH, and Reactome for signaling. |
| CellMarker | A curated database of human cell-type markers. | Single-cell type annotation, marker-to-cell queries, and marker interpretation within tissues. | Not declared in the table | produce | Prefer human markers; use CellMarker_mm or PanglaoDB for mouse markers. |
| CellMarker_mm | CellMarker mouse marker data. | Mouse single-cell annotation and marker queries in mouse tissues. | Not declared in the table | produce | Combine with PanglaoDB, CL, MA, and UBERON. |
| CellTalkDB | CellTalkDB, curated human/mouse ligand-receptor interaction pairs. | Cell communication, ligand-receptor pairs, and validation of intercellular signaling hypotheses. | Not declared in the table | interact_with | Complements CellChat; CellTalkDB edges usually look like `interact_with`. |
| CellTalkDB_mm | CellTalkDB mouse ligand-receptor interaction pairs. | Mouse cell communication and ligand-receptor inference. | Not declared in the table | interact_with | Combine with CellChat, Reactome, and STITCH. |
| CellTaxonomy | A cell-type classification and marker knowledge resource organized around single-cell identity. | Cell-type classification, marker-assisted annotation, and single-cell annotation candidates. | Not declared in the table | produce | Mainly `produce` edges; use CL to standardize cell names. |
| ClinVar | NCBI ClinVar, a public archive of the clinical significance of human variants and their disease/phenotype/drug-response associations. | Variant pathogenicity, variant-disease relations, and clinical-significance interpretation. | Mutation | causally_related_to | Primary mutation source; inspect submission details and review status for evidence strength. |
| ClinVar_gene | Gene-level disease/variant associations aggregated from ClinVar. | Finding diseases from genes or genes from diseases and obtaining a gene-level overview of variant evidence. | Not declared in the table | causally_related_to | Return to ClinVar for specific variants; add MONDO for disease normalization. |
| ClinicalTrails | Clinical-trial sources; the raw spelling is preserved and usually refers to ClinicalTrials.gov-like data. | Drug-disease clinical trials, intervention indications, and clinical evidence for repurposing. | Not declared in the table | pharmaco_effect_on | The raw spelling is `ClinicalTrails`; also try `ClinicalTrials` in searches. |
| DRUGCENTRAL | A KG import of DrugCentral covering drugs, targets, indications, mechanisms, and related information. | Drug actions, drug-target, drug-disease, and pharmacological-effect queries. | Not declared in the table | interact_with | DrugCentral/DrugCentral_hs/mm variants coexist; select by species and edge type. |
| DrugBank | DrugBank, a resource for drug structures, targets, mechanisms, interactions, and pharmacology. | Drug targets, mechanisms of action, drug-drug and drug-protein relations. | Chemical | pharmaco_effect_on, chemical_relation | Best for mechanistic interpretation; cross-check with DrugCentral, IUPHAR, and PharmGKB. |
| DrugCentral | DrugCentral, an open drug compendium containing active ingredients, products, indications, pharmacology, MoA, and bioactivity. | Drug indications, mechanisms of action, targets, and pharmacological classes. | Chemical | Not declared in the table | Node type is Chemical; some DrugCentral variants provide edges only. |
| DrugCentral_hs | DrugCentral human drug-target/drug-disease edges. | Human drug targets, human disease indications, and human pharmacological effects. | Not declared in the table | pharmaco_effect_on | For human drug questions, combine with DrugBank/IUPHAR/PHARMGKB. |
| DrugCentral_mm | DrugCentral mouse drug-target/drug-disease edges. | Drug action and pharmacological clues in mouse models. | Not declared in the table | pharmaco_effect_on | Be cautious about cross-species extrapolation; verify with human data. |
| ENSEMBL | Ensembl genome annotation resources for genes, transcripts, proteins, and genomic coordinates. | ID conversion, gene/transcript/protein localization, and genome annotation. | Protein | Not declared in the table | Node type is Protein; use NCBI_GENE as the primary gene database. |
| FOODB_COMPOUND | FooDB food-compound data, recording chemical components in foods and related pathways. | Food-compound relations, nutritional chemicals, food-derived metabolites, and pathway associations. | Chemical, Pathway | participate_in | Combine with FOODON, ChEBI, and HMDB. |
| FOODB_FOOD | FooDB food entity data, covering foods, ingredients, and source species. | Food classification, food-produced/contained compounds, and food-source species. | Other, Species | is_a, produce | Use FOODON for food normalization and ChEBI/PubChem for compound normalization. |
| FOODON | FoodOn, a food ontology covering products, sources, processing methods, and food classifications. | Food-name normalization, ingredient/food hierarchy, and food-to-biological-source relations. | Chemical | spatially_related_to, species_relation, is_a, functionally_related_to, temporally_related_to, has_part, participate_in, develope_into | Combine with FooDB for nutrition and food-chemistry questions. |
| NCBI_GENE | NCBI Gene, a cross-species gene information hub with GeneID, nomenclature, RefSeq, pathways, variants, and phenotype links. | Gene ID normalization, basic gene information, and cross-species gene queries. | Gene | Not declared in the table | Primary gene database; use UniProt for protein function and Gene2Protein for gene-protein mappings. |
| GENEPIO | Genomic Epidemiology Ontology, a vocabulary for genomic epidemiology, samples, experiments, sequencing, and public-health metadata. | Pathogen/genomic-epidemiology metadata, experimental workflows, and sample/sequencing semantics. | Other | spatially_related_to, produce, is_a, has_part, participate_in | Combine with IDO and NCBItaxon for infectious-disease/pathogen data. |
| GENO | Genotype Ontology, an ontology of genotypes, alleles, variants, and genetic composition. | Genotype/variant representation, allele semantics, and variant-to-phenotype structuring. | Other | interact_with, is_a, temporally_related_to, capable_of, sequentially_related_to, has_part, participate_in, develope_into | Use ClinVar for specific clinical variant evidence. |
| GOA | Gene Ontology Annotation, providing gene/protein-to-GO-term functional annotation. | Which functions a gene/protein has, which processes it participates in, and where it is located in cellular components. | Not declared in the table | spatially_related_to, no_relation, functionally_related_to, processly_related_to, capable_of, has_part, participate_in | Use GO for ontology hierarchies and GOA for annotation facts. |
| GOA_mm | Mouse Gene Ontology Annotation. | GO functions, processes, and cellular-component annotation for mouse genes/proteins. | Not declared in the table | spatially_related_to, no_relation, functionally_related_to, processly_related_to, capable_of, has_part, participate_in | Combine with Gene2Protein_mm, UNIPROT_rodents, MA, and EMAPA. |
| Gene2Protein | Gene-to-protein mapping tables. | Moving from gene symbols/GeneID to UniProt/protein nodes and connecting gene expression to protein function. | Not declared in the table | sequentially_related_to | Preferred when questions span gene-level and protein-level information. |
| Gene2Protein_mm | Mouse gene-to-protein mapping tables. | Mouse gene-to-protein mappings and mouse functional-annotation bridging. | Not declared in the table | sequentially_related_to | Combine with GOA_mm, UNIPROT_rodents, and BioGRID_mm. |
| HMDB | Human Metabolome Database, a resource for human metabolites, metabolic reactions, diseases, and biological samples. | Metabolite normalization, human metabolite classes, and metabolomics interpretation. | Chemical | is_a | Add ChEBI/PubChem for chemical normalization and RHEA/MetaCyc for reactions. |
| HPO_anno | HPO annotation data linking diseases, genes, and other entities to HPO phenotypes. | Which phenotypes a disease has, phenotype annotation facts, and phenotype enrichment/similarity. | Not declared in the table | participate_in | Use HPO for ontology structure and HPO_d2p_anno for disease-to-phenotype links. |
| HPO_d2p_anno | HPO disease-to-phenotype annotations. | Finding characteristic phenotypes for diseases or candidate diseases from symptom combinations. | Not declared in the table | Not declared in the table | Edge types are not declared, but semantically these are disease-phenotype annotations. |
| HPO_g2p_anno | HPO gene-to-phenotype annotations. | Finding phenotypes from genes, clinical manifestations of gene dysfunction, and candidate-gene prioritization. | Not declared in the table | causally_related_to | Combine with OMIM, ClinVar, and MONDO. |
| IDO | Infectious Disease Ontology, covering infectious diseases, hosts, pathogens, transmission, and pathogenesis. | Infectious disease normalization, pathogen-related diseases, and host-pathogen semantics. | Disease | is_a, participate_in | Appears as both Ontology and Public Database; use IDO + MONDO for infectious diseases. |
| IUPHAR | IUPHAR/BPS Guide to Pharmacology, an expert-curated resource for pharmacological targets and ligand interactions. | Receptors, ion channels, enzymes and other pharmacological targets, ligands, activities, and drug actions. | Chemical, Protein | pharmaco_effect_on | Preferred for target pharmacology; cross-check with DrugBank and DrugCentral. |
| IntAct | EMBL-EBI IntAct, a molecular interaction database curated from literature and submitted data. | Molecular interactions, PPI, interaction evidence, and protein-interaction networks. | Not declared in the table | interact_with | Similar to BioGRID in its curated focus; STRING is broader in association coverage. |
| KEGG | Kyoto Encyclopedia of Genes and Genomes, linking genomes to pathways, functions, diseases, drugs, and higher-level systems functions. | Pathway queries, gene/protein pathway participation, and functional modules. | Pathway, Protein | is_a | Note commercial/licensing restrictions; used in KG for Pathway, Protein, and is_a/annotation. |
| KEGG_anno | KEGG annotation linking genes/proteins to KEGG pathways. | Which KEGG pathways a gene/protein participates in. | Not declared in the table | participate_in | Combine with Reactome_anno and GOA for pathway/function interpretation. |
| KEGG_anno_mm | Mouse KEGG annotation. | Mouse gene/protein participation in KEGG pathways. | Not declared in the table | participate_in | Combine with Gene2Protein_mm and GOA_mm. |
| MEDGEN | NCBI MedGen, a clinical genetics and medical-concept resource aggregating diseases, phenotypes, syndromes, and mappings. | Clinical disease/phenotype terminology, CUI normalization, and genetic-disease concept mapping. | Not declared in the table | is_a | Combine with MONDO, HPO, OMIM, and ClinVar. |
| MEDRA | A KG source variant of MedDRA medical regulatory terminology. | Adverse events, clinical safety terminology, and drug side-effect normalization. | Disease | Not declared in the table | MedDRA also appears in the table; support case variants in searches. |
| METACYC | MetaCyc, a curated resource for metabolic pathways, enzymes, reactions, compounds, and genes. | Metabolic pathways, metabolic reactions, and enzyme-gene-compound relations. | Chemical, Gene, Pathway, Protein | produce | Combine with RHEA, Reactome, KEGG, and ChEBI. |
| MedDRA | Medical Dictionary for Regulatory Activities, a regulatory terminology commonly used for drug adverse events and clinical safety. | Side-effect/adverse-event normalization, clinical-term hierarchy, and MedDRA-to-other-vocabulary mappings. | Not declared in the table | is_a, equals_to | Combine with SIDER and BIOSNAP_ChSe-Decagon for drug side effects. |
| MedDRA2SNOMED_CT | MedDRA-to-SNOMED CT mappings. | Mapping adverse-event terms to clinical terminology for cross-database alignment. | Not declared in the table | equals_to | Use for equals_to mappings, not for causal claims. |
| Mesh | MeSH, the NLM medical subject headings, covering diseases, chemicals, proteins, processes, tissues, and many other concepts. | Literature-search term normalization, broad medical concept hierarchies, and pharmacological subject terms. | Cell, Chemical, Disease, Function, Other, Phenotype, Process, Protein, Species, Tissue | is_a, pharmaco_effect_on | Mesh is the primary Chemical database in the table; add ChEBI/PubChem for detailed chemical ontology. |
| NCBItaxon | NCBI Taxonomy, the classification and nomenclature system for organisms in public sequence databases. | Species normalization, taxon lineage, host/pathogen/model-organism identification. | Species | Not declared in the table | Primary Species database; useful for species relations in food, tissue, and disease questions. |
| OMIM | Online Mendelian Inheritance in Man, an authoritative catalog of human genes and genetic phenotypes/diseases. | Genetic diseases, gene-disease relations, Mendelian disorders, and phenotype descriptions. | Cell, Chemical, Disease, Other, Process, Protein, Species, Tissue | is_a, sequentially_related_to | Combine with MONDO, HPO, and ClinVar for genetic mechanisms and diagnostic candidates. |
| OMIM_anno | OMIM annotation data linking OMIM disease/gene records to other entities. | Genetic-disease annotation and supplementary disease-phenotype/gene relations. | Not declared in the table | participate_in | Use with HPO_d2p_anno and HPO_g2p_anno. |
| ORPHANET | The Orphanet rare-disease resource. | Rare-disease name normalization, classification, and rare-disease gene/phenotype clues. | Disease | Not declared in the table | Prefer combining with MONDO, HPO, and OMIM for rare diseases. |
| PANTHER_FAMILY | PANTHER protein-family data. | Protein families, homologous-function inference, and family hierarchies. | Protein | sequentially_related_to | Combine with UniProt, PANTHER_PATHWAY, and GOA. |
| PANTHER_PATHWAY | PANTHER pathway data. | Protein/gene pathway participation, especially function interpretation at family and pathway levels. | Pathway | participate_in | Reactome and KEGG are more common primary pathway entries. |
| PHARMGKB | PharmGKB, a pharmacogenomics knowledge base linking genes, variants, drugs, diseases, phenotypes, and pathways. | Drug-response genes, drug-variant associations, genotype effects on efficacy/toxicity, and PGx pathways. | Chemical, Gene, Mutation, Phenotype | causally_related_to, pharmaco_effect_on, transport_from_to, regulated, converse_from_to, participate_in | Preferred for drug-response and personalized-medicine questions; ClinVar focuses more on clinical variant significance. |
| PhosphoSitePlus | PhosphoSitePlus, a resource for post-translational modification sites and their regulatory/functional information. | Phosphorylation, acetylation, and other PTMs, kinase-substrate relations, and modification regulation. | Not declared in the table | regulated | Complements SIGNOR, Reactome, and BioGRID PTM evidence. |
| PSI_MI | PSI-MI molecular-interaction controlled vocabulary. | Interaction types, experimental methods, and standardization of molecular-interaction evidence terminology. | Other | is_a | Helps interpret evidence types in IntAct/BioGRID interaction edges. |
| PSI_MOD | PSI-MOD protein-modification ontology. | Standardization of protein-modification types such as phosphorylation and acetylation. | Other | has_part, is_a, functionally_related_to, develope_into | Combine with PhosphoSitePlus and UniProt PTM annotation. |
| PSI_MS | PSI-MS mass-spectrometry controlled vocabulary. | Mass-spectrometry experiments, instruments, and data-processing terminology. | Other | is_a, has_part | Best for proteomics/metabolomics metadata, not as a primary mechanism database. |
| PUBCHEMCID | PubChem Compound ID sources for compound ID/structure normalization. | Compound IDs, structures, synonym mappings, and PubChem CID alignment. | Chemical | Not declared in the table | Commonly used for Chemical entity mappings; KG edge types are not declared. |
| PW | Pathway Ontology, an ontology of pathway concepts and hierarchies. | Pathway-class normalization, pathway hierarchy, and pathway-term alignment. | Pathway | is_a, has_part | Use Reactome/KEGG/PANTHER/SMPDB for pathway facts and membership. |
| PanglaoDB | PanglaoDB, a human/mouse single-cell RNA-seq and cell-marker resource. | Cell-type markers, single-cell annotation, and mouse/human marker cross-validation. | Not declared in the table | produce | Combine with CellMarker, CL, and CellTaxonomy. |
| PrimeKG | Precision Medicine Knowledge Graph, integrating diseases, drugs, genes, proteins, pathways, and phenotypes from multiple sources. | Disease-mechanism overview, drug repurposing, and disease-drug-gene-phenotype path search. | Cell, Chemical, Disease, Function, Other, Phenotype, Process, Protein, Species, Tissue | is_a, causally_related_to, pharmaco_effect_on | Best as a broad-coverage network entry; verify key conclusions in original curated databases. |
| REPODB | repoDB, a drug-repositioning benchmark containing successful and failed drug-disease repositioning records. | Drug repurposing, approved/failed indications, and repositioning-method benchmarks. | Not declared in the table | pharmaco_effect_on | Combine with DrugCentral, ClinicalTrials, and DrugBank to assess evidence status. |
| Reactome | Reactome, a curated and peer-reviewed pathway and reaction database. | Pathway mechanisms, signaling/metabolic reactions, event ordering, and pathway interpretation. | Pathway | regulated, functionally_related_to, participate_in, temporally_related_to | High-quality curated pathway-mechanism source; KEGG is the KG primary pathway entry, while Reactome_anno provides member annotation. |
| Reactome_anno | Reactome annotation linking genes/proteins to Reactome pathway events. | Which Reactome pathways or reactions a gene/protein participates in. | Not declared in the table | participate_in | Combine with GOA and KEGG_anno for functional interpretation. |
| RHEA | Rhea, an expert-curated database of biochemical and transport reactions using ChEBI for participants. | Enzymatic reactions, transport reactions, metabolic reactions, and reaction-participant queries. | Not declared in the table | participate_in, is_a | Combine with ChEBI, UniProt, MetaCyc, and Reactome. |
| SCENIC | SCENIC regulon/TF-regulatory-network sources, commonly inferred from co-expression and motifs. | Single-cell TF activity, TF-target regulons, and regulatory-network interpretation. | Not declared in the table | regulated | An inferred network whose evidence strength differs from curated TF databases; add TRRUST/SIGNOR. |
| SCENIC_mm | Mouse SCENIC regulon/TF-regulatory-network sources. | Mouse single-cell TF activity and mouse TF-target regulation. | Not declared in the table | regulated | Combine with GOA_mm, Gene2Protein_mm, and TRRUST. |
| SIDER | SIDER, a drug side-effect/adverse-reaction resource. | Drug side effects, drug-adverse-event relations, and side-effect queries. | Not declared in the table | causally_related_to | Closely related to MedDRA terminology; combine with Decagon side-effect networks. |
| SIGNOR | SIGNOR, a curated signaling-causality database. | Activation, inhibition, regulation, and signaling-causal paths among proteins, chemicals, and phenotypes. | Chemical, Other, Phenotype, Protein | Not declared in the table | Best for mechanism-path interpretation; combine with Reactome, PhosphoSitePlus, and TRRUST. |
| SMPDB | Small Molecule Pathway Database, covering small-molecule metabolism, diseases, and drug pathways. | Metabolite/drug pathways, disease metabolic pathways, and small-molecule pathways. | Pathway | participate_in | Combine with HMDB, Reactome, KEGG, and RHEA. |
| SNOMED_CT | SNOMED CT, a clinical terminology system covering diseases, findings, procedures, anatomy, and substances. | Clinical terminology normalization, EHR concepts, and clinical expression of diseases/procedures/tissues. | Chemical, Disease, Other, Process, Protein, Species, Stage, Tissue | Not declared in the table | KG edge types are not declared; mainly used for clinical-entity mappings. |
| STITCH | STITCH, a chemical-protein interaction network integrating experiments, databases, text mining, and predictions. | Chemical-protein interactions, drug-target candidates, and chemical regulation of proteins. | Not declared in the table | interact_with, regulated, pharmaco_effect_on | Broad candidate coverage; verify mechanisms with DrugBank, IUPHAR, and literature. |
| StingDB | A KG StingDB source used for `interact_with`; the name may refer to a specific imported STRING/STING interaction table. | Protein/gene interaction candidates and network expansion. | Not declared in the table | interact_with | The name is not standard; cross-check with StringDB, BioGRID, and IntAct evidence. |
| StingDB_mm | A KG mouse StingDB source used for mouse `interact_with` relations. | Mouse protein/gene interaction-candidate expansion. | Not declared in the table | interact_with | Use the name cautiously; cross-check with BioGRID_mm, StringDB, and UNIPROT_rodents. |
| StringDB | STRING database, a protein-protein functional-association network integrating physical-interaction and functional-association evidence. | PPI/functional associations, network expansion, protein modules, and systems-biology analysis. | Protein | Not declared in the table | Broad coverage, but edges are not necessarily direct physical binding; use BioGRID/IntAct for direct interactions. |
| TRRUST | TRRUST, a curated transcription factor-target gene regulatory database. | TF-target relations, activation/inhibition, and transcriptional-regulation interpretation. | Not declared in the table | regulated | Distinguish from SCENIC: TRRUST is curated, while SCENIC is mainly inferred. |
| UNIPROT_human | UniProt human protein knowledge, including sequences, functions, domains, localization, PTM, and annotation. | Human protein IDs, functions, protein annotation, and protein-to-pathway/function queries. | Protein | sequentially_related_to, participate_in | Primary protein database; use Gene2Protein/UNIPROT_mapping for gene-to-protein mappings. |
| UNIPROT_mapping | UniProt ID mapping data. | Mapping among UniProt, GeneID, Ensembl, RefSeq, and other identifiers. | Not declared in the table | sequentially_related_to | Used for entity normalization and cross-database joins, not as a mechanism edge. |
| UNIPROT_rodents | UniProt rodent protein knowledge, covering mouse, rat, and other rodent proteins. | Mouse/rat protein functions, protein IDs, and species-specific annotation. | Protein | sequentially_related_to, participate_in | Combine with Gene2Protein_mm, GOA_mm, and BioGRID_mm. |

## 8. Query composition templates

### 8.1 Gene / protein function interpretation

**Recommended order**: `NCBI_GENE → Gene2Protein/UNIPROT_mapping → UNIPROT_human or UNIPROT_rodents → GOA/GOA_mm → KEGG_anno/Reactome_anno → BioGRID/IntAct/StringDB`.

Suitable for: Which protein does this gene encode? What does it do? Which pathways does it participate in? Which proteins does it interact with?

### 8.2 Disease / phenotype / variant interpretation

**Recommended order**: `DOID/MONDO → HPO → HPO_d2p_anno/HPO_g2p_anno → OMIM → ClinVar/ClinVar_gene → MEDGEN/NCIT/ORPHANET`.

Suitable for: Which phenotypes does a disease have? Is a gene/variant pathogenic? Which diseases match a combination of symptoms?

### 8.3 Drug / target / side effect / repurposing

**Recommended order**: `DrugBank/DrugCentral/IUPHAR → PHARMGKB → STITCH/BIOSNAP_ChG-Miner → SIDER/MedDRA/Decagon → REPODB/ClinicalTrails → PrimeKG`.

Suitable for: What is a drug's target? What is its mechanism? What are its side effects? Can it be repurposed for a disease?

### 8.4 Cell annotation / cell-cell communication

**Recommended order**: `CL → CellMarker/CellMarker_mm/PanglaoDB/CellTaxonomy → CellChat/CellTalkDB → Reactome/SIGNOR`.

Suitable for: Which cell types do these markers indicate? Which cells communicate through which ligand-receptor pairs?

### 8.5 Tissue / anatomy / development

**Recommended order**: `UBERON → BTO/EFO → MA/EMAPA/MMUSDV or HSAPDV → GO spatial/temporal relations`.

Suitable for: What is the parent structure of a tissue? Which tissues exist at a mouse developmental stage? Where does a process occur?

### 8.6 Metabolism / food / chemicals

**Recommended order**: `CHEBI/PUBCHEMCID → HMDB/FOODB_COMPOUND/FOODB_FOOD/FOODON → RHEA/METACYC/SMPDB/Reactome/KEGG`.

Suitable for: Which class does a compound belong to? Which foods contain it? Which metabolic reactions or pathways does it participate in?

## 9. Quality control and response guidelines

- For `equals_to`, state only “equivalent/mapping/synonymous”; do not interpret it as biological causality.
- For `causally_related_to`, always report the source; pathogenicity in ClinVar, side effects in SIDER, and network edges in BioSNAP have different evidence strengths.
- For `interact_with`, distinguish physical interactions, functional associations, predicted interactions, and ligand-receptor relations.
- For `produce`, interpret from source context: CellMarker may mean a cell produces/expresses a marker, while FooDB may mean a food/source produces a compound.
- Do not use `no_relation` as evidence for a positive answer.
- When sources conflict, the usual priority is: manually curated/official databases > specialized-domain databases > integrated KGs/predicted networks > mapping databases.
- For clinical, drug, and variant questions, state that the answer is not medical advice and should be checked against the original database and literature evidence.

## 10. External reference links

Use these links when agents need to verify database definitions, downloads, or official documentation.

| Resource | URL |
| --- | --- |
| GO | https://geneontology.org/docs/ontology-documentation/ |
| ChEBI | https://www.ebi.ac.uk/chebi/ |
| Cell Ontology | https://cell-ontology.github.io/ |
| Uberon | https://obofoundry.org/ontology/uberon.html |
| Reactome | https://reactome.org/ |
| KEGG | https://www.genome.jp/kegg/ |
| Rhea | https://www.rhea-db.org/ |
| DrugBank | https://go.drugbank.com/ |
| DrugCentral | https://drugcentral.org/ |
| PharmGKB | https://www.pharmgkb.org/ |
| IUPHAR Guide to Pharmacology | https://www.guidetopharmacology.org/ |
| ClinVar | https://www.ncbi.nlm.nih.gov/clinvar/ |
| HPO | https://hpo.jax.org/ |
| Mondo | https://mondo.monarchinitiative.org/ |
| OMIM | https://www.omim.org/ |
| CellMarker | https://xteam.xbio.top/CellMarker/ |
| PanglaoDB | https://panglaodb.se/ |
| CellChat | https://github.com/sqjin/CellChat |
| CellTalkDB | https://github.com/ZJUFanLab/CellTalkDB |
| BioGRID | https://thebiogrid.org/ |
| IntAct | https://www.ebi.ac.uk/intact/ |
| STRING | https://string-db.org/ |
| CORUM | https://mips.helmholtz-muenchen.de/corum/ |
| PrimeKG | https://github.com/mims-harvard/PrimeKG |
| BioSNAP | https://snap.stanford.edu/biodata/ |
| repoDB | https://unmtid-shinyapps.net/shiny/repodb/ |
| UniProt | https://www.uniprot.org/ |
| NCBI Gene | https://www.ncbi.nlm.nih.gov/gene/ |
| NCBI Taxonomy | https://www.ncbi.nlm.nih.gov/taxonomy |

---

Generated on 2026-07-01. The source table comes from `SupplementaryTable(1).xlsx` SupplementaryTable2; primary-database priorities come from SupplementaryTable3; NER coverage counts come from SupplementaryTable4.

## 11. Canonical source inventory and spelling notes (2026-07-07 audit)

This section is the machine-checkable KG source inventory aligned with `references/kg_source_policy.json`.
It preserves raw KG `info_source` spellings when the graph already uses them, while documenting canonical names used by BRICK routing.
Do not silently rewrite legacy raw tokens before querying; use canonical display names for reports and expand aliases for KG filters.

### 11.1 Complete policy inventory by evidence tier

| Tier | Source tokens |
| --- | --- |
| curated_clinical | `ClinVar`, `ClinVar_gene`, `OMIM`, `MEDGEN`, `ORPHANET`, `DrugBank`, `DrugCentral`, `DRUGCENTRAL`, `IUPHAR`, `PHARMGKB`, `ClinicalTrials`, `ClinicalTrails`, `REPODB` |
| curated_genetic_association | `GWAS_Catalog`, `GWAS_cataog` |
| curated_biological | `GOA`, `GOA_mm`, `GENEONTOLOGY`, `Reactome`, `Reactome_anno`, `KEGG`, `KEGG_anno`, `KEGG_anno_mm`, `RHEA`, `CORUM`, `IntAct`, `INTACT`, `BioGRID`, `BIOGRID`, `BioGRID_hs`, `BioGRID_mm`, `SIGNOR`, `UNIPROT`, `UNIPROT_human`, `UNIPROT_rodents` |
| single_cell_curated | `CellMarker`, `CellMarker_mm`, `PanglaoDB`, `CellTaxonomy`, `CellChat`, `CellTalkDB`, `CellTalkDB_mm`, `SCENIC`, `SCENIC_mm` |
| ontology | `CL`, `CLO`, `CELLLINEONTOLOGY`, `UBERON`, `GO`, `GO-PLUS`, `HPO`, `HUMAN PHENOTYPE ONTOLOGY`, `MONDO`, `DOID`, `BTO`, `CHEBI`, `EFO`, `PATO`, `NCIT`, `SNOMED_CT`, `MEDDRA`, `MESH`, `MPATH`, `NBO`, `UPHENO`, `UPHENO_ALL_WITH_RELATIONS`, `XREF`, `MA`, `EMAPA`, `FOODON` |
| integrated_network | `PrimeKG`, `Bioportal`, `AGRKB`, `AGRKB_Disease`, `AGRKB_EXPRESSION`, `AGRKB_GENE`, `BioSNAP_ChCh-Miner`, `BioSNAP_ChChSe-Decagon`, `BioSNAP_ChG-Miner`, `BioSNAP_ChSe-Decagon`, `BioSNAP_DCh-Miner`, `BioSNAP_DG-Miner` |
| predicted_or_text_mined | `StringDB`, `STITCH`, `StingDB`, `StingDB_mm`, `SIDER`, `PubMed text-mined` |

### 11.2 Canonical aliases and known raw spellings

| Canonical display name | Raw KG / accepted query aliases | Note |
| --- | --- | --- |
| `GWAS_Catalog` | `GWAS_Catalog`, `GWAS Catalog`, `GWAS_cataog`, `gwas_catalog`, `GWASCatalog`, `EBI GWAS Catalog` | `GWAS_cataog` is a legacy raw KG spelling; display and reports should use `GWAS_Catalog`. |
| `ClinicalTrials` | `ClinicalTrials`, `ClinicalTrails`, `ClinicalTrials.gov` | `ClinicalTrails` is a legacy raw KG spelling; query filters should expand to both names. |
| `StringDB` | `StringDB`, `StingDB`, `StingDB_mm` | `StingDB`/`StingDB_mm` are preserved raw KG source tokens and should be cross-checked with STRING/BioGRID/IntAct evidence. |
| `DrugCentral` | `DrugCentral`, `DRUGCENTRAL`, `DrugCentral_hs`, `DrugCentral_mm` | Uppercase and species-specific imports coexist in KG outputs. |
| `IntAct` | `IntAct`, `INTACT` | Case variants should be treated as the same source family for quality tiering. |
| `BioGRID` | `BioGRID`, `BIOGRID`, `BioGRID_hs`, `BioGRID_mm` | Generic and species-specific imports coexist in KG outputs. |
| `Gene Ontology` | `GO`, `GO-PLUS`, `GENEONTOLOGY`, `GOA`, `GOA_mm` | `GENEONTOLOGY` is a raw source token, while `GOA` is annotation evidence. |
| `HPO` | `HPO`, `HPO_anno`, `HPO_d2p_anno`, `HPO_g2p_anno`, `HUMAN PHENOTYPE ONTOLOGY` | Long-form raw source names should stay query-compatible. |
| `Cell Ontology` | `CL`, `CLO`, `CELLLINEONTOLOGY` | `CELLLINEONTOLOGY` is a raw all-caps import token. |
| `BioSNAP` | `BioSNAP_ChCh-Miner`, `BioSNAP_ChChSe-Decagon`, `BioSNAP_ChG-Miner`, `BioSNAP_ChSe-Decagon`, `BioSNAP_DCh-Miner`, `BioSNAP_DG-Miner`, `BIOSNAP_ChG-Miner`, `BIOSNAP_DCh-Miner`, `BIOSNAP_DG-Miner`, `BIOSNAP_GG-EnhancedTissue` | Mixed `BioSNAP`/`BIOSNAP` casing appears in docs and imports; preserve exact KG tokens when filtering. |
| `AGRKB` | `AGRKB`, `AGRKB_Disease`, `AGRKB_EXPRESSION`, `AGRKB_GENE`, `AGRKB:PMID:*`, `AGRKB:MGI:*` | Prefix-qualified values are AGRKB provenance IDs, not separate source families. |
| `PubMed text-mined` | `PMID:*`, `PubMed text-mined`, `mentioned`, `co_mentioned` | PMID-only provenance is weak literature trace unless snippet-supported by the evidence packet policy. |
| `GSEAPY` | `GSEAPY` | Observed in cached KG rows as enrichment/tool provenance; keep as unknown/diagnostic unless promoted by a future policy update. |

### 11.3 Spelling audit decisions

- Fixed for routing: `GWAS_cataog` canonicalizes to `GWAS_Catalog`, but query filters expand back to the raw legacy spelling.
- Fixed for routing: `ClinicalTrails` canonicalizes to `ClinicalTrials`, but query filters expand back to the raw legacy spelling.
- Preserved as raw KG tokens: `StingDB`, `StingDB_mm`, `GENEONTOLOGY`, `CELLLINEONTOLOGY`, `HUMAN PHENOTYPE ONTOLOGY`, `UPHENO_ALL_WITH_RELATIONS`, `XREF`, `DRUGCENTRAL`, `BIOGRID`, `INTACT`.
- Case-only variants such as `BioGRID`/`BIOGRID`, `IntAct`/`INTACT`, and `DrugCentral`/`DRUGCENTRAL` should share quality tiers.
- External REST sources such as GWAS Catalog API are not D-main KG-only evidence unless records are ingested into a frozen BiOmics-KG snapshot and inference queries only that snapshot.

---

Generated on 2026-07-01. Source inventory: SupplementaryTable(1).xlsx SupplementaryTable2; primary database priority: SupplementaryTable3; NER coverage counts: SupplementaryTable4.
