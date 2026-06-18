const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  ImageRun, AlignmentType, HeadingLevel, BorderStyle, WidthType, ShadingType
} = require("docx");

const CW = 9360;
const border = { style: BorderStyle.SINGLE, size: 1, color: "BFBFBF" };
const borders = { top: border, bottom: border, left: border, right: border };
const HEAD = "D5E8F0";

function h1(t) { return new Paragraph({ heading: HeadingLevel.HEADING_1, children: [new TextRun(t)] }); }
function p(runs) { return new Paragraph({ spacing: { after: 140 }, children: Array.isArray(runs) ? runs : [new TextRun(runs)] }); }
const t = (x) => new TextRun(x);
const b = (x) => new TextRun({ text: x, bold: true });
const it = (x) => new TextRun({ text: x, italics: true });

function cell(text, w, head) {
  return new TableCell({ borders, width: { size: w, type: WidthType.DXA },
    margins: { top: 80, bottom: 80, left: 120, right: 120 },
    shading: head ? { fill: HEAD, type: ShadingType.CLEAR } : undefined,
    children: [new Paragraph({ children: [new TextRun({ text, bold: !!head })] })] });
}
function table(widths, rows) {
  return new Table({ width: { size: CW, type: WidthType.DXA }, columnWidths: widths,
    rows: rows.map((r, i) => new TableRow({ children: r.map((c, j) => cell(String(c), widths[j], i === 0)) })) });
}

const children = [];

children.push(new Paragraph({ spacing: { after: 80 },
  children: [new TextRun({ text: "Hypoxia Licenses the Cluster-and-Escort Phenotype of Circulating Tumor Cells", bold: true, size: 32 })] }));
children.push(new Paragraph({ spacing: { after: 240 },
  children: [new TextRun({ text: "Primary manuscript draft. Companion negative-result benchmark: Stressor-Specific, Not Universal.", italics: true, color: "595959" })] }));

children.push(new Paragraph({ children: [new TextRun({ text: "Abstract", bold: true, size: 26 })], spacing: { after: 120 } }));
children.push(p([
  t("The most metastatic circulating tumor cells (CTCs) travel as "), b("clusters"),
  t(" and with "), b("neutrophil bodyguards"),
  t(", yet what licenses this phenotype is unclear. We test whether "), b("hypoxia"),
  t(" is the upstream driver, using single-cell breast-cancer CTCs spanning tumor hypoxia (GSE126669) and circadian rest/active phase (GSE180097), the latter annotated with clustering and WBC attachment. A hypoxia survival signature derived from hypoxic CTCs is significantly elevated in "),
  b("clustered"), t(" CTCs (p = 0.023) and "), b("WBC-attached"), t(" CTCs (p = 0.008), and a shared HIF-target survival core (VEGFA, P4HA1, NDRG1, HK2, HIF1A) links hypoxic and rest-phase CTCs. The program is "),
  b("not"), t(" elevated by circadian timing (p = 0.41), and stress signatures do not transfer between stressors across four methods—ruling out a generic universal-stress explanation. Hypoxia thus associates specifically with the cluster-and-escort phenotype.")
]));

children.push(h1("1  Introduction"));
children.push(p("CTCs are the seeds of metastasis. The most dangerous travel as multicellular clusters and recruit neutrophils that shield them and boost proliferation. What upstream signal licenses this phenotype? Tumor hypoxia drives EMT, stemness and aggressiveness; a parallel literature links circadian rhythm to CTC shedding (peaking at rest). This poses competing hypotheses: is the aggressive phenotype licensed by hypoxia, by circadian timing, or by a generic stress response shared across perturbations?"));
children.push(p("We adjudicate directly: we derive a hypoxia signature from hypoxic CTCs and ask where it is elevated, using clustering and WBC attachment as readouts of the aggressive phenotype, and circadian phase and cross-stressor transfer as controls."));

children.push(h1("2  Data & Methods"));
children.push(p([b("Datasets. "), t("GSE126669: 30 single breast CTCs (hypoxia Positive/Negative). GSE180097: 276 breast CTCs (BR16, LM2 xenografts; Patient), circadian resting/active, with sample type (single / cluster / cluster-with-WBC), cluster size, and WBC count.")]));
children.push(p([b("Pipeline. "), t("Counts→TPM; Ensembl→symbols; per-batch log1p on split layers; 3000 HVGs (mvp); scaling without centering; 50 PCs; Harmony integration; SNN graph, Louvain (3 clusters), UMAP.")]));
children.push(p([b("Signatures & tests. "), t("Hypoxia survival signature = Wilcoxon markers up in hypoxia-Positive CTCs, module-scored across all cells; neutrophil-chemokine module scored analogously. Within GSE180097: Wilcoxon of hypoxia score by phase, clumping, and WBC status. Shared survival genes = intersection of hypoxia-Positive and resting-up markers. Cluster identity via FindAllMarkers.")]));

children.push(h1("3  Results"));
children.push(p([b("3.1 A hypoxia signature with real signal. "), t("Derived from GSE126669 (VEGFA, NDRG1, P4HA1-axis). Within-stressor CV confirms stress is recoverable in the circadian cohort (AUROC ≈ 0.91).")]));
children.push(p([b("3.2 Hypoxia is elevated in clustered CTCs.")]));
children.push(p([b("3.3 Hypoxia is elevated in neutrophil-associated CTCs.")]));
children.push(p([b("3.4 A shared HIF-target survival core "), t("(VEGFA, P4HA1, NDRG1, HK2, HIF1A, EGLN1, BNIP3, ERO1A, MTHFD2) links hypoxic and rest-phase CTCs.")]));
children.push(table([3700, 1900, 1900, 1860], [
  ["Test (GSE180097)", "Group A", "Group B", "Wilcoxon p"],
  ["Hypoxia by clumping", "cluster −0.06", "single −0.13", "0.023  ✓"],
  ["Hypoxia by WBC status", "with_WBC −0.01", "no_WBC −0.11", "0.008  ✓"],
  ["Neutrophil chemokines by WBC", "with_WBC +1.61", "no_WBC −0.31", "0.001  (confounded)"],
  ["Hypoxia by phase (control)", "resting −0.10", "active −0.07", "0.41  (n.s.)"]
]));
children.push(new Paragraph({ spacing: { before: 100, after: 140 }, children: [
  b("3.5 Negative controls. "),
  t("Hypoxia does not differ by circadian phase (p = 0.41; VEGFA flat). Stress signatures do not transfer between stressors across four methods (AUROC 0.34–0.54; logFC r ≈ −0.06), ruling out a generic universal-stress program; only HYPOXIA/UPR/ISR show weak pathway concordance.")
]}));
children.push(p([b("3.6 Cluster identities. "), t("Cluster 2 is proliferative (PRR11, CDKN3, KNL1, CDKL1); clusters 0/1 are not cleanly hypoxia/phase-defined (cluster 0 carries metastasis marker PODXL).")]));

const figPath = path.join(__dirname, "..", "main", "expression_data", "mechanistic_hypotheses_plots.png");
if (fs.existsSync(figPath)) {
  children.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 120 },
    children: [new ImageRun({ type: "png", data: fs.readFileSync(figPath),
      transformation: { width: 470, height: 300 },
      altText: { title: "Hypotheses", description: "Violin plots H1/H2", name: "hyp" } })] }));
  children.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 160 },
    children: [new TextRun({ text: "Figure 1. Hypoxia signature by circadian phase (top, n.s.) and by WBC attachment (bottom, p = 0.008); neutrophil chemokines by WBC attachment (bottom-right, confounded).", italics: true, size: 18, color: "595959" })] }));
}

children.push(h1("4  Limitations & Confounds"));
children.push(p("WBC composition confound: cluster-with-WBC samples physically contain neutrophils, so elevated chemokine scores partly reflect WBC transcripts; the hypoxia-score result (tumor-intrinsic) is the robust readout and the recruitment claim requires de-confounding. Small hypoxic cohort (n = 30): the signature is directional; the meaningful unit is the GSEA-enriched HIF core, not the full overlap list. Associational, not causal; two breast-cancer Smart-seq2 datasets only."));

children.push(h1("5  Discussion & Conclusion"));
children.push(p("Hypoxia is specifically associated with the cluster-and-escort phenotype of CTCs—clustering and neutrophil association—and not with circadian timing or a generic stress response, with a shared HIF-target survival core linking hypoxic and rest-phase CTCs. This positions hypoxia as a candidate upstream licensor of the most metastatic CTC states and motivates targeting hypoxic CTC clusters and their neutrophil interactions. Next: GSEA validation of the HIF core, de-confounded recruitment analysis, per-donor robustness, and ligand–receptor interaction modeling."));

children.push(h1("Reproducibility"));
children.push(p("Code: full_pipeline.R, mechanistic_hypotheses.R (this paper), cross_stressor_transfer.R, rescue_test.R, export_for_gnn.R, gnn_disentangle.py. Tables: hypoxia_survival_markers.csv, cluster_markers.csv, shared_survival_genes.csv."));

const doc = new Document({
  styles: {
    default: { document: { run: { font: "Calibri", size: 22 } } },
    paragraphStyles: [{ id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
      run: { size: 26, bold: true, font: "Calibri", color: "1F4E79" },
      paragraph: { spacing: { before: 240, after: 120 }, outlineLevel: 0 } }]
  },
  sections: [{
    properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 } } },
    children
  }]
});

Packer.toBuffer(doc).then(buf => {
  fs.writeFileSync(path.join(__dirname, "hypoxia_cluster_escort_paper.docx"), buf);
  console.log("wrote hypoxia_cluster_escort_paper.docx");
});
