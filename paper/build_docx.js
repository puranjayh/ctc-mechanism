const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  ImageRun, AlignmentType, HeadingLevel, BorderStyle, WidthType, ShadingType
} = require("docx");

const CW = 9360; // content width (US Letter, 1" margins)
const border = { style: BorderStyle.SINGLE, size: 1, color: "BFBFBF" };
const borders = { top: border, bottom: border, left: border, right: border };
const HEAD = "D5E8F0";

function h1(t) { return new Paragraph({ heading: HeadingLevel.HEADING_1, children: [new TextRun(t)] }); }
function p(runs) { return new Paragraph({ spacing: { after: 140 }, children: Array.isArray(runs) ? runs : [new TextRun(runs)] }); }
function t(text) { return new TextRun(text); }
function b(text) { return new TextRun({ text, bold: true }); }
function it(text) { return new TextRun({ text, italics: true }); }

function cell(text, w, opts = {}) {
  return new TableCell({
    borders, width: { size: w, type: WidthType.DXA },
    margins: { top: 80, bottom: 80, left: 120, right: 120 },
    shading: opts.head ? { fill: HEAD, type: ShadingType.CLEAR } : undefined,
    children: [new Paragraph({ children: [new TextRun({ text, bold: !!opts.head })] })]
  });
}
function table(widths, rows) {
  return new Table({
    width: { size: CW, type: WidthType.DXA }, columnWidths: widths,
    rows: rows.map((r, i) => new TableRow({
      children: r.map((c, j) => cell(String(c), widths[j], { head: i === 0 }))
    }))
  });
}

const children = [];

// Title
children.push(new Paragraph({
  spacing: { after: 80 },
  children: [new TextRun({ text: "Stressor-Specific, Not Universal: A Cross-Stressor Transfer Benchmark for Circulating Tumor Cells", bold: true, size: 32 })]
}));
children.push(new Paragraph({
  spacing: { after: 240 },
  children: [new TextRun({ text: "Workshop manuscript draft — target: NeurIPS LMRL / ML4H (negative-results & benchmarks track)", italics: true, color: "595959" })]
}));

// Abstract
children.push(new Paragraph({ children: [new TextRun({ text: "Abstract", bold: true, size: 26 })], spacing: { after: 120 } }));
children.push(p([
  t("Circulating tumor cells (CTCs) endure multiple physiological stresses in the bloodstream, motivating a popular assumption: that CTCs activate a "),
  it("conserved, universal stress program"),
  t(" usable as a liquid-biopsy biomarker. This assumption is rarely tested directly. We introduce a "),
  b("cross-stressor transfer benchmark"),
  t(" that operationalizes it as a falsifiable generalization question—does a stress signature learned from one stressor predict a biologically orthogonal stressor?—and a "),
  b("disentanglement graph neural network (GNN)"),
  t(" that separates shared from stressor-specific transcriptional axes. Applying four independent methods (sparse-linear transfer, genome-wide fold-change concordance, pathway-activity transfer, and the disentanglement GNN) to single-cell breast-cancer CTCs under hypoxic stress (GSE126669) and circadian-rest stress (GSE180097), we find that "),
  b("cross-stressor transfer does not exceed chance"),
  t(" (AUROC 0.34–0.54; concordance r ≈ −0.06). The CTC stress response is largely stressor-specific; the only shared component is a weak concordance confined to metabolic-stress pathways (hypoxia/UPR/ISR). A sensitivity analysis shows the result is robust to, though partly limited by, the small hypoxic cohort (n = 30). Our benchmark and convergent negative finding caution against single-signature stress biomarkers for CTCs and provide a reusable protocol for testing biological generalization in rare-cell regimes.")
]));

// 1 Introduction
children.push(h1("1  Introduction"));
children.push(p("CTCs are the seeds of metastasis and the substrate of liquid biopsy. In circulation they face hypoxia, shear stress, anchorage loss, oxidative stress and immune attack, and are hypothesized to mount adaptive stress responses promoting survival and dormancy. A natural corollary, increasingly assumed in the liquid-biopsy literature, is that these responses converge on a universal stress signature usable as a biomarker."));
children.push(p("We argue this corollary is a generalization claim and should be tested as one. If a universal stress program exists, a classifier trained on one stressor should transfer, zero-shot, to a different stressor. We turn two orthogonal perturbations—hypoxia and circadian rest phase—into source and target domains and ask whether stress transfers between them."));
children.push(p(b("Contributions.")));
children.push(new Paragraph({ numbering: { reference: "n", level: 0 }, children: [t("A "), b("cross-stressor transfer benchmark"), t(" for CTCs: a leave-one-stressor-out protocol with permutation and power controls, evaluated by four complementary methods.")] }));
children.push(new Paragraph({ numbering: { reference: "n", level: 0 }, children: [t("A "), b("disentanglement GNN"), t(" factoring each cell's embedding into a domain-invariant shared-stress axis and a stressor-specific axis.")] }));
children.push(new Paragraph({ numbering: { reference: "n", level: 0 }, children: [t("A "), b("convergent negative result"), t(": transfer fails across all four methods, overturning the universal-stress assumption and surfacing a weaker, pathway-restricted concordance for future work.")] }));

// 2 Related work
children.push(h1("2  Related Work"));
children.push(p([b("Universal stress responses. "), t("Conserved multi-stress programs are classic (e.g., the yeast Environmental Stress Response) and inspired literature-derived stress gene sets in toxicology and multi-stress classifiers in plants. These establish the concept but do not test cross-stressor transfer in CTCs.")]));
children.push(p([b("CTC stress and heterogeneity. "), t("Single-cell studies show CTCs activate stress-tolerance and immune-evasion programs and that hypoxia increases CTC aggressiveness; none test whether distinct stressors share a transferable program.")]));
children.push(p([b("Graph and transfer learning for single cells. "), t("GNN label transfer across single-cell datasets (scGCN) and CTC-specific transfer learning (CTC-Tracer) are established. We repurpose graph-transfer machinery for generalization testing rather than label imputation, adding explicit shared/specific disentanglement.")]));

// 3 Data
children.push(h1("3  Data"));
children.push(table([1100, 2200, 3360, 1200, 1500], [
  ["Domain", "Accession", "Stressor", "Cells", "Stressed/Baseline"],
  ["A", "GSE126669", "Hypoxia (Positive/Negative)", "30", "14 / 16"],
  ["B", "GSE180097", "Circadian (resting/active)", "276", "139 / 137"]
]));
children.push(new Paragraph({ spacing: { before: 100, after: 140 }, children: [t("Both are Smart-seq2 breast-cancer CTC datasets (GSE180097 spans BR16, LM2 xenograft and patient samples). We define a unified binary label where resting phase and hypoxia-positive = Stressed.")] }));

// 4 Methods
children.push(h1("4  Methods"));
children.push(p([b("Preprocessing & integration. "), t("Counts → TPM, Ensembl → symbols. To avoid the dominant batch driving feature selection, log1p was applied per batch on split layers; 3000 HVGs via mean-variance-plot; scaling without centering (for zero-inflated TPM); 50 PCs; batch integration with Harmony. A shared-nearest-neighbor (SNN) cell graph was built on the Harmony embedding.")]));
children.push(p([b("Transfer protocol. "), t("For each (source, target) we train on the source stressor's labels and evaluate zero-shot on the target (AUROC). Controls: within-domain cross-validation, a 200-shuffle label-permutation null, and a power/sensitivity analysis given n_A = 30.")]));
children.push(p([b("Methods. "), t("(1) L1-logistic transfer on HVG expression. (2) Per-gene Wilcoxon log2FC within each domain, correlating the two genome-wide vectors. (3) Module scoring of six curated stress programs and transfer on the 6-D pathway representation. (4) A two-layer GCN whose embedding splits into shared and specific blocks, trained with source-domain stress prediction, a gradient-reversal domain adversary on the shared block, and a shared/specific orthogonality penalty; only the shared block is used for transfer (10 seeds).")]));

// 5 Results
children.push(h1("5  Results"));
children.push(p([b("A within-stressor signal exists. "), t("Within-domain cross-validation recovers stress state in the larger circadian cohort (CV AUROC ≈ 0.91), confirming the labels carry real signal. The hypoxic cohort (n = 30) is too small for a stable within-domain estimate.")]));
children.push(p([b("Cross-stressor transfer fails across all four methods.")]));
children.push(table([3360, 2000, 2000, 2000], [
  ["Method", "Hyp→Circ", "Circ→Hyp", "Significance"],
  ["Sparse-linear (genes)", "0.504", "0.540", "perm p = 0.41 / 0.68"],
  ["Pathway-score transfer", "0.488", "0.446", "—"],
  ["Disentanglement GNN", "0.492 ± 0.08", "0.342 ± 0.11", "—"],
  ["Genome-wide log2FC corr.", "r = −0.059", "ρ = −0.070", "negligible effect"]
]));
children.push(new Paragraph({ spacing: { before: 100, after: 140 }, children: [t("No method exceeds chance; the GNN's Circadian→Hypoxia direction is below chance, echoing the slightly negative genome-wide concordance—the two responses are, if anything, mildly opposed. The sparse-linear model selected zero genes shared and sign-concordant across both domains.")] }));
children.push(p([b("The one shared component is pathway-restricted.")]));
children.push(table([3360, 2000, 2000, 2000], [
  ["Pathway", "Hyp shift", "Circ shift", "Concordant?"],
  ["HYPOXIA", "+0.519", "+0.136", "Yes"],
  ["UPR", "+0.009", "+0.196", "Yes"],
  ["ISR", "+0.006", "+0.148", "Yes"],
  ["P53", "−0.120", "+0.148", "No"],
  ["NF-κB", "+0.700", "−0.096", "No"],
  ["DORMANCY", "+0.220", "−0.038", "No"]
]));
children.push(new Paragraph({ spacing: { before: 100, after: 140 }, children: [t("A shared adaptive metabolic-stress core (hypoxia/UPR/ISR) moves the same direction under both stressors but is too weak to support transfer; inflammatory, p53 and dormancy responses are stressor-specific.")] }));

// Figure
const figPath = path.join(__dirname, "..", "main", "expression_data", "all_stress_conditions.png");
if (fs.existsSync(figPath)) {
  children.push(new Paragraph({
    alignment: AlignmentType.CENTER, spacing: { before: 120 },
    children: [new ImageRun({ type: "png", data: fs.readFileSync(figPath),
      transformation: { width: 460, height: 300 },
      altText: { title: "UMAP", description: "Harmony UMAP colored by stress condition", name: "umap" } })]
  }));
  children.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 160 },
    children: [new TextRun({ text: "Figure 1. Harmony-integrated UMAP. Stressed and baseline cells are intermixed across both stressors, consistent with the absence of a transferable shared stress axis.", italics: true, size: 18, color: "595959" })] }));
}

// 6 Limitations
children.push(h1("6  Limitations"));
children.push(p("The hypoxic cohort is small (n = 30): both training on it and predicting it (the unstable 0.342 ± 0.11) are underpowered, so part of the null reflects power, not only biology. The rigorous claim is therefore that, with available data, a shared transferable CTC stress axis is not detectable. We deliberately did not tune the GNN until a metric crossed an arbitrary threshold, to avoid optimistic bias. Both datasets are breast-cancer Smart-seq2; generality to other cancers/platforms is untested. Pathway concordance is hypothesis-generating, not confirmatory."));

// 7 Conclusion
children.push(h1("7  Conclusion"));
children.push(p("Across four independent methods, hypoxic and circadian-rest stress responses in CTCs do not share a transferable signature: stress in CTCs is stressor-specific, not universal, with only a weak pathway-restricted metabolic-stress overlap. We contribute a reusable cross-stressor transfer benchmark and a disentanglement GNN for probing biological generalization in rare-cell regimes. Our results caution against single-signature stress biomarkers and motivate larger matched multi-stressor cohorts and methods that model stressor-specific adaptation."));

// Reproducibility
children.push(h1("Reproducibility"));
children.push(p("Code: full_pipeline.R (preprocessing + Harmony), cross_stressor_transfer.R (Method 1), rescue_test.R (Methods 2–3), export_for_gnn.R + gnn_disentangle.py (Method 4). Result tables: concordant_stress_signature.csv."));

const doc = new Document({
  numbering: { config: [{ reference: "n", levels: [{ level: 0, format: "decimal", text: "%1.", alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 720, hanging: 360 } } } }] }] },
  styles: {
    default: { document: { run: { font: "Calibri", size: 22 } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 26, bold: true, font: "Calibri", color: "1F4E79" },
        paragraph: { spacing: { before: 240, after: 120 }, outlineLevel: 0 } }
    ]
  },
  sections: [{
    properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 } } },
    children
  }]
});

Packer.toBuffer(doc).then(buf => {
  fs.writeFileSync(path.join(__dirname, "cross_stressor_ctc_paper.docx"), buf);
  console.log("wrote cross_stressor_ctc_paper.docx");
});
