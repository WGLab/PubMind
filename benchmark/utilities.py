import pandas as pd
import re
import os
import glob
import csv
import xml.etree.ElementTree as ET
from typing import List, Dict, Tuple, Optional

##########################################################################
########helper function to load different format other than csv/tsv#######
##########################################################################
def read_pubtator_format(path):
    """Input file with PubTator format, output pd dataframe with all entities with PMID."""
    
    with open(path, 'r') as file:
        corpus = file.readlines()

    # Regular expression to identify rows starting with a PubMed ID followed by tab-separated values
    pattern = r'^\d+\t.+'

    # Filter rows matching the pattern
    table_rows = [line.strip() for line in corpus if re.match(pattern, line)]
    data = [row.split('\t') for row in table_rows]
    df = pd.DataFrame(data)
    #df.columns=['PMID','start','end','Mutation','Type']
    return df

def load_seth_corpus_xml_labels(
    xml_path: str,
    sep: str = "\n",
    keep_doc_text: bool = False,
) -> Tuple[pd.DataFrame, Optional[pd.DataFrame]]:
    """
    Parse SETH/resources/.../corpus.xml and extract <gene> and <variant> mentions with offsets.

    Offsets are computed on doc_text = Title_text + sep + Abstract_text (if both exist).

    Returns:
      labels_df: DataFrame with columns:
        PMID, start, end, text, Type, section, source, attrs...
      docs_df (optional): DataFrame with PMID, title, abstract, doc_text
    """

    tree = ET.parse(xml_path)
    root = tree.getroot()

    label_rows: List[Dict] = []
    doc_rows: List[Dict] = []

    def _reconstruct_with_spans(node: ET.Element, pmid: str, section: str, base_offset: int) -> Tuple[str, List[Dict]]:
        """
        Reconstruct plain text from an element containing mixed content (text + child tags),
        while recording spans for child tags <gene> and <variant>.
        Offsets are relative to doc_text, so we use base_offset + running length.
        """
        parts: List[str] = []
        spans: List[Dict] = []
        cur_len = 0

        # node.text occurs before any first child
        if node.text:
            parts.append(node.text)
            cur_len += len(node.text)

        for child in list(node):
            # child element text content (the mention itself)
            mention_text = child.text or ""
            tag = child.tag  # "gene" or "variant" (expected)

            # record span if gene/variant
            if tag in ("gene", "variant") and mention_text:
                start = base_offset + cur_len
                end = start + len(mention_text)
                row = {
                    "PMID": pmid,
                    "start": start,
                    "end": end,
                    "text": mention_text,
                    "Type": "Gene" if tag == "gene" else "Variant",
                    "section": section,
                    "source": "SETH_XML",
                }
                # include attributes (g_id/g_lex/v_lex/v_norm etc.)
                for k, v in child.attrib.items():
                    row[k] = v
                spans.append(row)

            # append the mention text to reconstructed plain text
            parts.append(mention_text)
            cur_len += len(mention_text)

            # child.tail occurs after the child element before next sibling
            if child.tail:
                parts.append(child.tail)
                cur_len += len(child.tail)

        return "".join(parts), spans

    # The file structure is <Articles><Article>...
    for art in root.findall(".//Article"):
        pmid_el = art.find("Pmid")
        if pmid_el is None or pmid_el.text is None:
            continue
        pmid = pmid_el.text.strip()

        title_el = art.find("Title")
        abstract_el = art.find("Abstract")

        title_text = ""
        abstract_text = ""
        title_spans: List[Dict] = []
        abstract_spans: List[Dict] = []

        # Title part: base_offset = 0
        if title_el is not None:
            title_text, title_spans = _reconstruct_with_spans(
                node=title_el,
                pmid=pmid,
                section="Title",
                base_offset=0
            )

        # Abstract part: base_offset = len(title_text) + len(sep) (if title exists), else 0
        if abstract_el is not None:
            base = len(title_text) + (len(sep) if title_text else 0)
            abstract_text, abstract_spans = _reconstruct_with_spans(
                node=abstract_el,
                pmid=pmid,
                section="Abstract",
                base_offset=base
            )

        # finalize doc_text (used for offset reference)
        doc_text = title_text + (sep + abstract_text if title_text and abstract_text else abstract_text)

        # collect spans
        label_rows.extend(title_spans)
        label_rows.extend(abstract_spans)

        if keep_doc_text:
            doc_rows.append({
                "PMID": pmid,
                "title": title_text,
                "abstract": abstract_text,
                "doc_text": doc_text
            })

    labels_df = pd.DataFrame(label_rows).sort_values(["PMID", "start", "end"]).reset_index(drop=True)
    docs_df = pd.DataFrame(doc_rows) if keep_doc_text else None
    return labels_df, docs_df

######################################################################################
########helper function for comparison with tools' reulst with Ground Truth(GT)#######
######################################################################################

####################################
######DNA
####################################

def clean_dna_mut_for_benchmark(df):
    df_dna = df.copy()
    df_dna['Mutation']=df_dna['Mutation'].apply(lambda x: x[2:] if x[:2]=='c.'else x)
    df_dna=df_dna[df_dna['Mutation'].str.contains(r'\d', na=False)]
    #remove the space, _, () for better comparison
    df_dna['Mutation']=df_dna['Mutation'].str.replace(r"[ _()-]", "", regex=True)
    df_dna['Mutation']=df_dna['Mutation'].str.replace(r"incDNA", "", regex=True)
    df_dna['Mutation']=df_dna['Mutation'].str.replace(r"by", "", regex=True)
    df_dna['Mutation']=df_dna['Mutation'].str.replace(r"promoter", "", regex=True)
    return df_dna

def dna_benchmark_tool_with_gt(
    df_tool_dna: pd.DataFrame,
    df_gt_dna: pd.DataFrame,
    tool: str,
    gt_col: str = "Mutation",
    tool_col: str = "Mutation",
    match_mode: str = "substring",   # "exact" or "substring"
    macro: bool = True,
    verbose: bool = True
):
    """
    Compare tool vs GT for DNA variants using per-PMID sets.

    Returns:
      merged_df: per-PMID sets + TP/FP/FN + per-PMID metrics
      summary: dict of micro + macro metrics

    Metrics:
      precision = TP / (TP + FP)
      recall    = TP / (TP + FN)
      f1        = 2PR/(P+R)
      accuracy  = Jaccard = TP / (TP + FP + FN)   (set-based accuracy)
    """

    tool_mut_col = f"{tool}_Mutations"
    tool_total_col = f"Total_{tool}"
    
    df_gt_dna=df_gt_dna.drop_duplicates()
    df_tool_dna=df_tool_dna.drop_duplicates()
    
    # ---- group into sets per PMID ----
    gt_grouped = df_gt_dna.groupby("PMID")[gt_col].apply(lambda x: set(x.dropna())).reset_index(name="GT_Mutations")
    tool_grouped = df_tool_dna.groupby("PMID")[tool_col].apply(lambda x: set(x.dropna())).reset_index(name=tool_mut_col)

    gt_grouped["PMID"] = gt_grouped["PMID"].astype(str)
    tool_grouped["PMID"] = tool_grouped["PMID"].astype(str)

    merged = pd.merge(gt_grouped, tool_grouped, on="PMID", how="outer")
    merged["GT_Mutations"] = merged["GT_Mutations"].apply(lambda x: x if isinstance(x, set) else set())
    merged[tool_mut_col] = merged[tool_mut_col].apply(lambda x: x if isinstance(x, set) else set())

    # exclude PMIDs where GT is empty
    merged = merged[merged["GT_Mutations"].apply(len) > 0].copy()

    # ---- matching functions ----
    def _hit(gt_item: str, tool_item: str) -> bool:
        if match_mode == "exact":
            return gt_item == tool_item
        # substring mode: treat GT as matched if GT string occurs in tool string OR vice versa
        # (robust when one has extra context)
        return (gt_item in tool_item) or (tool_item in gt_item)

    def _count_tp_fp_fn(gt_set: set, tool_set: set):
        # TP: how many GT items are matched by any tool item
        tp = 0
        for g in gt_set:
            if any(_hit(str(g), str(t)) for t in tool_set):
                tp += 1

        fn = len(gt_set) - tp

        # FP: tool items that do not match any GT item
        fp = 0
        for t in tool_set:
            if not any(_hit(str(g), str(t)) for g in gt_set):
                fp += 1

        return tp, fp, fn

    # ---- compute TP/FP/FN per PMID ----
    merged[["TP", "FP", "FN"]] = merged.apply(
        lambda row: pd.Series(_count_tp_fp_fn(row["GT_Mutations"], row[tool_mut_col])),
        axis=1
    )

    merged["Total_GT"] = merged["GT_Mutations"].apply(len)
    merged[tool_total_col] = merged[tool_mut_col].apply(len)

    # per-PMID metrics
    def _prf(tp, fp, fn):
        p = tp / (tp + fp) if (tp + fp) else 0.0
        r = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = (2*p*r/(p+r)) if (p+r) else 0.0
        acc = tp / (tp + fp + fn) if (tp + fp + fn) else 0.0  # Jaccard accuracy
        return p, r, f1, acc

    merged[["Precision", "Recall", "F1", "Accuracy"]] = merged.apply(
        lambda row: pd.Series(_prf(row["TP"], row["FP"], row["FN"])),
        axis=1
    )

    # ---- micro aggregation ----
    TP = int(merged["TP"].sum())
    FP = int(merged["FP"].sum())
    FN = int(merged["FN"].sum())

    micro_p, micro_r, micro_f1, micro_acc = _prf(TP, FP, FN)

    summary = {
        "tool": tool,
        "match_mode": match_mode,
        "micro": {
            "TP": TP, "FP": FP, "FN": FN,
            "precision": micro_p,
            "recall": micro_r,
            "f1": micro_f1,
            "accuracy_jaccard": micro_acc,
        }
    }

    # ---- macro aggregation (optional) ----
    if macro:
        df_nonempty = merged[merged["Total_GT"] > 0].copy()
        summary["macro"] = {
            "precision": float(df_nonempty["Precision"].mean()) if len(df_nonempty) else 0.0,
            "recall": float(df_nonempty["Recall"].mean()) if len(df_nonempty) else 0.0,
            "f1": float(df_nonempty["F1"].mean()) if len(df_nonempty) else 0.0,
            "accuracy_jaccard": float(df_nonempty["Accuracy"].mean()) if len(df_nonempty) else 0.0,
        }

    if verbose:
        m = summary["micro"]
        print(
            f"## {tool} benchmark ({match_mode}) | "
            f"TP={m['TP']} FP={m['FP']} FN={m['FN']} | "
            f"P={m['precision']:.4f} R={m['recall']:.4f} F1={m['f1']:.4f} Acc(Jacc)={m['accuracy_jaccard']:.4f}"
        )

    # return per-PMID details + summary
    keep_cols = ["PMID", "GT_Mutations", tool_mut_col, "TP", "FP", "FN", "Total_GT", tool_total_col,
                 "Precision", "Recall", "F1", "Accuracy"]
    return merged[keep_cols], summary


####################################
#######protein
####################################
def replace_amino_acids(text):
    def replace_match(match):
        return new_three_to_one_letter_map.get(match.group(0), match.group(0))

    # Replace three-letter amino acid codes found in the string
    # should allow wrong upper/lower case spelling, e.g., {R753Q, R753GLn}	{R753Q}
    return re.sub(r'[A-Za-z]{3}', replace_match, text)

def clean_pro_mut_for_benchmark(df):
    df_tmVar_pro = df.copy()
    df_tmVar_pro.loc[:, 'Mutation']=df_tmVar_pro['Mutation'].apply(lambda x: x[2:] if x[:2]==r'p.'else x)
    df_tmVar_pro=df_tmVar_pro[df_tmVar_pro['Mutation'].str.contains(r'\d', na=False)]
    #df_tmVar_pro['Mutation']=df_tmVar_pro['Mutation'].apply(replace_amino_acids)
    df_tmVar_pro['Mutation']=df_tmVar_pro['Mutation'].str.replace(r"[ _()-]", "", regex=True)
    df_tmVar_pro['Mutation']=df_tmVar_pro['Mutation'].str.replace(r"by", "", regex=True)
    df_tmVar_pro['Mutation']=df_tmVar_pro['Mutation'].str.replace(r"to", "", regex=True)
    df_tmVar_pro['Mutation']=df_tmVar_pro['Mutation'].str.replace(r"stop", "*", regex=True)
    return df_tmVar_pro

def protein_benchmark_tool_with_gt(
    df_tool_pro: pd.DataFrame,
    df_gt_pro: pd.DataFrame,
    tool: str,
    gt_col: str = "Mutation",
    tool_col: str = "Mutation",
    match_mode: str = "substring",   # "exact" or "substring"
    macro: bool = True,
    verbose: bool = True
):
    """
    Compare tool vs GT for protein variants using per-PMID sets.

    Returns:
      detail_df: per-PMID sets + TP/FP/FN + per-PMID metrics
      summary: dict of micro + (optional) macro metrics

    Metrics:
      precision = TP / (TP + FP)
      recall    = TP / (TP + FN)
      f1        = 2PR/(P+R)
      accuracy  = Jaccard = TP / (TP + FP + FN)   (set-based accuracy)
    """

    tool_mut_col = f"{tool}_Mutations"
    tool_total_col = f"Total_{tool}"
    
    df_gt_pro=df_gt_pro.drop_duplicates()
    df_tool_pro=df_tool_pro.drop_duplicates()
    
    # ---- group into sets per PMID ----
    gt_grouped = df_gt_pro.groupby("PMID")[gt_col].apply(lambda x: set(x.dropna())).reset_index(name="GT_Mutations")
    tool_grouped = df_tool_pro.groupby("PMID")[tool_col].apply(lambda x: set(x.dropna())).reset_index(name=tool_mut_col)

    gt_grouped["PMID"] = gt_grouped["PMID"].astype(str)
    tool_grouped["PMID"] = tool_grouped["PMID"].astype(str)

    merged = pd.merge(gt_grouped, tool_grouped, on="PMID", how="outer")
    merged["GT_Mutations"] = merged["GT_Mutations"].apply(lambda x: x if isinstance(x, set) else set())
    merged[tool_mut_col] = merged[tool_mut_col].apply(lambda x: x if isinstance(x, set) else set())

    # exclude PMIDs where GT is empty
    merged = merged[merged["GT_Mutations"].apply(len) > 0].copy()

    # ---- matching ----
    def _hit(gt_item: str, tool_item: str) -> bool:
        if match_mode == "exact":
            return gt_item == tool_item
        # substring (relaxed): treat as match if either contains the other
        # robust when tool adds context or formatting differences remain
        return (gt_item in tool_item) or (tool_item in gt_item)

    def _count_tp_fp_fn(gt_set: set, tool_set: set):
        # TP: #GT variants matched by any tool variant
        tp = 0
        for g in gt_set:
            if any(_hit(str(g), str(t)) for t in tool_set):
                tp += 1
        fn = len(gt_set) - tp

        # FP: #tool variants that match NO GT variant
        fp = 0
        for t in tool_set:
            if not any(_hit(str(g), str(t)) for g in gt_set):
                fp += 1

        return tp, fp, fn

    merged[["TP", "FP", "FN"]] = merged.apply(
        lambda row: pd.Series(_count_tp_fp_fn(row["GT_Mutations"], row[tool_mut_col])),
        axis=1
    )

    merged["Total_GT"] = merged["GT_Mutations"].apply(len)
    merged[tool_total_col] = merged[tool_mut_col].apply(len)

    # ---- per-PMID metrics ----
    def _prf(tp, fp, fn):
        p = tp / (tp + fp) if (tp + fp) else 0.0
        r = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = (2*p*r/(p+r)) if (p+r) else 0.0
        acc = tp / (tp + fp + fn) if (tp + fp + fn) else 0.0  # Jaccard
        return p, r, f1, acc

    merged[["Precision", "Recall", "F1", "Accuracy"]] = merged.apply(
        lambda row: pd.Series(_prf(row["TP"], row["FP"], row["FN"])),
        axis=1
    )

    # ---- micro ----
    TP = int(merged["TP"].sum())
    FP = int(merged["FP"].sum())
    FN = int(merged["FN"].sum())
    micro_p, micro_r, micro_f1, micro_acc = _prf(TP, FP, FN)

    summary = {
        "tool": tool,
        "match_mode": match_mode,
        "micro": {
            "TP": TP, "FP": FP, "FN": FN,
            "precision": micro_p,
            "recall": micro_r,
            "f1": micro_f1,
            "accuracy_jaccard": micro_acc,
        }
    }

    # ---- macro ----
    if macro:
        df_nonempty = merged[merged["Total_GT"] > 0].copy()
        summary["macro"] = {
            "precision": float(df_nonempty["Precision"].mean()) if len(df_nonempty) else 0.0,
            "recall": float(df_nonempty["Recall"].mean()) if len(df_nonempty) else 0.0,
            "f1": float(df_nonempty["F1"].mean()) if len(df_nonempty) else 0.0,
            "accuracy_jaccard": float(df_nonempty["Accuracy"].mean()) if len(df_nonempty) else 0.0,
        }

    if verbose:
        m = summary["micro"]
        print(
            f"## {tool} benchmark ({match_mode}) | "
            f"TP={m['TP']} FP={m['FP']} FN={m['FN']} | "
            f"P={m['precision']:.4f} R={m['recall']:.4f} F1={m['f1']:.4f} Acc(Jacc)={m['accuracy_jaccard']:.4f}"
        )

    keep_cols = ["PMID", "GT_Mutations", tool_mut_col, "TP", "FP", "FN",
                 "Total_GT", tool_total_col, "Precision", "Recall", "F1", "Accuracy"]
    return merged[keep_cols], summary

####################################
######DNA and protein teogether 
####################################
#(some corpus or tools only provide "Variant" level annotation without specifying DNA/protein mutation.
def clean_dnaNpro_mut_for_benchmark(df):
    df_dna = df.copy()
    df_dna['Mutation']=df_dna['Mutation'].apply(lambda x: x[2:] if x[:2]=='c.'else x)
    df_dna.loc[:, 'Mutation']=df_dna['Mutation'].apply(lambda x: x[2:] if x[:2]==r'p.'else x)
    df_dna=df_dna[df_dna['Mutation'].str.contains(r'\d', na=False)]
    #remove the space, _, () for better comparison
    df_dna['Mutation']=df_dna['Mutation'].str.replace(r"[ _()-]", "", regex=True)
    df_dna['Mutation']=df_dna['Mutation'].str.replace(r"incDNA", "", regex=True)
    df_dna['Mutation']=df_dna['Mutation'].str.replace(r"by", "", regex=True)
    df_dna['Mutation']=df_dna['Mutation'].str.replace(r"promoter", "", regex=True)
    df_dna['Mutation']=df_dna['Mutation'].str.replace(r"to", "", regex=True)
    df_dna['Mutation']=df_dna['Mutation'].str.replace(r"stop", "*", regex=True)
    return df_dna

def dnaNpro_benchmark_tool_with_gt(
    df_tool_dna: pd.DataFrame,
    df_gt_dna: pd.DataFrame,
    tool: str,
    gt_col: str = "Mutation",
    tool_col: str = "Mutation",
    match_mode: str = "substring",   # "exact" or "substring"
    macro: bool = True,
    verbose: bool = True
):
    """
    Compare tool vs GT for DNA variants using per-PMID sets.

    Returns:
      merged_df: per-PMID sets + TP/FP/FN + per-PMID metrics
      summary: dict of micro + macro metrics

    Metrics:
      precision = TP / (TP + FP)
      recall    = TP / (TP + FN)
      f1        = 2PR/(P+R)
      accuracy  = Jaccard = TP / (TP + FP + FN)   (set-based accuracy)
    """

    tool_mut_col = f"{tool}_Mutations"
    tool_total_col = f"Total_{tool}"

    df_gt_dna=df_gt_dna.drop_duplicates()
    df_tool_dna=df_tool_dna.drop_duplicates()
    
    # ---- group into sets per PMID ----
    gt_grouped = df_gt_dna.groupby("PMID")[gt_col].apply(lambda x: set(x.dropna())).reset_index(name="GT_Mutations")
    tool_grouped = df_tool_dna.groupby("PMID")[tool_col].apply(lambda x: set(x.dropna())).reset_index(name=tool_mut_col)

    gt_grouped["PMID"] = gt_grouped["PMID"].astype(str)
    tool_grouped["PMID"] = tool_grouped["PMID"].astype(str)

    merged = pd.merge(gt_grouped, tool_grouped, on="PMID", how="outer")
    merged["GT_Mutations"] = merged["GT_Mutations"].apply(lambda x: x if isinstance(x, set) else set())
    merged[tool_mut_col] = merged[tool_mut_col].apply(lambda x: x if isinstance(x, set) else set())

    # exclude PMIDs where GT is empty
    merged = merged[merged["GT_Mutations"].apply(len) > 0].copy()

    # ---- matching functions ----
    def _hit(gt_item: str, tool_item: str) -> bool:
        if match_mode == "exact":
            return gt_item == tool_item
        # substring mode: treat GT as matched if GT string occurs in tool string OR vice versa
        # (robust when one has extra context)
        return (gt_item in tool_item) or (tool_item in gt_item)

    def _count_tp_fp_fn(gt_set: set, tool_set: set):
        # TP: how many GT items are matched by any tool item
        tp = 0
        for g in gt_set:
            if any(_hit(str(g), str(t)) for t in tool_set):
                tp += 1

        fn = len(gt_set) - tp

        # FP: tool items that do not match any GT item
        fp = 0
        for t in tool_set:
            if not any(_hit(str(g), str(t)) for g in gt_set):
                fp += 1

        return tp, fp, fn

    # ---- compute TP/FP/FN per PMID ----
    merged[["TP", "FP", "FN"]] = merged.apply(
        lambda row: pd.Series(_count_tp_fp_fn(row["GT_Mutations"], row[tool_mut_col])),
        axis=1
    )

    merged["Total_GT"] = merged["GT_Mutations"].apply(len)
    merged[tool_total_col] = merged[tool_mut_col].apply(len)

    # per-PMID metrics
    def _prf(tp, fp, fn):
        p = tp / (tp + fp) if (tp + fp) else 0.0
        r = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = (2*p*r/(p+r)) if (p+r) else 0.0
        acc = tp / (tp + fp + fn) if (tp + fp + fn) else 0.0  # Jaccard accuracy
        return p, r, f1, acc

    merged[["Precision", "Recall", "F1", "Accuracy"]] = merged.apply(
        lambda row: pd.Series(_prf(row["TP"], row["FP"], row["FN"])),
        axis=1
    )

    # ---- micro aggregation ----
    TP = int(merged["TP"].sum())
    FP = int(merged["FP"].sum())
    FN = int(merged["FN"].sum())

    micro_p, micro_r, micro_f1, micro_acc = _prf(TP, FP, FN)

    summary = {
        "tool": tool,
        "match_mode": match_mode,
        "micro": {
            "TP": TP, "FP": FP, "FN": FN,
            "precision": micro_p,
            "recall": micro_r,
            "f1": micro_f1,
            "accuracy_jaccard": micro_acc,
        }
    }

    # ---- macro aggregation (optional) ----
    if macro:
        df_nonempty = merged[merged["Total_GT"] > 0].copy()
        summary["macro"] = {
            "precision": float(df_nonempty["Precision"].mean()) if len(df_nonempty) else 0.0,
            "recall": float(df_nonempty["Recall"].mean()) if len(df_nonempty) else 0.0,
            "f1": float(df_nonempty["F1"].mean()) if len(df_nonempty) else 0.0,
            "accuracy_jaccard": float(df_nonempty["Accuracy"].mean()) if len(df_nonempty) else 0.0,
        }

    if verbose:
        m = summary["micro"]
        print(
            f"## {tool} benchmark ({match_mode}) | "
            f"TP={m['TP']} FP={m['FP']} FN={m['FN']} | "
            f"P={m['precision']:.4f} R={m['recall']:.4f} F1={m['f1']:.4f} Acc(Jacc)={m['accuracy_jaccard']:.4f}"
        )

    # return per-PMID details + summary
    keep_cols = ["PMID", "GT_Mutations", tool_mut_col, "TP", "FP", "FN", "Total_GT", tool_total_col,
                 "Precision", "Recall", "F1", "Accuracy"]
    return merged[keep_cols], summary



####################################
######DNA or protein benchmark using the manual adjusted result 
####################################

def dnaorpro_benchmark_from_counts(
    df_reviewed: pd.DataFrame,
    tool: str,
    tp_col: str = "TP",
    fp_col: str = "FP",
    fn_col: str = "FN",
    pmid_col: str = "PMID",
    macro: bool = True,
    verbose: bool = True,
):
    """
    Compute benchmark metrics directly from a reviewed table that already contains
    TP / FP / FN per PMID.

    Returns:
      reviewed_df: original table + per-PMID Precision/Recall/F1/Accuracy
      summary: dict with micro + macro metrics

    Metrics:
      precision = TP / (TP + FP)
      recall    = TP / (TP + FN)
      f1        = 2PR / (P + R)
      accuracy  = Jaccard = TP / (TP + FP + FN)
    """

    df = df_reviewed.copy()

    # make sure columns exist
    required_cols = [pmid_col, tp_col, fp_col, fn_col]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # numeric cleanup
    for col in [tp_col, fp_col, fn_col]:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)

    # optional total GT and total tool counts reconstructed from TP/FP/FN
    df["Total_GT"] = df[tp_col] + df[fn_col]
    df[f"Total_{tool}"] = df[tp_col] + df[fp_col]

    # per-PMID metrics
    def _prf(tp, fp, fn):
        p = tp / (tp + fp) if (tp + fp) else 0.0
        r = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = (2 * p * r / (p + r)) if (p + r) else 0.0
        acc = tp / (tp + fp + fn) if (tp + fp + fn) else 0.0
        return pd.Series([p, r, f1, acc], index=["Precision", "Recall", "F1", "Accuracy"])

    df[["Precision", "Recall", "F1", "Accuracy"]] = df.apply(
        lambda row: _prf(row[tp_col], row[fp_col], row[fn_col]),
        axis=1
    )

    # ---- micro aggregation ----
    TP = int(df[tp_col].sum())
    FP = int(df[fp_col].sum())
    FN = int(df[fn_col].sum())

    micro_p = TP / (TP + FP) if (TP + FP) else 0.0
    micro_r = TP / (TP + FN) if (TP + FN) else 0.0
    micro_f1 = (2 * micro_p * micro_r / (micro_p + micro_r)) if (micro_p + micro_r) else 0.0
    micro_acc = TP / (TP + FP + FN) if (TP + FP + FN) else 0.0

    summary = {
        "tool": tool,
        "match_mode": "substring+manual_adjusted",
        "micro": {
            "TP": TP,
            "FP": FP,
            "FN": FN,
            "precision": micro_p,
            "recall": micro_r,
            "f1": micro_f1,
            "accuracy_jaccard": micro_acc,
        }
    }

    # ---- macro aggregation ----
    if macro:
        df_nonempty = df[df["Total_GT"] > 0].copy()
        summary["macro"] = {
            "precision": float(df_nonempty["Precision"].mean()) if len(df_nonempty) else 0.0,
            "recall": float(df_nonempty["Recall"].mean()) if len(df_nonempty) else 0.0,
            "f1": float(df_nonempty["F1"].mean()) if len(df_nonempty) else 0.0,
            "accuracy_jaccard": float(df_nonempty["Accuracy"].mean()) if len(df_nonempty) else 0.0,
        }

    if verbose:
        m = summary["micro"]
        print(
            f"## {tool} benchmark (from reviewed TP/FP/FN) | "
            f"TP={m['TP']} FP={m['FP']} FN={m['FN']} | "
            f"P={m['precision']:.4f} R={m['recall']:.4f} "
            f"F1={m['f1']:.4f} Acc(Jacc)={m['accuracy_jaccard']:.4f}"
        )

    return df, summary
