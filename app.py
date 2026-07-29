import streamlit as st
import requests
import pandas as pd
import time

st.set_page_config(page_title="Gene & Variant Explorer", layout="wide")

st.title("Gene & Variant Explorer")
st.write(
    "Look up a protein's sequence from UniProt and its known disease-linked "
    "variants from ClinVar. Export a FASTA file ready for structure prediction "
    "(e.g. ColabFold)."
)

# ---------------------------------------------------------------------------
# Core functions (same logic as the notebook version)
# ---------------------------------------------------------------------------

@st.cache_data(show_spinner=False)
def get_uniprot_info(gene_name, organism_id=9606):
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"gene:{gene_name} AND organism_id:{organism_id} AND reviewed:true",
        "format": "json",
        "fields": "accession,gene_names,protein_name,sequence,length",
        "size": 1,
    }
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    if not data.get("results"):
        raise ValueError(f"No reviewed UniProt entry found for gene '{gene_name}'.")

    entry = data["results"][0]
    return {
        "accession": entry["primaryAccession"],
        "protein_name": entry.get("proteinDescription", {})
                              .get("recommendedName", {})
                              .get("fullName", {})
                              .get("value", "unknown"),
        "gene_names": [g["geneName"]["value"] for g in entry.get("genes", []) if "geneName" in g],
        "sequence": entry["sequence"]["value"],
        "length": entry["sequence"]["length"],
    }


@st.cache_data(show_spinner=False)
def get_clinvar_variants(gene_name, retmax=50):
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    search_params = {
        "db": "clinvar",
        "term": f"{gene_name}[gene]",
        "retmax": retmax,
        "retmode": "json",
    }
    search_resp = requests.get(f"{base}/esearch.fcgi", params=search_params, timeout=30)
    search_resp.raise_for_status()
    id_list = search_resp.json()["esearchresult"].get("idlist", [])

    if not id_list:
        return pd.DataFrame()

    time.sleep(0.4)
    summary_params = {
        "db": "clinvar",
        "id": ",".join(id_list),
        "retmode": "json",
    }
    summary_resp = requests.get(f"{base}/esummary.fcgi", params=summary_params, timeout=30)
    summary_resp.raise_for_status()
    summary_data = summary_resp.json().get("result", {})

    rows = []
    for uid in summary_data.get("uids", []):
        rec = summary_data.get(uid, {})
        germline = rec.get("germline_classification", {})
        rows.append({
            "variant_name": rec.get("title"),
            "clinical_significance": germline.get("description"),
            "last_evaluated": germline.get("last_evaluated"),
            "conditions": "; ".join(
                t.get("trait_name", "") for t in germline.get("trait_set", [])
            ) if germline.get("trait_set") else None,
        })

    return pd.DataFrame(rows)


def to_fasta(accession, gene_name, sequence):
    lines = [f">{accession}_{gene_name}"]
    for i in range(0, len(sequence), 60):
        lines.append(sequence[i:i + 60])
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

col1, col2 = st.columns([3, 1])
with col1:
    gene_name = st.text_input("Gene name", value="CFTR")
with col2:
    st.write("")
    st.write("")
    run = st.button("Run", type="primary")

if run and gene_name:
    try:
        with st.spinner("Fetching protein info from UniProt..."):
            protein_info = get_uniprot_info(gene_name)

        st.subheader("Protein info")
        c1, c2, c3 = st.columns(3)
        c1.metric("Accession", protein_info["accession"])
        c2.metric("Length", f"{protein_info['length']} aa")
        c3.metric("Gene name(s)", ", ".join(protein_info["gene_names"]))
        st.caption(protein_info["protein_name"])
        st.code(protein_info["sequence"], language=None)

        with st.spinner("Fetching variants from ClinVar..."):
            variants_df = get_clinvar_variants(gene_name)

        st.subheader(f"ClinVar variants ({len(variants_df)} found)")
        if not variants_df.empty:
            st.dataframe(variants_df, use_container_width=True)

            pathogenic_mask = variants_df["clinical_significance"].str.contains(
                "pathogenic", case=False, na=False
            )
            pathogenic_df = variants_df[pathogenic_mask]
            st.subheader(f"Pathogenic variants ({len(pathogenic_df)})")
            st.dataframe(pathogenic_df, use_container_width=True)

            csv_bytes = variants_df.to_csv(index=False).encode("utf-8")
            st.download_button(
                "Download variant table (CSV)",
                data=csv_bytes,
                file_name=f"{gene_name}_clinvar_variants.csv",
                mime="text/csv",
            )
        else:
            st.info("No ClinVar variants found for this gene.")

        fasta_text = to_fasta(protein_info["accession"], gene_name, protein_info["sequence"])
        st.download_button(
            "Download sequence (FASTA)",
            data=fasta_text,
            file_name=f"{gene_name}_{protein_info['accession']}.fasta",
            mime="text/plain",
        )

    except Exception as e:
        st.error(f"Something went wrong: {e}")
