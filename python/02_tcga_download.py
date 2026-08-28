"""Package: pandas. Prerequisite: AURORA candidate CSV and internet access to cBioPortal."""

import json
import time
import urllib.parse
import urllib.request

import pandas as pd

from common import AURORA_DIR, TCGA_RAW_DIR, ensure_output_dirs, require_columns


API_BASE = "https://www.cbioportal.org/api"


def get_json(url, timeout=60):
    request = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def post_json(url, body, timeout=120):
    request = urllib.request.Request(
        url,
        data=json.dumps(body).encode("utf-8"),
        method="POST",
        headers={"Accept": "application/json", "Content-Type": "application/json"},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def paginated_clinical(clinical_type, attribute_id=None):
    records = []
    page = 0
    while True:
        query = {
            "clinicalDataType": clinical_type,
            "projection": "DETAILED",
            "pageSize": 2000,
            "pageNumber": page,
        }
        if attribute_id:
            query["clinicalAttributeId"] = attribute_id
        url = f"{API_BASE}/studies/brca_tcga/clinical-data?{urllib.parse.urlencode(query)}"
        batch = get_json(url)
        records.extend(batch)
        if len(batch) < 2000:
            return records
        page += 1
        time.sleep(0.5)


def main():
    ensure_output_dirs()
    candidates = pd.read_csv(AURORA_DIR / "aurora_candidate_genes.csv")
    require_columns(candidates, ["gene"], AURORA_DIR / "aurora_candidate_genes.csv")
    genes = candidates.gene.tolist()
    if len(genes) != 10:
        raise ValueError(f"Expected 10 AURORA candidates, found {len(genes)}")

    clinical_long = paginated_clinical("PATIENT")
    clinical_records = {}
    for row in clinical_long:
        patient_id = row.get("patientId", "")
        clinical_records.setdefault(patient_id, {})[row.get("clinicalAttributeId", "")] = row.get("value", "")
    clinical = pd.DataFrame.from_dict(clinical_records, orient="index")
    clinical.index.name = "patientId"
    clinical.reset_index().to_csv(TCGA_RAW_DIR / "tcga_brca_clinical_raw.csv", index=False)

    gene_ids = {}
    for gene in genes:
        matches = get_json(f"{API_BASE}/genes?keyword={urllib.parse.quote(gene)}&pageSize=10")
        exact = [item for item in matches if item.get("hugoGeneSymbol") == gene]
        if not exact:
            raise ValueError(f"cBioPortal did not return an exact gene match for {gene}")
        gene_ids[gene] = exact[0]["entrezGeneId"]
        time.sleep(0.2)

    mutation_url = (
        f"{API_BASE}/molecular-profiles/brca_tcga_mutations/mutations/fetch"
        "?projection=DETAILED&pageSize=5000&pageNumber=0"
    )
    mutation_json = post_json(
        mutation_url,
        {"sampleListId": "brca_tcga_sequenced", "entrezGeneIds": list(gene_ids.values())},
    )
    mutation_rows = []
    for item in mutation_json:
        mutation_rows.append(
            {
                "sampleId": item.get("sampleId", ""),
                "patientId": item.get("patientId", ""),
                "gene": item.get("gene", {}).get("hugoGeneSymbol", ""),
                "entrezGeneId": item.get("gene", {}).get("entrezGeneId", ""),
                "mutationType": item.get("mutationType", ""),
                "proteinChange": item.get("proteinChange", ""),
                "mutationStatus": item.get("mutationStatus", ""),
            }
        )
    mutations = pd.DataFrame(mutation_rows)
    mutations = mutations[mutations.mutationStatus.eq("Somatic")]
    mutations.to_csv(TCGA_RAW_DIR / "tcga_brca_mutations_raw.csv", index=False)

    tmb_long = paginated_clinical("SAMPLE", "TMB_NONSYNONYMOUS")
    tmb = pd.DataFrame(
        [
            {"sampleId": row.get("sampleId", ""), "tmb": pd.to_numeric(row.get("value", ""), errors="coerce")}
            for row in tmb_long
            if row.get("clinicalAttributeId") == "TMB_NONSYNONYMOUS"
        ]
    ).dropna(subset=["tmb"]).drop_duplicates("sampleId")
    tmb.to_csv(TCGA_RAW_DIR / "tcga_brca_tmb.csv", index=False)
    print(f"Clinical patients: {len(clinical)}")
    print(f"Somatic mutation records: {len(mutations)} across {len(genes)} candidate genes")
    print(f"TMB samples: {len(tmb)}")


if __name__ == "__main__":
    main()
