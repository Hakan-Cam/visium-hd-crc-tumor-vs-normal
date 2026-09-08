#!/usr/bin/env bash
# download_data.sh
#
# Fetches the two Visium HD samples (P5 tumor + P5 normal-adjacent colon,
# from Oliveira et al. 2025, Nature Genetics, 10x Genomics -- GEO series
# GSE280315) and the Pelka et al. 2021 (Cell) human CRC scRNA-seq atlas
# (GEO series GSE178341) used as the RCTD deconvolution reference.
#
# Run this from the repo root:
#   bash download_data.sh
#
# Requires: curl, gunzip. Needs a normal (unrestricted) internet
# connection -- this will NOT work from a network-locked-down sandbox.
#
# Total download: ~530 MB (Visium HD, both samples, matrix+coords+images
# only -- NOT the multi-GB full-resolution whole-slide images, which
# aren't needed for this analysis) + ~1.1 GB (scRNA-seq reference).

set -euo pipefail

DATA_DIR="data"
VHD_DIR="${DATA_DIR}/visium_hd"
REF_DIR="${DATA_DIR}/scref"

mkdir -p "${VHD_DIR}/P5_Tumor/binned_outputs/square_008um/spatial"
mkdir -p "${VHD_DIR}/P5_NormalAdjacent/binned_outputs/square_008um/spatial"
mkdir -p "${REF_DIR}"

# GEO's stable per-file download convention: drop the accession's last 3
# digits, append "nnn", to get the bucket directory.
#   GSM8594569 -> GSM8594nnn   GSM8594571 -> GSM8594nnn   GSE178341 -> GSE178nnn
GEO_SAMPLES="https://ftp.ncbi.nlm.nih.gov/geo/samples"
GEO_SERIES="https://ftp.ncbi.nlm.nih.gov/geo/series"

fetch() {
  # fetch <url> <destination_file>
  local url="$1" dest="$2"
  if [ -s "$dest" ]; then
    echo "  already have $(basename "$dest"), skipping"
    return
  fi
  echo "  downloading $(basename "$dest") ..."
  curl -sSL --fail -o "$dest" "$url"
}

# ---------------------------------------------------------------------------
echo "== P5 Tumor (GSM8594569) =="
S="${GEO_SAMPLES}/GSM8594nnn/GSM8594569/suppl"
D="${VHD_DIR}/P5_Tumor/binned_outputs/square_008um"
fetch "${S}/GSM8594569_P5CRC_filtered_feature_bc_matrix.h5" \
      "${D}/filtered_feature_bc_matrix.h5"
fetch "${S}/GSM8594569_P5CRC_tissue_positions.parquet.gz" \
      "${D}/spatial/tissue_positions.parquet.gz"
fetch "${S}/GSM8594569_P5CRC_scalefactors_json.json.gz" \
      "${D}/spatial/scalefactors_json.json.gz"
fetch "${S}/GSM8594569_P5CRC_tissue_lowres_image.png.gz" \
      "${D}/spatial/tissue_lowres_image.png.gz"
gunzip -f "${D}/spatial/"*.gz

echo "== P5 Normal-Adjacent (GSM8594571) =="
S="${GEO_SAMPLES}/GSM8594nnn/GSM8594571/suppl"
D="${VHD_DIR}/P5_NormalAdjacent/binned_outputs/square_008um"
fetch "${S}/GSM8594571_P5NAT_filtered_feature_bc_matrix.h5" \
      "${D}/filtered_feature_bc_matrix.h5"
fetch "${S}/GSM8594571_P5NAT_tissue_positions.parquet.gz" \
      "${D}/spatial/tissue_positions.parquet.gz"
fetch "${S}/GSM8594571_P5NAT_scalefactors_json.json.gz" \
      "${D}/spatial/scalefactors_json.json.gz"
fetch "${S}/GSM8594571_P5NAT_tissue_lowres_image.png.gz" \
      "${D}/spatial/tissue_lowres_image.png.gz"
gunzip -f "${D}/spatial/"*.gz

echo "== Pelka et al. 2021 CRC scRNA-seq reference (GSE178341) =="
S="${GEO_SERIES}/GSE178nnn/GSE178341/suppl"
fetch "${S}/GSE178341_crc10x_full_c295v4_submit.h5" \
      "${REF_DIR}/crc10x_full_c295v4_submit.h5"
fetch "${S}/GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz" \
      "${REF_DIR}/crc10x_full_c295v4_submit_cluster.csv.gz"
fetch "${S}/GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz" \
      "${REF_DIR}/crc10x_full_c295v4_submit_metatables.csv.gz"
gunzip -fk "${REF_DIR}/"*.csv.gz   # -k: keep the .gz too, harmless either way

echo ""
echo "Done. Verify with:"
echo "  find ${DATA_DIR} -type f | sort"
echo ""
echo "NOTE on bin size: GEO depositors typically upload one representative"
echo "Space Ranger bin resolution per sample rather than all three (2/8/16um)."
echo "01_reference_prep.Rmd / 02_visium_hd_analysis.Rmd both include a quick"
echo "sanity check (spot spacing from scalefactors_json.json) to confirm what"
echo "you actually got lines up with the 8um bins this analysis assumes --"
echo "read that check's output before trusting downstream results."
