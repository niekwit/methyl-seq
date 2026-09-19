import os

# Checked-in static resource used directly (no download rule needed) by the
# "mm39" branch below -- workflow/resources/ is always reachable via this
# module's own file location regardless of the analysis --directory, unlike
# a bare "workflow/resources/..." relative path which would only resolve by
# accident (see the same concern documented for WORKFLOW_SCRIPTS/
# WORKFLOW_RESOURCES in resources.smk).
_ICR_REGIONS_MM39 = os.path.normpath(
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "resources",
        "icr_regions_mm39.bed",
    )
)

# Same idea for hg38 -- the 25 canonical human ICRs (Skaar et al. 2012, ILAR
# J, PMID 23744971), with hg38 coordinates and gene-symbol annotation from
# humanicr.org (Sanchez-Delgado et al. 2022, Epigenetics, PMID 35786392) via
# UCSC's REST API. Regenerate with workflow/scripts/get_icr_regions_hg38.py;
# see workflow/resources/icr_regions_hg38.references.tsv for full per-region
# citations.
_ICR_REGIONS_HG38 = os.path.normpath(
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "resources",
        "icr_regions_hg38.bed",
    )
)


class Resources:
    """Gets URLs and file names of fasta and GTF files for a given genome and build"""

    # Create genome directory
    os.makedirs("resources/downloaded_fasta", exist_ok=True)

    def __init__(self, genome, build):
        self.genome = genome
        self.build = build

        # Paternally-imprinted-region (ICR) heatmap: mm39 (4 hand-curated
        # ICRs, see workflow/resources/icr_regions_mm39.bed) and hg38 (25
        # canonical ICRs, see workflow/resources/icr_regions_hg38.bed and
        # its .references.tsv) for now, plus a single-region (Rasgrf1 only)
        # subset for the "test" genome in its own shifted mini-genome
        # coordinates, so this code path is exercised in CI. None (feature
        # off) for every other genome; set below.
        self.icr_regions_url = None
        self.icr_regions = None

        # Base URLs
        base_url_ens = f"https://ftp.ensembl.org/pub/release-{build}/"

        if "hg" in genome:
            if genome == "hg19":
                name = "GRCh37"
                ucsc_build = "hg19"
            elif genome == "hg38":
                name = "GRCh38"
                ucsc_build = "hg38"

            # Create URLs for genome files
            self.fasta_url = f"{base_url_ens}fasta/homo_sapiens/dna/Homo_sapiens.{name}.dna.primary_assembly.fa.gz"
            self.gtf_url = (
                f"{base_url_ens}gtf/homo_sapiens/Homo_sapiens.{name}.{build}.gtf.gz"
            )
            self.regulatory_gtf_url = f"{base_url_ens}regulation/homo_sapiens/{name}/annotation/Homo_sapiens.{name}.regulatory_features.v{build}.gff3.gz"
            self.cpg_islands_url = f"https://hgdownload.soe.ucsc.edu/goldenPath/{ucsc_build}/database/cpgIslandExt.txt.gz"
            self.repeat_mask_url = f"https://hgdownload.soe.ucsc.edu/goldenPath/{ucsc_build}/database/rmsk.txt.gz"
            if genome == "hg38":
                # Static, already-checked-in file -- no download rule needed
                self.icr_regions = _ICR_REGIONS_HG38

        elif "mm" in genome:
            if genome == "mm38":
                name = "GRCm38"
                ucsc_build = "mm10"
            elif genome == "mm39":
                name = "GRCm39"
                ucsc_build = "mm39"

            # Create URLs for genome files
            self.fasta_url = f"{base_url_ens}fasta/mus_musculus/dna/Mus_musculus.{name}.dna.primary_assembly.fa.gz"
            self.gtf_url = (
                f"{base_url_ens}gtf/mus_musculus/Mus_musculus.{name}.{build}.gtf.gz"
            )
            self.regulatory_gtf_url = f"{base_url_ens}regulation/mus_musculus/{name}/annotation/Mus_musculus.{name}.regulatory_features.v{build}.gff3.gz"
            self.cpg_islands_url = f"https://hgdownload.soe.ucsc.edu/goldenPath/{ucsc_build}/database/cpgIslandExt.txt.gz"
            self.repeat_mask_url = f"https://hgdownload.soe.ucsc.edu/goldenPath/{ucsc_build}/database/rmsk.txt.gz"
            if genome == "mm39":
                # Static, already-checked-in file -- no download rule needed
                self.icr_regions = _ICR_REGIONS_MM39

        elif "dm" in genome:
            if genome == "dm6":
                name = "BDGP6.46"

            self.fasta_url = f"{base_url_ens}fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.{name}.dna.toplevel.fa.gz"
            self.gtf_url = f"{base_url_ens}gtf/drosophila_melanogaster/Drosophila_melanogaster.{name}.{build}.gtf.gz"
            # Ensembl does not provide regulatory build data for Drosophila melanogaster
            self.regulatory_gtf_url = None
            # UCSC does not provide a CpG Islands track for Drosophila melanogaster
            self.cpg_islands_url = None
            # UCSC does not provide a RepeatMasker track for Drosophila melanogaster
            self.repeat_mask_url = None

        elif "test" in genome:
            # The test reads (.test/reads/) are real Bismark-aligned reads
            # subsetted to a small mouse (mm39) locus around Rasgrf1 -- see
            # .test/make_test_data.py. They are NOT human, so this genome
            # must be built from the matching mini reference carved out
            # alongside them (committed to this repo), not downloaded from
            # Ensembl/UCSC.
            # Pinned to a commit SHA, not "main" -- a floating branch ref
            # means every CI run refetches whatever "main" happens to be at
            # request time, which can race against this same session's own
            # pushes (or a stale raw.githubusercontent.com CDN response) and
            # is not reproducible. Update this SHA when .test/resources/
            # content actually changes.
            test_base_url = "https://github.com/niekwit/methyl-seq/raw/1d9706dd59004b62ecb969c9c55198afdbf391e6/.test/resources"
            self.fasta_url = f"{test_base_url}/genome.fa.gz"
            self.gtf_url = f"{test_base_url}/genes.gtf.gz"
            # Subsetted regulatory build / CpG island / RepeatMasker data
            # all exist (carved from the real mm39 tracks by
            # make_test_data.py) so the promoter/cpg_islands/repeat_mask
            # code paths are all exercised in CI
            self.regulatory_gtf_url = f"{test_base_url}/regulatory_features.gff3.gz"
            self.cpg_islands_url = f"{test_base_url}/cpgIslandExt.txt.gz"
            self.repeat_mask_url = f"{test_base_url}/rmsk.txt.gz"
            # Single-region (Rasgrf1 only) subset, in the mini genome's own
            # shifted coordinates -- see .test/resources/icr_regions.bed
            self.icr_regions_url = f"{test_base_url}/icr_regions.bed"
            self.icr_regions = "resources/icr_regions.bed"

        else:
            raise ValueError("Genome {genome} not supported")

        # Downloaded unzipped file names
        self.fasta = self._file_from_url(self.fasta_url)
        self.filtered_fasta = self._filtered_fasta_from_url(self.fasta_url)
        self.gtf = self._file_from_url(self.gtf_url)
        self.regulatory_gtf = (
            self._file_from_url(self.regulatory_gtf_url)
            if self.regulatory_gtf_url
            else None
        )

        # Control fasta from NEB GitHub
        # https://github.com/FelixKrueger/Bismark/issues/166#issuecomment-378349782
        self.control_fasta_url = "https://raw.githubusercontent.com/nebiolabs/EM-seq/refs/heads/master/assets/methylation_controls.fa"
        self.control_fasta = self._file_from_url(self.control_fasta_url)

        # Prepare CpG island BED file
        # Equivalent to downloading cpg_islands.bed.gz from
        # https://genome.ucsc.edu/cgi-bin/hgTables
        #   Group: Expression and Regulation
        #   Track: CpG Islands
        #   Table: cpgIslandExt
        #   Region: genome
        #   Output format: BED
        #   File: cpg_islands.bed.gz
        # but fetched non-interactively from UCSC's goldenPath database dump
        # (same underlying cpgIslandExt table) and converted to BED by the
        # prepare_cpg_islands rule, since hgTables itself has no stable
        # scriptable download URL.
        self.cpg_islands = "resources/cpg_islands.bed" if self.cpg_islands_url else None

        # Prepare RepeatMasker BED file (transposable elements only)
        # Equivalent to downloading repeat_mask.txt.gz from
        # https://genome.ucsc.edu/cgi-bin/hgTables
        #   Group: Variation and Repeats
        #   Track: RepeatMasker
        #   Table: rmsk
        #   Region: genome
        #   Output format: All fields from selected table
        #   File: repeat_mask.txt.gz
        # but fetched non-interactively from UCSC's goldenPath database dump
        # (same underlying rmsk table), converted to BED, and filtered to
        # drop non-transposable-element repeat classes/families (see
        # workflow/resources/nonTE_repClasses.txt) by the prepare_repeat_mask
        # rule, since hgTables itself has no stable scriptable download URL.
        self.repeat_mask = "resources/repeat_mask.bed" if self.repeat_mask_url else None

    def _file_from_url(self, url):
        """Returns file path for unzipped downloaded file"""
        return f"resources/downloaded_fasta/{os.path.basename(url).replace('.gz','')}"

    def _filtered_fasta_from_url(self, url):
        """Returns file path for unzipped downloaded file"""
        return f"resources/downloaded_fasta/{os.path.basename(url).replace('.fa.gz','_filtered.fa')}"
