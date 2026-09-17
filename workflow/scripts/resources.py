import os


class Resources:
    """Gets URLs and file names of fasta and GTF files for a given genome and build"""

    # Create genome directory
    os.makedirs("resources/downloaded_fasta", exist_ok=True)

    def __init__(self, genome, build):
        self.genome = genome
        self.build = build

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
            test_base_url = (
                "https://github.com/niekwit/methyl-seq/raw/main/.test/resources"
            )
            self.fasta_url = f"{test_base_url}/genome.fa.gz"
            self.gtf_url = f"{test_base_url}/genes.gtf.gz"
            # Subsetted regulatory build / CpG island / RepeatMasker data
            # all exist (carved from the real mm39 tracks by
            # make_test_data.py) so the promoter/cpg_islands/repeat_mask
            # code paths are all exercised in CI
            self.regulatory_gtf_url = f"{test_base_url}/regulatory_features.gff3.gz"
            self.cpg_islands_url = f"{test_base_url}/cpgIslandExt.txt.gz"
            self.repeat_mask_url = f"{test_base_url}/rmsk.txt.gz"

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
