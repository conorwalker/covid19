# some code re-used from:
# https://github.com/galaxyproject/SARS-CoV-2/blob/master/genomics/4-Variation/fetch_sra_acc.sh
IGNORE_ACCESSIONS="ignore_accessions.txt"
SEEN_ACCESSIONS="seen_accessions.txt"

if [ ! -f "$IGNORE_ACCESSIONS" ]; then
    touch "$IGNORE_ACCESSIONS";
fi

# fetch maintained list of SARS-CoV-2 sequences from SRA
echo "Fetching from genbank..."
curl -s https://www.ncbi.nlm.nih.gov/core/assets/genbank/files/ncov-sequences.yaml --output ncov-sequences.yaml
grep -o "SRR[[:digit:]]\+" ncov-sequences.yaml > genbank.txt

# get acccessions from SRA, likely redundant from previous step
echo "Fetching from SRA..."
python scripts/get_sra_accessions.py

# get accessions from ENA
echo "Fetching from ENA..."
curl -s 'https://www.ebi.ac.uk/ena/browser/api/xml/links/taxon?accession=2697049&result=read_run&download=true' | grep -o "SRR[[:digit:]]\+" > ena.txt

# get concatted list of SARS-CoV-2 tagged sequencing runs
cat ena.txt genbank.txt sra.txt | sort | uniq > current.txt

# clean up
rm ena.txt sra.txt genbank.txt

# get metadata for each new sample
# this step could be sped up with a bulk download, didn't change yet
echo "Fetching metadata..."
for accession in $(<current.txt); do
    if ! grep -q "$accession" "$IGNORE_ACCESSIONS"; then
        # check if already processed
        if ! grep -q "$accession" current_metadata.txt; then
            # only get header if no metadata file exists
            if [ -f current_metadata.txt ]; then
                output=$(esearch -db sra -q "$accession" | efetch -format runinfo | tail -n 2);
            else
                output=$(esearch -db sra -q "$accession" | efetch -format runinfo);
            fi
            # adding accessions with no valid metadata to a list of accessions to be ignored
            # in future runs
            if [ -z "$output" ]; then
                echo "$accession" >> "$IGNORE_ACCESSIONS";
            fi
            echo "$output" >> current_metadata.txt;
        fi
    fi
done

# clean up
#rm current.txt

# remove blank lines, some SRA accessions don't link to valid metadata
sed -i '/^$/d' current_metadata.txt

# generate new files containing most recent illumina and ONT sequencing run accession IDs
grep -e 'Illumina\|NextSeq' current_metadata.txt | cut -d"," -f1 > current_illumina.txt
grep -e 'MinION\|GridION' current_metadata.txt | cut -d"," -f1 > current_ONT.txt

# generate new list of previously processed accession IDs
cat current_illumina.txt current_ONT.txt > "$SEEN_ACCESSIONS"
