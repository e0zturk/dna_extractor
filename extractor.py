import os
import subprocess
import argparse
import json
import glob
import shutil
import gzip
import tarfile
import platform
import sys
try:
    from pyfaidx import Fasta
except ImportError:
    raise ImportError("pyfaidx is not installed. Please check in requirements.txt"
                      "recommended run 'pip install -r requirements.txt'")


SYSTEM = platform.system()
IS_WINDOWS = SYSTEM == "Windows"
IS_LINUX = SYSTEM == "Linux"
IS_MACOS = SYSTEM == "Darwin"

def OS_check():
    if IS_WINDOWS:
        return "python"
    else:
        return "python3"
def OS_terminal_download():
    if IS_WINDOWS:
        return ["powershell", "-Command"]
    else:
        return ["curl", "-L"]

ModeUsages = " Mode options:\n\n" \
             "  canonical: Only the canonical transcript sequence is included in the output.\n\n" \
             "  no-repeats: Repetitive elements (e.g., Alu, LINE-1) are masked from the extended TSS sequence.\n\n" \
             "  no-exons: Exon sequences from all transcripts are masked from the extended TSS sequence.\n\n" \
             "  all: The full extended TSS sequence is included without any filtering (default)\n\n"

arg = argparse.ArgumentParser(description="Get extended TSS sequences from a FASTA file",formatter_class=argparse.RawDescriptionHelpFormatter, epilog=ModeUsages)
arg.add_argument("-i", "--input", required=False, help="Input FASTA file")
arg.add_argument("-up", "--upstream", type=int, required=False, help="Number of base pairs upstream of the transcription start site to include in the extended TSS sequence")
arg.add_argument("-down", "--downstream", type=int, required=False, help="Number of base pairs downstream of the transcription start site to include in the extended TSS sequence")
arg.add_argument("-o", "--output", required=False, help="Output file's name, default is <gene_id>_extended-TSS.fasta")
arg.add_argument("-mail", "--email", required=False, help="Email address to use for NCBI Entrez API (optional, but recommended to avoid rate limits)")

arg.add_argument("-mod", "--mode", required=False, help="Choose processing mode")
args = vars(arg.parse_args())

basepath = os.path.dirname(os.path.abspath(__file__))
raw_dir = os.path.join(basepath, "raw_dir")

print("INFO : Checking for JSON files in the raw directory.")

if not args["input"]:
    json_data = glob.glob(os.path.join(raw_dir, f"*.json"))
    if not json_data:
        raise OSError("No JSON files found in the raw directory")
else:
    if os.path.exists(raw_dir):
        shutil.rmtree(raw_dir)
    os.makedirs(raw_dir, exist_ok=True)
    print("INFO : Running 'get_tss.py' to get the transcription start site information.\n   Creating a JSON file could take 3-5 minutes\n\n" \
          "     <|If long runtime is observed, please try to run with your personal email address using the -mail option|>\n(Ctrl+C to cancel)\n\n" \
          "     <|To more info about visit and sign up for NCBI Entrez API: https://www.ncbi.nlm.nih.gov/account/|>\n\n")
    tss_cmd = [OS_check(), os.path.join(basepath, "get_tss.py"), "-i", args["input"], "-r", raw_dir]
    if args["email"]:
        tss_cmd.extend(["-mail", args["email"]])
    subprocess.run(tss_cmd)
    json_data = glob.glob(os.path.join(raw_dir, f"*.json"))
    if not json_data:
        raise OSError("Even if run get_tss.py, no JSON files found in the raw directory")

print("INFO : JSON file is ready, now checking for the genome file. If not found, it will be downloaded from Ensembl")

if not os.path.exists(os.path.join(basepath, "genome", "Homo_sapiens.GRCh38.fa")):
    os.makedirs(os.path.join(basepath, "genome"), exist_ok=True)
    gz_file = os.path.join(basepath, "genome", "Homo_sapiens.GRCh38.fa.gz")
    fa_file = os.path.join(basepath, "genome", "Homo_sapiens.GRCh38.fa")
    URL = "https://ftp.ensembl.org/pub/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
    OS_terminal_run = OS_terminal_download()
    if IS_WINDOWS:
        cmd_download = [f"$ProgressPreference = 'SilentlyContinue'; Invoke-WebRequest -Uri {URL} -OutFile {gz_file}"]
    else:
        cmd_download = [URL, "-o", gz_file]
    print("INFO : Genome (Homo_sapiens.GRCh38.fa.gz) downloading and extracting. It takes 2-3 minutes, please wait...")
    subprocess.run(OS_terminal_run + cmd_download)
    with gzip.open(gz_file, 'rb') as f_in:
        with open(fa_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
else:
    pass

print("INFO : Genome is opening.")

with open(json_data[0],"r") as js:
    gene_location = json.load(js)

if gene_location["species"] != "homo_sapiens":
    print("ERROR: Out of scope species was detected")
    raise ValueError("Only Homo sapiens is supported")

chrom = gene_location['seq_region_name']
strand = gene_location['strand']


if not args["downstream"] or not args["upstream"]:
    print("INFO : Downstream and Upstream values not provided or missing, using default value of 500(downstream), 1000(upstream)")

upstream = int(args["upstream"]) if args["upstream"] else 1000
downstream = int(args["downstream"]) if args["downstream"] else 500

if os.path.exists(os.path.join(basepath, "genome", "Homo_sapiens.GRCh38.fa")):
    genome = Fasta(os.path.join(basepath, "genome", "Homo_sapiens.GRCh38.fa"))
else:
    raise OSError("Genome file not found")

if strand == -1:
    transcript_starts = [int(t["start"]) for t in gene_location["Transcript"]]

    TSS_start = min(transcript_starts)

    transcript_ends = [int(t["end"]) for t in gene_location["Transcript"]]

    TSS_end = max(transcript_ends)

    prom_start = TSS_start - upstream
    prom_end = TSS_end + downstream

    target = genome[chrom][prom_start-1:prom_end]

else:

    transcript_starts = [int(t["start"]) for t in gene_location["Transcript"]]
    TSS_start = max(transcript_starts)

    transcript_ends = [int(t["end"]) for t in gene_location["Transcript"]]
    TSS_end = min(transcript_ends)

    prom_start = TSS_start + upstream
    
    prom_end = TSS_end - downstream

    target = genome[chrom][prom_start-1:prom_end]


fasta_seq = str(target) if strand == 1 else str(target.reverse.complement)

if not args["output"]:
    args["output"] = f"{gene_location['id']}_extended-TSS.fasta"

"""
mask_exons used to mask exon sequences from the extended TSS sequence.
a gunzip file made as an output file includes overall transcript 'masked-exon' fasta sequence
"""
def mask_exons(js_data, prom_start, prom_end, full_seq, sstrand):
    transcripts = {i: {"exons": js_data["Transcript"][i]["Exon"]} for i in range(len(js_data["Transcript"]))}
    
    Exons = {
        f"transcript_{i+1}": {
            f"exon_{j+1}": {'start': exon['start'], 'end': exon['end']}
            for j, exon in enumerate(transcripts[i]['exons'])
        }
        for i in range(len(transcripts))
    }
    
    
    raw_targets = {}
    for i in range(len(Exons)):
        t_key = f"target{i+1}"
        raw_targets[t_key] = {}
        transcript_key = f"transcript_{i+1}"
        t_exons = Exons[transcript_key]
        for j in range(len(t_exons)):
            exon_out_key = f"exon{j+1}"
            exon_coords_key = f"exon_{j+1}"
            coords = t_exons[exon_coords_key]
            
            raw_targets[t_key][exon_out_key] = coords 
    
    masked_seqs = {}
    for t_key, exons in raw_targets.items():
        seq_list = list(full_seq.upper())
        
        for exon_key, coords in exons.items():
            if sstrand == 1:
                idx_start = coords['start'] - 1 - prom_start
                idx_end = coords['end'] - prom_start
            else:
                idx_start = prom_end - coords['end']
                idx_end = prom_end - coords['start'] + 1
            
            idx_start = max(0, idx_start)
            idx_end = min(len(seq_list), idx_end)
            
            if idx_start < idx_end:
                seq_list[idx_start:idx_end] = ["N"] * (idx_end - idx_start)
        masked_seqs[t_key] = "".join(seq_list)

    return masked_seqs


masked_exon_seqs = mask_exons(gene_location, prom_start, prom_end, fasta_seq, strand)

if not args["output"].endswith(".fasta"):
    args["output"] += ".fasta"

def write_mskexon(masked_seqs_dict, js_dict, chrom, strand, prom_start, prom_end):
    sorted_keys = sorted(
        masked_seqs_dict.keys(),
        key=lambda x: int(x.replace("target", ""))
    )
    gene_id = js_dict.get("id")
    written_files = []
    
    for t_key in sorted_keys:
        transcript_index = int(t_key.replace("target", "")) - 1
        transcript_id = js_dict["Transcript"][transcript_index]["id"]
        chrom_local = js_dict.get("seq_region_name", chrom)
        strand_local = js_dict.get("strand", strand)
        seq = masked_seqs_dict[t_key]
        
        if not seq or len(seq) == 0:
            continue
        
        valid_chars = set("ATGCatgcNn")
        invalid_chars = set(seq) - valid_chars
        if invalid_chars:
            print(f"WARNING : Transcript {transcript_id} contains invalid characters: {invalid_chars}")
        
        out_path = f"masked-exons_{transcript_id}.fasta"
        
        try:
            with open(out_path, "w", encoding="utf-8") as out:
                out.write(f">{transcript_id}|{chrom_local}:{prom_start}-{prom_end}|strand:{strand_local}|masked_exons\n")
                for i in range(0, len(seq), 80):
                    out.write(seq[i:i+80] + "\n")
            written_files.append(out_path)
        except IOError as e:
            raise IOError(f"Error writing file {out_path}: {e}")
    
    tar_gz_path = f"{gene_id}_masked-exons.tar.gz"
    
    with tarfile.open(tar_gz_path, "w:gz") as tar:
        for fasta_file in written_files:
            tar.add(fasta_file, arcname=os.path.basename(fasta_file))
    
    for fasta_file in written_files:
        os.remove(fasta_file)
    
    print(f"INFO : Exon-masked extended TSS sequences saved to {tar_gz_path}")


def write_norepeats(full_seq,chrom, strand, prom_start, prom_end):
    cleaned_seq = ''.join([c for c in full_seq if not c.islower()])
    out_path = f"{gene_location['id']}_norepeats.fasta"
    with open(out_path, "w") as out:
        out.write(f">{gene_location['id']}|{chrom}:{prom_start}-{prom_end}|strand:{strand}|norepeats\n")
        for i in range(0, len(cleaned_seq), 80):
            out.write(cleaned_seq[i:i+80] + "\n")
    print(f"INFO : No-repeat extended TSS sequence saved to {out_path}")


def write_canonical(json_data, genome, chromomosome, strand, upstream, downstream):
    canonical_seq = None
    for transcript in json_data["Transcript"]:
        if transcript.get("is_canonical", False):
            canonical_seq = ""            
            pr_start = transcript['start'] - 1
            pr_end = transcript['end']
            if strand == -1:
                seq_start = pr_start - upstream
                seq_end = pr_end + downstream
            else:
                seq_start = pr_start + upstream
                seq_end = pr_end - downstream
              
            exon_seq_obj = genome[chromomosome][seq_start:seq_end]
            exon_seq = exon_seq_obj.seq if strand == 1 else exon_seq_obj.reverse.complement.seq
            canonical_seq += exon_seq.upper()
            break

    if canonical_seq:
        out_path = f"{gene_location['id'].replace('.fasta', '')}_canonical.fasta"
        with open(out_path, "w") as out:
            out.write(f">{json_data['id']}|{chrom}:{prom_start}-{prom_end}|strand:{strand}|canonical\n")
            for i in range(0, len(canonical_seq), 80):
                out.write(canonical_seq[i:i+80] + "\n")
        print(f"INFO : Canonical extended TSS sequence saved to {out_path}")
    else:
        print("INFO : No canonical transcript found, skipping canonical output.")

def write_all(full_seq,chrom, strand, prom_start, prom_end):
    out_path = f"{gene_location['id']}_all.fasta"
    with open(out_path, "w") as out:
        out.write(f">{gene_location['id']}|{chrom}:{prom_start}-{prom_end}|strand:{strand}\n{full_seq.upper()}\n")
    print(f"INFO : Full extended TSS sequence successfully saved to {out_path}")



if args["mode"] == "no-exons":
    write_mskexon(masked_exon_seqs, gene_location, chrom, strand, prom_start, prom_end)
elif args["mode"] == "no-repeats":
    write_norepeats(fasta_seq, chrom, strand, prom_start, prom_end)
elif args["mode"] == "canonical":
    write_canonical(gene_location, genome, chrom, strand, upstream, downstream)
elif args["mode"] == "all":
    write_all(fasta_seq, chrom, strand, prom_start, prom_end)
else:
    print("BatchError: Invalid mode selected. Please choose from <canonical>, <no-repeats>, <no-exons>, or <all>.")

if not args["mode"]:
    write_all(fasta_seq, chrom, strand, prom_start, prom_end)
    print("INFO : No mode selected, defaulting to 'all' mode. Full extended TSS sequence saved.")



genome.close()
