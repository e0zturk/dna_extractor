import os
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio import BiopythonWarning
import mygene
import requests
import json
import argparse
import shutil
import time
import warnings

arg = argparse.ArgumentParser(description="Get TSS information for a given protein sequence")
arg.add_argument("-i", "--input", required=False, help="Input FASTA file")
arg.add_argument("-r", "--raw", required=False, help="Raw json and xml files otherwise default to the current directory")
arg.add_argument("-mail", "--email", required=False, help="Email address to use for NCBI Entrez API (optional, but recommended to avoid rate limits)")
args = vars(arg.parse_args())

if args["email"]:
    Entrez.email = args["email"]
else:
    Entrez.email = "eo11003@mail2.gantep.edu.tr" # Default email address, replace with your own if desired


fasta_file = args["input"]
output_dir = args["raw"]

if not fasta_file or not output_dir:
    print("Error: Missing required arguments. Both --input (-i) and --raw (-r) must be provided.")
    exit(1)

    
record = SeqIO.read(fasta_file, format="fasta")


def WW_handle(seq, output_dir):

    if len(seq) < 25:
        print("Sequence is too short for BLAST search. Please provide a sequence of at least 25 amino acids.")
        return None
    invalid_chars = set("BJOUXZ")
    if any(char in invalid_chars for char in seq):
        print("Sequence contains invalid characters (B, J, O, U, X, Z). Please provide a valid protein sequence.")
        return None
    
    # BiopythonWarning'leri gizle (taking longer uyarısı)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        result_handle = NCBIWWW.qblast("blastp", "refseq_protein", seq, hitlist_size=5)

    with open(os.path.join(output_dir, f"{record.id}_blast_result.xml"), "w") as out_handle:
        out_handle.write(result_handle.read())

    result_handle.close()
    
    output_file = os.path.join(output_dir, f"{record.id}_blast_result.xml")

    return output_file


def Xml_parser(output_dir):
    result_file = os.path.join(output_dir, f"{record.id}_blast_result.xml")
    
    with open(result_file) as result:
        blast_record = NCBIXML.read(result)
        
        if not blast_record.alignments:
            print("ERROR: BLAST search returned no alignments!")
            print("Possible reasons:\n")
            print("  1. Protein may not be from human (Homo sapiens)")
            print("  2. Protein not found in NCBI RefSeq database")
            print("  3. Sequence may contain errors or be invalid\n\n")
            raise ValueError("No BLAST matches found. Verify that your protein sequence is from a human gene or change the sequence")
        
        alignment = blast_record.alignments[0]
        hsp = alignment.hsps[0]
        accession = alignment.hit_id

    return accession


def get_nt(accession):
    max_retries = 3
    for attempt in range(max_retries):
        try:
            handle = Entrez.esearch(db="protein", term=accession)
            result = Entrez.read(handle)
            handle.close()
            
            if not result["IdList"]:
                print(f"WARNING: No protein UID found for accession {accession}")
                return None
            protein_uid = result["IdList"][0]
            
            handle = Entrez.elink(dbfrom="protein", db="nucleotide", id=protein_uid)
            result = Entrez.read(handle)
            handle.close()

            if not result[0].get("LinkSetDb"):
                print(f"WARNING: No nucleotide links found")
                return None
            
            nt_ids = [link["Id"] for link in result[0]["LinkSetDb"][0]["Link"]]
            
            for nt_id in nt_ids:
                handle = Entrez.efetch(db="nucleotide", id=nt_id, rettype="gb", retmode="text")
                gb_text = handle.read()
                handle.close()
                
                for line in gb_text.split("\n"):
                    if "NM_" in line:
                        for word in line.split():
                            if word.startswith("NM_"):
                                return word
            return None
            
        except RuntimeError as e:
            print(f"WARNING: Attempt {attempt+1}/{max_retries} failed: {str(e)[:60]}...")
            if attempt < max_retries - 1:
                print("INFO: Retrying in 5 seconds...")
                time.sleep(5)
            else:
                print("ERROR: All retry attempts failed. NCBI server may be temporarily unavailable.")
                return None


def get_ensembl(ensgene, result_dir):

    server = "https://rest.ensembl.org"
    ext = f"/lookup/id/{ensgene}?expand=1"
    
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        return None
    
    decoded = r.json()

    with open(os.path.join(result_dir, f"{ensgene}_ensembl.json"), "w") as js:
        json.dump(decoded, js, indent=4, ensure_ascii=False)

    json_file = os.path.join(result_dir, f"{ensgene}_ensembl.json")

    return json_file

basepath = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":

    xml_file = os.path.join(output_dir, f"{record.id}_blast_result.xml")
    if not os.path.exists(xml_file):
        WW_handle(record.seq, output_dir)

    pr_accession = Xml_parser(output_dir).split("|")[1]

    nm_accession = get_nt(pr_accession)
    
    if not nm_accession:
        print("ERROR: Could not retrieve nucleotide accession from NCBI Entrez.")
        print("NCBI servers may be temporarily unavailable. Please try again later.")
        print("https://www.ncbi.nlm.nih.gov/home/about/status/")
        import sys
        sys.exit(1)

    mg = mygene.MyGeneInfo()

    gene_info = mg.query(scopes = "refseq", q = str(nm_accession.split(".")[0]),
                        fields="ensembl.gene", species = "human" )

    try:
        ensg = gene_info["hits"][0]["ensembl"]["gene"]
    except (IndexError, KeyError):
        ensg = None

    # JSON kontrolü
    if not os.path.exists(os.path.join(output_dir, f"{ensg}_ensembl.json")):
        ensembl_json = get_ensembl(ensg, output_dir)
    else:
        with open(os.path.join(output_dir, f"{ensg}_ensembl.json"), "r") as js:
            ensembl_json = json.load(js)
    if ensembl_json:
        jsoncache_dir = os.path.join(basepath, "jsoncache")
        os.makedirs(jsoncache_dir, exist_ok=True)
        
        json_filename = os.path.basename(ensembl_json)
        jsoncache_path = os.path.join(jsoncache_dir, json_filename)
        shutil.copy(ensembl_json, jsoncache_path)

    print(f"INFO : JSON also cached to: {jsoncache_path}")
    print("Running BLAST search and parsing results...")
    print(f"Protein Accession: {pr_accession}")
    print("-"*50)
    print(f"NM Accession: {nm_accession}")
    print("-"*50)
    print(f"Ensembl Gene ID: {ensg}")
    print("-"*50)
    print(f"XML file saved at: {os.path.join(output_dir, f'{record.id}_blast_result.xml')}")
    print("_"*10)
    print(f"Ensembl JSON file saved at: {os.path.join(output_dir, f'{ensg}_ensembl.json')}")
    print("="*50)