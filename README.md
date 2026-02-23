# DNA Extractor

## Overview

DNA Extractor is a bioinformatics tool designed to extract extended transcription start site (TSS) sequences from the human genome. While originally developed for researchers working with G4-quadruplex structural analysis, this tool provides flexible sequence extraction capabilities with multiple filtering modes that make it useful for various genomic research applications.

The tool takes a protein sequence as input, identifies the corresponding gene through BLAST search and NCBI Entrez queries, and extracts customizable genomic regions surrounding the transcription start site. Output sequences are generated in FASTA format with comprehensive metadata.

### Key Features

- Supports human genome sequences only (Homo sapiens, GRCh38)
- Automatic genome file management with Ensembl integration
- Multiple filtering modes for sequence extraction
- Customizable upstream and downstream region sizes
- Generates reference datasets for downstream analysis
- Handles multiple transcript variants

### Important Note

This tool is designed exclusively for the human genome (Homo sapiens). Non-human species input will be rejected with an error message.

---

## Quick Start: Complete Setup from Scratch

### Step 1: Download and Install Miniconda

**Windows PowerShell:**
```powershell
# Download Miniconda installer
Invoke-WebRequest -Uri "https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe" -OutFile "Miniconda3-installer.exe"

# Run installer
.\Miniconda3-installer.exe
```
During installation, make sure to check "Add Miniconda to PATH". After installation, restart PowerShell.

**Linux and macOS Terminal:**
```bash
# Download Miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Run installer
bash Miniconda3-latest-Linux-x86_64.sh
```
Follow the prompts and accept the license agreement. After installation, restart your terminal.

If `conda --version` returns an error, manually add Miniconda to PATH:
- **Windows:** Add `C:\Users\YourUsername\miniconda3\Scripts` to System Environment Variables PATH
- **Linux/macOS:** Run `echo 'export PATH="$HOME/miniconda3/bin:$PATH"' >> ~/.bashrc && source ~/.bashrc`

### Step 2: Download DNA Extractor Repository

**Option A: Using Git (if Git is installed)**
```bash
# Windows PowerShell, Linux, and macOS
git clone https://github.com/yourusername/dna_extractor.git
cd dna_extractor
```

**Option B: Download ZIP (no Git required)**
1. Go to: https://github.com/yourusername/dna_extractor
2. Click **Code** → **Download ZIP**
3. Extract the ZIP file to your desired location
4. Open terminal/PowerShell and navigate to folder:
```bash
# Windows PowerShell
cd path\to\dna_extractor

# Linux and macOS
cd path/to/dna_extractor
```

### Step 3: Create Conda Environment

```bash
# Windows PowerShell
conda env create -f environment.yml
conda activate dna_extractor

# Linux and macOS
conda env create -f environment.yml
conda activate dna_extractor
```

### Step 4: Verify Installation

```bash
# Windows
python extractor.py -h

# Linux and macOS
python3 extractor.py -h
```

If you see the help message with argument descriptions, the installation was successful!

---

## Installation

### Prerequisites

- Miniconda (downloaded and installed via Quick Start above)
- Python 3.12 or higher (included with Miniconda)
- Internet connection (for genome file download and NCBI queries)

### Important Note on Input File (-i argument)

When running the script for the first time with a new input protein sequence, use the `-i` argument to provide the FASTA file. This will trigger the BLAST search and Ensembl data retrieval process. Once the JSON file is generated and stored in the `raw_dir` directory, subsequent runs do not require the `-i` argument. Using `-i` again will regenerate the data, which is unnecessary and time-consuming. Only use `-i` when analyzing a different protein sequence.

---

## Usage

### Basic Command Structure

```bash
# Windows
python extractor.py -i <input_fasta> -o <output_name> -up <upstream_bp> -down <downstream_bp> -mod <mode> -mail <email>

# Linux and macOS
python3 extractor.py -i <input_fasta> -o <output_name> -up <upstream_bp> -down <downstream_bp> -mod <mode> -mail <email>
```

### Arguments

- `-i, --input`: Input protein sequence in FASTA format (required for first run)
- `-o, --output`: Output filename (default: `<gene_id>_extended-TSS.fasta`)
- `-up, --upstream`: Base pairs upstream of TSS (default: 1000)
- `-down, --downstream`: Base pairs downstream of TSS (default: 500)
- `-mod, --mode`: Filtering mode - see modes section below (default: all)
- `-mail, --email`: NCBI email address (optional but recommended to avoid rate limiting)

### Extraction Modes

**all**: Returns the complete extended TSS sequence without filtering. This is the default mode.

**canonical**: Extracts only the canonical transcript sequence, excluding variant transcripts.

**no-repeats**: Removes repetitive elements (e.g., Alu, LINE-1) from the sequence. These elements are typically represented in lowercase in the reference genome.

**no-exons**: Removes exon sequences from all transcripts in the region. A separate FASTA file is generated for each transcript to ensure precise removal according to each transcript's coordinates. Output is provided as a tar.gz archive containing individual files for easy per-transcript examination.

### Example Commands

Extract extended TSS sequence with default parameters:
```bash
# Windows
python extractor.py -i sample.fasta -o my_gene_sequences.fasta

# Linux and macOS
python3 extractor.py -i sample.fasta -o my_gene_sequences.fasta
```

Extract with custom upstream/downstream sizes:
```bash
# Windows
python extractor.py -i sample.fasta -o my_gene_sequences.fasta -up 2000 -down 1000

# Linux and macOS
python3 extractor.py -i sample.fasta -o my_gene_sequences.fasta -up 2000 -down 1000
```

Extract canonical transcript only:
```bash
# Windows
python extractor.py -i sample.fasta -mod canonical -mail your.email@example.com

# Linux and macOS
python3 extractor.py -i sample.fasta -mod canonical -mail your.email@example.com
```

Remove repetitive elements:
```bash
# Windows
python extractor.py -i sample.fasta -mod no-repeats

# Linux and macOS
python3 extractor.py -i sample.fasta -mod no-repeats
```

Extract intronic sequences (remove exons):
```bash
# Windows
python extractor.py -i sample.fasta -mod no-exons -mail your.email@example.com

# Linux and macOS
python3 extractor.py -i sample.fasta -mod no-exons -mail your.email@example.com
```

---

## Output Files

The tool generates output files in FASTA format with descriptive headers containing metadata about the extracted sequence:

- **all mode**: `<gene_id>_all.fasta` - Complete extended TSS sequence
- **canonical mode**: `<gene_id>_canonical.fasta` - Canonical transcript only
- **no-repeats mode**: `<gene_id>_norepeats.fasta` - Sequence with repetitive elements removed
- **no-exons mode**: `<gene_id>_removed-exons.tar.gz` - Compressed archive containing individual FASTA files for each transcript with exons removed

FASTA headers include the following information:
```
>gene_id|chromosome:start-end|strand:+/-|mode_applied
```

---

## Requirements

The following packages are automatically installed via conda:

- pyfaidx >= 0.7.2
- biopython >= 1.81
- mygene >= 3.2.2
- requests >= 2.31.0

---

## Troubleshooting

**Error: "No JSON files found in the raw directory"**
- Ensure your input FASTA file is valid and contains a proper protein sequence
- Check internet connection for NCBI and Ensembl access
- Verify BLAST query meets minimum length requirements (at least 25 amino acids)

**Error: "Only Homo sapiens is supported"**
- This tool is designed exclusively for human genome analysis
- Verify that your input sequence maps to a human gene
- Check your NCBI email configuration

**Error: "Genome file not found" after download attempt**
- Verify you have sufficient disk space (approximately 3 GB for genome file)
- Check internet connection stability
- Ensure write permissions in the application directory

**BLAST search takes very long or times out**

This is normal for the first run as NCBI queries require time to process.

If the search exceeds 5 minutes or appears stuck:
1. Cancel the process using Ctrl+C
2. Run again with your personal NCBI email address:
   ```bash
   # Windows
   python extractor.py -i sample.fasta -mail your.email@example.com
   
   # Linux and macOS
   python3 extractor.py -i sample.fasta -mail your.email@example.com
   ```
3. If the issue persists, create a free NCBI account for improved API priority:
   - Visit: https://www.ncbi.nlm.nih.gov/account/
   - Create an account and use your registered email with the `-mail` parameter
   - This significantly improves query speed and reduces timeouts

Once data is cached locally (JSON file stored in `raw_dir`), subsequent runs will be much faster and will not require the BLAST search step.

---

## Citation

If you use DNA Extractor in your research, please cite this repository.

---

## License

This project is provided for research use.

---

---

# DNA Extractor (Türkçe)

## Genel Bakış

DNA Extractor, insan genomu üzerinden transkripsiyon başlangıç sitesi (TSS) etrafındaki genişletilmiş DNA sekanslarını çıkarmak için tasarlanmış bir biyoinformatik aracıdır. Araç başta G4-quadruplex yapı analizi yapan araştırmacılar için geliştirilmiş olsa da, sunduğu esnek sekans çıkarma yetenekleri ve farklı filtreleme modları sayesinde çeşitli genomik araştırma uygulamalarında da kullanışlıdır.

Araç, protein sekansını girdi olarak alır, BLAST araması ve NCBI Entrez sorguları aracılığıyla karşılık gelen geni tanımlar ve transkripsiyon başlangıç sitesi etrafındaki genomik bölgeleri özelleştirilebilir şekilde çıkarır. Çıktı sekansları kapsamlı meta veriler ile FASTA formatında oluşturulur.

### Temel Özellikler

- Yalnızca insan genomu sekanslarını destekler (Homo sapiens, GRCh38)
- Ensembl entegrasyonu ile otomatik genom dosya yönetimi
- Sekans çıkarımı için birden fazla filtreleme modu
- Upstream ve downstream bölge boyutları özelleştirilebilir
- Alt analizler için referans veri setleri oluşturur
- Birden fazla transkript varyantını işler

### Önemli Not

Bu araç yalnızca insan genomu (Homo sapiens) için tasarlanmıştır. İnsan dışı türlerden gelen girdiler hata mesajı ile reddedilecektir.

---

## Hızlı Başlangıç: Sıfırdan Kurulum

### Adım 1: Miniconda İndirme ve Kurulumu

**Windows PowerShell:**
```powershell
# Miniconda yükleyicisini indir
Invoke-WebRequest -Uri "https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe" -OutFile "Miniconda3-installer.exe"

# Yükleyiciyi çalıştır
.\Miniconda3-installer.exe
```
Kurulum sırasında "Add Miniconda to PATH" seçeneğini işaretlediğinizden emin olun. Kurulumdan sonra PowerShell'i yeniden başlatın.

**Linux ve macOS Terminal:**
```bash
# Miniconda yükleyicisini indir
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Yükleyiciyi çalıştır
bash Miniconda3-latest-Linux-x86_64.sh
```
Talimatları izleyin ve lisans sözleşmesini kabul edin. Kurulumdan sonra terminalinizi yeniden başlatın.

Eğer `conda --version` hata veriyorsa, Miniconda'yı manuel olarak PATH'e ekleyin:
- **Windows:** `C:\Users\YourUsername\miniconda3\Scripts` adresini Sistem Ortam Değişkenleri PATH'ine ekleyin
- **Linux/macOS:** `echo 'export PATH="$HOME/miniconda3/bin:$PATH"' >> ~/.bashrc && source ~/.bashrc` komutunu çalıştırın

### Adım 2: DNA Extractor Deposunu İndirme

**Seçenek A: Git Kullanarak (Git yüklüyse)**
```bash
# Windows PowerShell, Linux ve macOS
git clone https://github.com/yourusername/dna_extractor.git
cd dna_extractor
```

**Seçenek B: ZIP İndir (Git gerekli değil)**
1. Şu adrese gidin: https://github.com/yourusername/dna_extractor
2. **Code** → **Download ZIP** seçeneğini tıklayın
3. ZIP dosyasını istediğiniz konuma çıkarın
4. Terminal/PowerShell açın ve klasöre gidin:
```bash
# Windows PowerShell
cd path\to\dna_extractor

# Linux ve macOS
cd path/to/dna_extractor
```

### Adım 3: Conda Ortamını Oluşturma

```bash
# Windows PowerShell
conda env create -f environment.yml
conda activate dna_extractor

# Linux ve macOS
conda env create -f environment.yml
conda activate dna_extractor
```

### Adım 4: Kurulumu Doğrulama

```bash
# Windows
python extractor.py -h

# Linux ve macOS
python3 extractor.py -h
```

Parametrelerin açıklandığı yardım mesajını görürseniz, kurulum başarılı olmuştur!

---

## Kurulum

### Ön Koşullar

- Miniconda (yukarıdaki Hızlı Başlangıç bölümü ile indirilmiş ve kurulmuş)
- Python 3.12 veya daha üstü (Miniconda ile birlikte gelir)
- İnternet bağlantısı (genom dosyası indirme ve NCBI sorguları için)

### Giriş Dosyası Hakkında Önemli Not (-i argümanı)

Betiği yeni bir protein sekansı ile ilk kez çalıştırırken `-i` argümanını kullanarak FASTA dosyasını sağlayın. Bu, BLAST aramasını ve Ensembl veri alımını tetikleyecektir. JSON dosyası oluşturulduktan ve `raw_dir` dizinine kaydedildikten sonra, sonraki çalıştırmalarda `-i` argümanına gerek yoktur. `-i` argümanını tekrar kullanmak verileri yeniden oluşturacaktır ki bu gereksiz ve zaman alıcıdır. `-i` argümanını yalnızca farklı bir protein sekansı analiz ederken kullanın.

---

## Kullanım

### Temel Komut Yapısı

```bash
# Windows
python extractor.py -i <girdi_fasta> -o <çıktı_adı> -up <upstream_bp> -down <downstream_bp> -mod <mod> -mail <email>

# Linux ve macOS
python3 extractor.py -i <girdi_fasta> -o <çıktı_adı> -up <upstream_bp> -down <downstream_bp> -mod <mod> -mail <email>
```

### Parametreler

- `-i, --input`: FASTA formatında protein sekansı (ilk çalıştırma için gerekli)
- `-o, --output`: Çıktı dosya adı (varsayılan: `<gene_id>_extended-TSS.fasta`)
- `-up, --upstream`: TSS'in upstream tarafındaki baz çifti sayısı (varsayılan: 1000)
- `-down, --downstream`: TSS'in downstream tarafındaki baz çifti sayısı (varsayılan: 500)
- `-mod, --mode`: Filtreleme modu - modlar bölümüne bakınız (varsayılan: all)
- `-mail, --email`: NCBI e-posta adresi (isteğe bağlı ancak hız sınırlandırmasını önlemek için önerilir)

### Çıkarım Modları

**all**: Tam genişletilmiş TSS sekansını filtreleme olmaksızın döndürür. Bu varsayılan moddur.

**canonical**: Yalnızca kanonik transkript sekansını çıkarır, varyant transkriptleri hariç tutar.

**no-repeats**: Sekansdan tekrarlayan elemanları (örneğin Alu, LINE-1) kaldırır. Bu elemanlar tipik olarak referans genomda küçük harflerle temsil edilir.

**no-exons**: Bölgedeki tüm transkriptlerden ekzon sekanslarını kaldırır. Her transkript için ayrı FASTA dosyası oluşturulur ve ekzon koordinatlarına göre hassas şekilde kaldırılır. Çıktı, her transkript için ayrı dosyalar içeren tar.gz arşivi olarak sağlanır.

### Örnek Komutlar

Varsayılan parametrelerle genişletilmiş TSS sekansı çıkarın:
```bash
# Windows
python extractor.py -i sample.fasta -o my_gene_sequences.fasta

# Linux ve macOS
python3 extractor.py -i sample.fasta -o my_gene_sequences.fasta
```

Özel upstream/downstream boyutları ile çıkarın:
```bash
# Windows
python extractor.py -i sample.fasta -o my_gene_sequences.fasta -up 2000 -down 1000

# Linux ve macOS
python3 extractor.py -i sample.fasta -o my_gene_sequences.fasta -up 2000 -down 1000
```

Yalnızca kanonik transkripti çıkarın:
```bash
# Windows
python extractor.py -i sample.fasta -mod canonical -mail your.email@example.com

# Linux ve macOS
python3 extractor.py -i sample.fasta -mod canonical -mail your.email@example.com
```

Tekrarlayan elemanları kaldırın:
```bash
# Windows
python extractor.py -i sample.fasta -mod no-repeats

# Linux ve macOS
python3 extractor.py -i sample.fasta -mod no-repeats
```

İntron sekanslarını çıkarın (ekzonları kaldırın):
```bash
# Windows
python extractor.py -i sample.fasta -mod no-exons -mail your.email@example.com

# Linux ve macOS
python3 extractor.py -i sample.fasta -mod no-exons -mail your.email@example.com
```

---

## Çıktı Dosyaları

Araç, çıkarılan sekans hakkında meta veriler içeren açıklayıcı başlıklar ile FASTA formatında çıktı dosyaları oluşturur:

- **all modu**: `<gene_id>_all.fasta` - Tam genişletilmiş TSS sekansı
- **canonical modu**: `<gene_id>_canonical.fasta` - Yalnızca kanonik transkript
- **no-repeats modu**: `<gene_id>_norepeats.fasta` - Tekrarlayan elemanlar kaldırılmış sekans
- **no-exons modu**: `<gene_id>_removed-exons.tar.gz` - Her bir transkript için ayrı FASTA dosyaları içeren sıkıştırılmış arşiv (ekzonlar kaldırılmış)

FASTA başlıkları aşağıdaki bilgileri içerir:
```
>gene_id|chromosome:start-end|strand:+/-|applied_mode
```

---

## Gereksinimler

Aşağıdaki paketler conda aracılığıyla otomatik olarak kurulur:

- pyfaidx >= 0.7.2
- biopython >= 1.81
- mygene >= 3.2.2
- requests >= 2.31.0

---

## Sorun Giderme

**Hata: "No JSON files found in the raw directory"**
- Girdi FASTA dosyanızın geçerli olduğundan ve uygun bir protein sekansı içerdiğinden emin olun
- NCBI ve Ensembl erişimi için internet bağlantısını kontrol edin
- BLAST sorgusu minimum uzunluk gereksinimini karşıladığından emin olun (en az 25 amino asit)

**Hata: "Only Homo sapiens is supported"**
- Bu araç yalnızca insan genomu analizi için tasarlanmıştır
- Girdi sekansının insan geni olarak eşlendiğini doğrulayın
- NCBI e-posta yapılandırmanızı kontrol edin

**Hata: "Genome file not found" indirme girişiminden sonra**
- Yeterli disk alanına sahip olduğunuzdan emin olun (genom dosyası için yaklaşık 3 GB)
- İnternet bağlantısı stabilitesini kontrol edin
- Uygulama dizininde yazma izinlerini doğrulayın

**BLAST araması çok uzun sürüyor veya zaman aşımı oluyor**

Bu, ilk çalıştırmada NCBI sorguları zaman aldığı için normaldir.

Arama 5 dakikayı aşarsa veya takılı kalırsa:
1. Ctrl+C kullanarak işlemi iptal edin
2. Kişisel NCBI e-posta adresiniz ile tekrar çalıştırın:
   ```bash
   # Windows
   python extractor.py -i sample.fasta -mail your.email@example.com
   
   # Linux ve macOS
   python3 extractor.py -i sample.fasta -mail your.email@example.com
   ```
3. Sorun devam ederse, geliştirilmiş API önceliği için ücretsiz NCBI hesabı oluşturun:
   - Ziyaret edin: https://www.ncbi.nlm.nih.gov/account/
   - Hesap oluşturun ve kayıtlı e-postanızı `-mail` parametresi ile kullanın
   - Bu, sorgu hızını önemli ölçüde iyileştirir ve zaman aşımlarını azaltır

## Alıntı

DNA Extractor'ı araştırmanızda kullanıyorsanız lütfen bu depoyu alıntı yapınız.

---

## Lisans

Bu proje araştırma amaçlı kullanım için sağlanmıştır.
