import os, csv, gzip
from snakemake.shell import shell

def is_fastq_file(filepath):
    """Check if file is a valid FASTQ file by reading first few lines"""
    try:
        # Handle both gzipped and uncompressed files
        opener = gzip.open if filepath.endswith('.gz') else open
        mode = 'rt' if filepath.endswith('.gz') else 'r'
        
        with opener(filepath, mode) as f:
            first_line = f.readline().strip()
            if not first_line.startswith('@'):
                return False
            
            # Read sequence line
            seq_line = f.readline().strip()
            if not seq_line:
                return False
                
            # Read plus line
            plus_line = f.readline().strip()
            if not plus_line.startswith('+'):
                return False
                
            # Read quality line
            qual_line = f.readline().strip()
            if not qual_line or len(qual_line) != len(seq_line):
                return False
                
        return True
    except:
        return False

def get_file_size_mb(filepath):
    """Get file size in MB"""
    try:
        return round(os.path.getsize(filepath) / (1024 * 1024), 1)
    except:
        return 0

def main(input_file, output_file):
    errors = []
    warnings = []
    seen = set()
    total_samples = 0
    total_size_mb = 0

    with open(input_file, newline="") as f:
        r = csv.reader(f)

        # Header überspringen falls vorhanden
        header = next(r, None)
        if not header or len(header) < 3 or not header[0].lower().startswith("sample"):
            # kein valider Header -> zurück zum Start
            f.seek(0)
            r = csv.reader(f)

        for i, row in enumerate(r, start=2):
            if len(row) < 3:
                errors.append(f"Line {i}: expected 3 columns, got {len(row)} -> {row}")
                continue

            sample, left, right = (c.strip() for c in row[:3])

            if not sample or not left or not right:
                errors.append(f"{sample or f'Line {i}'}: empty field(s)")
                continue

            if sample in seen:
                errors.append(f"{sample}: duplicate sample_id")
            else:
                seen.add(sample)
                total_samples += 1

            if left == right:
                errors.append(f"{sample}: left and right are identical -> {left}")

            # Performance optimized: use os.path.exists instead of os.path.isfile
            if not os.path.exists(left):
                errors.append(f"{sample}: left file missing -> {left}")
            elif not is_fastq_file(left):
                errors.append(f"{sample}: left file is not valid FASTQ -> {left}")
            else:
                left_size = get_file_size_mb(left)
                total_size_mb += left_size
                if left_size == 0:
                    warnings.append(f"{sample}: left file is empty -> {left}")
                elif left_size < 1:
                    warnings.append(f"{sample}: left file very small ({left_size}MB) -> {left}")
                    
            if not os.path.exists(right):
                errors.append(f"{sample}: right file missing -> {right}")
            elif not is_fastq_file(right):
                errors.append(f"{sample}: right file is not valid FASTQ -> {right}")
            else:
                right_size = get_file_size_mb(right)
                total_size_mb += right_size
                if right_size == 0:
                    warnings.append(f"{sample}: right file is empty -> {right}")
                elif right_size < 1:
                    warnings.append(f"{sample}: right file very small ({right_size}MB) -> {right}")

    with open(output_file, "w") as out:
        if errors:
            out.write("❌ Input checks failed:\n" + "\n".join(f"- {e}" for e in errors) + "\n")
        else:
            out.write(f"✅ All sample checks passed.\n")
            out.write(f"📊 Summary: {total_samples} samples, {total_size_mb:.1f}MB total\n")
            
        if warnings:
            out.write("\n⚠️  Warnings:\n" + "\n".join(f"- {w}" for w in warnings) + "\n")

if __name__ == "__main__":
    main(snakemake.input.samples_csv, snakemake.output[0])
