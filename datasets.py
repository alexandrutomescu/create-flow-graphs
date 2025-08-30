import os, sys
from math import ceil
import gzip
from typing import List, Tuple

# Loader returns (genome_dir, genome_files, fixed_abundances_or_None)
# fixed_abundances_or_None is a list of floats (length = ngenomes) or None if distribution should be used later.

ALLOWED_FASTA_SUFFIXES = ('.fa', '.fna', '.fasta')
ALLOWED_FASTA_SUFFIXES_GZ = tuple(s + '.gz' for s in ALLOWED_FASTA_SUFFIXES) + ALLOWED_FASTA_SUFFIXES

def _list_fasta_files(directory: str) -> List[str]:
    return sorted([f for f in os.listdir(directory) if f.endswith(ALLOWED_FASTA_SUFFIXES_GZ)])

def _resolve_fname_with_gz(directory: str, fname: str) -> Tuple[str, str|None]:
    """Return existing filename (possibly switching to .gz) and optional info message."""
    base_path = os.path.join(directory, fname)
    if os.path.isfile(base_path):
        return fname, None
    # try gzip variant if not already gz
    if not fname.endswith('.gz') and os.path.isfile(base_path + '.gz'):
        return fname + '.gz', f'INFO: Using compressed file {fname}.gz instead of missing {fname}'
    # try removing .gz if provided
    if fname.endswith('.gz'):
        uncompressed = fname[:-3]
        if os.path.isfile(os.path.join(directory, uncompressed)):
            return uncompressed, f'INFO: Using uncompressed file {uncompressed} instead of missing {fname}'
    return fname, None  # will be caught by caller if truly missing

def load_dataset(name: str, ngenomes: int, user_abundances, allow_distribution: bool=True):
    """Load dataset metadata.

    name: one of ecoli, labmix, complex32, medium20, helicobacter-hepaticus
    ngenomes: number requested (must not exceed available)
    user_abundances: list[float] or None from -A
    allow_distribution: if False distribution will be ignored (when fixed abundances provided)

    Returns (genome_dir, genome_files, abundances_list_or_None, info_messages)
    info_messages: list of strings to print (INFO / WARNING messages)
    """
    msgs = []
    if name == 'ecoli':
        genome_dir = os.path.join('datasets','ecoli')
        if not os.path.isdir(genome_dir):
            print(f'ERROR: missing directory {genome_dir}')
            sys.exit(1)
        genome_files = _list_fasta_files(genome_dir)
        if ngenomes > len(genome_files):
            print(f'ERROR: --ngenomes ({ngenomes}) exceeds available ecoli genomes ({len(genome_files)})')
            sys.exit(1)
        genome_files = genome_files[:ngenomes]
        # If user provided abundances, we pass them back (validated in caller)
        return genome_dir, genome_files, user_abundances, msgs

    if name == 'helicobacter-hepaticus':
        genome_dir = os.path.join('datasets','helicobacter-hepaticus')
        if not os.path.isdir(genome_dir):
            print(f'ERROR: missing directory {genome_dir}')
            sys.exit(1)
        # Accept common compressed FASTA suffixes
        allowed_suffixes = ('.fa.gz', '.fna.gz', '.fasta.gz')
        genome_files = sorted([f for f in os.listdir(genome_dir) if f.endswith(allowed_suffixes)])
        if len(genome_files) == 0:
            print(f'ERROR: No compressed FASTA (.fa.gz/.fna.gz/.fasta.gz) files found in {genome_dir}')
            sys.exit(1)
        if ngenomes > len(genome_files):
            print(f'ERROR: --ngenomes ({ngenomes}) exceeds available helicobacter-hepaticus genomes ({len(genome_files)})')
            sys.exit(1)
        genome_files = genome_files[:ngenomes]
        # Behaves like ecoli: use user-provided abundances if given, else distribution later
        return genome_dir, genome_files, user_abundances, msgs

    if name == 'labmix':
        genome_dir = os.path.join('datasets','labmix')
        labmix_order = ["896.fasta", "HXB2.fasta", "JRCSF.fasta", "NL43.fasta", "YU2.fasta"]
        # Resolve possible gzip variants
        resolved = []
        for f in labmix_order:
            rf, info = _resolve_fname_with_gz(genome_dir, f)
            if info:
                msgs.append(info)
            resolved.append(rf)
        missing = [f for f in resolved if not os.path.isfile(os.path.join(genome_dir, f))]
        if missing:
            print(f'ERROR: Missing expected labmix genome files: {missing}')
            sys.exit(1)
        if ngenomes > len(resolved):
            print(f'ERROR: --ngenomes ({ngenomes}) exceeds available labmix genomes ({len(resolved)})')
            sys.exit(1)
        genome_files = resolved[:ngenomes]
        fixed_abundances = [17.0, 6.0, 27.0, 36.0, 15.0][:ngenomes]
        if user_abundances is not None:
            msgs.append('INFO: Ignoring --abundances for labmix; using fixed abundances 17,6,27,36,15 (truncated).')
        return genome_dir, genome_files, fixed_abundances, msgs

    if name == 'complex32':
        genome_dir = os.path.join('datasets','complex32')
        abundance_file = os.path.join(genome_dir, 'nanosim.abundances.tsv')
        if not os.path.isfile(abundance_file):
            print(f'ERROR: Missing abundance file {abundance_file}')
            sys.exit(1)
        parsed = []
        with open(abundance_file) as af:
            for line in af:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) != 2:
                    msgs.append(f'WARNING: Malformed line in abundances TSV ignored: {line}')
                    continue
                fname, abund = parts
                try:
                    abund_val = float(abund)
                except ValueError:
                    msgs.append(f'WARNING: Non-numeric abundance for {fname}: {abund}; skipping')
                    continue
                actual_fname, info = _resolve_fname_with_gz(genome_dir, fname)
                if info:
                    msgs.append(info)
                if not os.path.isfile(os.path.join(genome_dir, actual_fname)):
                    msgs.append(f'WARNING: Listed genome missing in directory (even after .gz fallback): {fname}')
                    continue
                parsed.append((actual_fname, abund_val))
        if len(parsed) == 0:
            print('ERROR: No valid genomes found for complex32 dataset.')
            sys.exit(1)
        if ngenomes > len(parsed):
            print(f'ERROR: --ngenomes ({ngenomes}) exceeds available complex32 genomes ({len(parsed)})')
            sys.exit(1)
        genome_files = [fname for fname, _ in parsed][:ngenomes]
        # Scale abundances by 10 and round to nearest int as requested
        abundances = [int(round(abund * 10)) for _, abund in parsed][:ngenomes]
        if user_abundances is not None:
            msgs.append('INFO: Ignoring --abundances for complex32; using (scaled x10, rounded) abundances from nanosim.abundances.tsv (truncated).')
        return genome_dir, genome_files, abundances, msgs

    if name == 'medium20':
        genome_dir = os.path.join('datasets','medium20')
        abundance_file = os.path.join(genome_dir, 'nanosim.abundances.tsv')
        if not os.path.isfile(abundance_file):
            print(f'ERROR: Missing abundance file {abundance_file}')
            sys.exit(1)
        parsed = []
        with open(abundance_file) as af:
            for line in af:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) != 2:
                    msgs.append(f'WARNING: Malformed line in abundances TSV ignored: {line}')
                    continue
                fname, abund = parts
                try:
                    abund_val = float(abund)
                except ValueError:
                    msgs.append(f'WARNING: Non-numeric abundance for {fname}: {abund}; skipping')
                    continue
                actual_fname, info = _resolve_fname_with_gz(genome_dir, fname)
                if info:
                    msgs.append(info)
                if not os.path.isfile(os.path.join(genome_dir, actual_fname)):
                    msgs.append(f'WARNING: Listed genome missing in directory (even after .gz fallback): {fname}')
                    continue
                parsed.append((actual_fname, abund_val))
        if len(parsed) == 0:
            print('ERROR: No valid genomes found for medium20 dataset.')
            sys.exit(1)
        if ngenomes > len(parsed):
            print(f'ERROR: --ngenomes ({ngenomes}) exceeds available medium20 genomes ({len(parsed)})')
            sys.exit(1)
        genome_files = [fname for fname, _ in parsed][:ngenomes]
        abundances = [int(round(abund * 10)) for _, abund in parsed][:ngenomes]
        if user_abundances is not None:
            msgs.append('INFO: Ignoring --abundances for medium20; using (scaled x10, rounded) abundances from nanosim.abundances.tsv (truncated).')
        return genome_dir, genome_files, abundances, msgs

    print('ERROR: unknown dataset')
    sys.exit(1)
