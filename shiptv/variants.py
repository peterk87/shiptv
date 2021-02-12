import re

aa_codes = dict(
    ALA='A',
    ARG='R',
    ASN='N',
    ASP='D',
    CYS='C',
    GLU='E',
    GLN='Q',
    GLY='G',
    HIS='H',
    ILE='I',
    LEU='L',
    LYS='K',
    MET='M',
    PHE='F',
    PRO='P',
    SER='S',
    THR='T',
    TRP='W',
    TYR='Y',
    VAL='V',
)


def parse_aa(gene: str,
             ref: str,
             alt: str,
             nt_pos: int,
             aa_pos: int,
             snpeff_aa: str) -> str:
    m = re.match(r'p\.([a-zA-Z]+)(\d+)([a-zA-Z]+)', snpeff_aa)
    if snpeff_aa == '.' or m is None:
        return f'{ref}{nt_pos}{alt}'
    ref_aa, aa_pos_str, alt_aa = m.groups()
    ref_aa = get_aa(ref_aa)
    alt_aa = get_aa(alt_aa)
    return f'{ref}{nt_pos}{alt}({gene}:{ref_aa}{aa_pos}{alt_aa})'


def get_aa(s: str) -> str:
    out = ''
    for i in range(0, len(s), 3):
        aa = s[i: i + 3]
        aa_code = aa_codes[aa.upper()]
        out += aa_code
    return out
