import os
import datetime
import dataclasses
from typing import List
import math

@dataclasses.dataclass
class GenBankLabel:
    name: str
    start: int
    end: int
    is_forward: bool

    def __post__init__(self):
        assert(self.start <= self.end)

@dataclasses.dataclass
class GenBank:
    seq: str
    labels: List[GenBankLabel] = []
    
    def __post_init__(self):
        self.seq = self.seq.strip()

    def save(self, fpath, overwrite=False):
        if os.path.exists(fpath) and overwrite is False:
            raise FileExistsError(f"{fpath}")
        
        today = datetime.date.today()

        with open(fpath, 'w') as f:
            lines = []
            lines.append(f"LOCUS       genbank.py                     {len(self.seq)} bp    DNA        linear       {today.strftime(r'%d-%b-%Y').upper()}")
            lines.append("FEATURES             Location/Qualifiers")
            # label
            for label in self.labels:
                if label.is_forward:
                    lines.append(f"     misc_feature    {label.start}..{label.end}")
                else:
                    lines.append(f"     misc_feature    complement({label.start}..{label.end})")
                lines.append(f'                     /label="{label.name}"')
            # seq
            lines.append("ORIGIN")
            letters_per_line, letters_per_chunk = 60, 10
            num_lines = math.ceil(len(self.seq) / letters_per_line)
            for ln in range(num_lines):
                line = f"{ln * letters_per_line + 1:>9d} "
                left_seq = self.seq[ln * letters_per_line:]
                num_letters = min(letters_per_line, len(left_seq))
                line += " ".join([left_seq[i * letters_per_chunk: (i + 1) * letters_per_chunk] for i in range(math.ceil(num_letters / letters_per_chunk))])
                lines.append(line)
            lines.append("//")
            f.write("\n".join(lines))
