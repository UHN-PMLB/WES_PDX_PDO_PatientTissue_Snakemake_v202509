#!/usr/bin/env python3
import sys

vcf_in = sys.argv[1]
driver_file = sys.argv[2]
vcf_out = sys.argv[3]

# Load driver genes
drivers = set(line.strip() for line in open(driver_file) if line.strip())

functional_terms = {
    "missense_variant", "stop_gained", "frameshift_variant", "start_lost",
    "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant",
    "inframe_insertion", "inframe_deletion"
}

csq_idx = {}
seen = set()

with open(vcf_in) as fin, open(vcf_out, "w") as fout:
    for line in fin:
        if line.startswith("##INFO=<ID=CSQ"):
            desc = line.split("Format: ")[1].rstrip('"> \n')
            fields = desc.split("|")
            csq_idx = {name: i for i, name in enumerate(fields)}

        if line.startswith("#"):
            fout.write(line)
            continue

        fields = line.rstrip("\n").split("\t")
        info = fields[7]

        if "CSQ=" not in info:
            continue

        csq_entries = info.split("CSQ=")[1].split(";")[0].split(",")

        chrom, pos, _, ref, alt = fields[:5]
        key = f"{chrom}:{pos}:{ref}:{alt}"

        for entry in csq_entries:
            parts = entry.split("|")

            gene = parts[csq_idx["SYMBOL"]]
            consequences = parts[csq_idx["Consequence"]].split("&")

            if gene in drivers and any(c in functional_terms for c in consequences):
                if key not in seen:
                    fout.write(line)
                    seen.add(key)
                break

