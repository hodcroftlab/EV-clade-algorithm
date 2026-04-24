import json, os, argparse, yaml

def extract_clades(json_file, output_dir, clade_key="new-clade"):
    with open(json_file) as fh:
        data = json.load(fh)

    os.makedirs(output_dir, exist_ok=True)

    def walk(n):
        if "branch_attrs" in n and "labels" in n["branch_attrs"]:
            label = n["branch_attrs"]["labels"].get(clade_key)
            if label:
                muts = []
                for gene, changes in n["branch_attrs"].get("mutations", {}).items():
                    if gene == "nuc": continue
                    for m in changes:
                        try:
                            pos = int(m[1:-1])
                            state = m[-1]
                            muts.append({"locus": gene, "position": pos, "state": state})
                        except:
                            continue
                yml = {
                    "name": label,
                    "unaliased_name": label,
                    "parent": ".".join(label.split(".")[:-1]) or "none",
                    "representatives": [],
                    "defining_mutations": muts
                }
                with open(os.path.join(output_dir, f"{label}.yml"), "w") as out:
                    yaml.dump(yml, out, sort_keys=False)
        for c in n.get("children", []):
            walk(c)

    walk(data["tree"])

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", required=True, help="Path to suggested_*.json")
    parser.add_argument("--outdir", default="subclades", help="Output directory for .yml files")
    args = parser.parse_args()
    extract_clades(args.json, args.outdir)
