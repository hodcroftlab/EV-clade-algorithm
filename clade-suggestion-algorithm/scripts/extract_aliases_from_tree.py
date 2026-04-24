"""
Extract starter aliases from Nextstrain tree.
Converts flat clade names to hierarchical tuples using regex splitting.
Universal approach: converts ordered letter suffixes (a, b, c) to numeric indices.
"""
import json
import argparse
import re


def is_ordered_letter_suffix(letter):
    """
    Check if letter is a simple alphabetic progression (a, b, c, d, ...).
    Returns the numeric index (1 for 'a', 2 for 'b', etc.) or None.
    """
    if len(letter) == 1 and letter.islower():
        return ord(letter) - ord('a') + 1
    return None


def parse_clade_string(clade_str, existing_clades=None):
    """
    Convert clade string to hierarchical list.
        
    Args:
        clade_str: clade name string
        existing_clades: set of tuples already assigned (for conflict avoidance)
    """
    if existing_clades is None:
        existing_clades = set()
    
    if not clade_str or clade_str == "unassigned":
        return "unassigned"
    
    # Skip "pre-" prefix clades
    if "pre" in clade_str.lower():
        return "unassigned"
    
    if "rfs" in clade_str.lower():
        return "unassigned"
    
    # Handle dot-separated format first (e.g., "A.1.2")
    if "." in clade_str:
        parts = clade_str.split(".")
        hierarchy = [parts[0]]
        
        for part in parts[1:]:
            # Check for suffix like "2.r" → treat "r" as special case
            if part.lower() == "r":
                # Recombinant: assign next available number
                return normalize_special_suffix(hierarchy, "r", existing_clades)
            elif re.match(r'^[a-z]$', part):
                # Single letter suffix: check if it's ordered
                letter_idx = is_ordered_letter_suffix(part)
                if letter_idx:
                    hierarchy.append(letter_idx)
                else:
                    hierarchy.append(part)
            else:
                try:
                    hierarchy.append(int(part))
                except ValueError:
                    hierarchy.append(part)
        
        return hierarchy
    
    # Check for special suffixes first ("-like", ".r")
    suffix_match = re.search(r'(-like|\.?r)$', clade_str, re.IGNORECASE)
    if suffix_match:
        suffix = suffix_match.group().lower()
        base_str = clade_str[:suffix_match.start()]
        
        base_hierarchy = parse_clade_string(base_str, existing_clades)
        if base_hierarchy:
            return normalize_special_suffix(base_hierarchy, suffix, existing_clades)
        else:
            return None
    
    # Check for ordered letter suffix (a, b, c, etc.) - but only if preceded by number
    letter_suffix_match = re.search(r'([a-z])$', clade_str)
    if letter_suffix_match:
        letter = letter_suffix_match.group(1)
        base_str = clade_str[:letter_suffix_match.start()]
        
        # Only convert if base ends with a digit (e.g., "B1a", not "Ba")
        if re.search(r'\d$', base_str):
            letter_idx = is_ordered_letter_suffix(letter)
            if letter_idx:
                base_hierarchy = parse_clade_string(base_str, existing_clades)
                if base_hierarchy:
                    # Check for collision; if it exists, skip conversion
                    candidate = tuple(base_hierarchy + [letter_idx])
                    if candidate in existing_clades:
                        # Collision; keep as string
                        return base_hierarchy + [letter]
                    else:
                        return base_hierarchy + [letter_idx]
    
    # Strip non-alphanumeric suffixes (e.g., "A2/D" → "A2")
    clade_clean = re.match(r'^[a-zA-Z0-9]+', clade_str)
    if not clade_clean:
        return None
    
    clade_clean = clade_clean.group()
    
    # Split on letter→number boundaries only (not number→letter)
    tokens = re.split(r'(?<=[a-zA-Z])(?=\d)', clade_clean)
    tokens = [t for t in tokens if t]
    
    if not tokens:
        return None
    
    # Convert numeric tokens to int, keep letters as strings
    hierarchy = []
    for token in tokens:
        try:
            hierarchy.append(int(token))
        except ValueError:
            hierarchy.append(token)
    
    return hierarchy if hierarchy else None


def normalize_special_suffix(base_clade, suffix, existing_clades):
    """
    Convert special suffixes ("-like", "r") to next available numeric index.
    
    Args:
        base_clade: list (e.g., ["C", 1])
        suffix: string ("-like", "r", ".r")
        existing_clades: set of tuples to avoid collisions
    
    Returns:
        base_clade + [next_available_idx]
    """
    idx = 1
    candidate = tuple(base_clade + [idx])
    
    while candidate in existing_clades:
        idx += 1
        candidate = tuple(base_clade + [idx])
    
    return base_clade + [idx]


def extract_clade_hierarchy(tree, clade_key="clade_membership"):
    """
    Extract clade labels from tree and convert to hierarchical format.
    Two-pass: first collect all clades, then normalize with collision awareness.
    
    Returns:
        dict: {short_name: hierarchical_list, ...}
    """
    # First pass: collect all clades
    raw_clades = {}
    existing_hierarchies = set()
    
    def walk_tree(node):
        if "node_attrs" in node and clade_key in node["node_attrs"]:
            clade_value = node["node_attrs"][clade_key]["value"]
            
            if clade_value and clade_value != "unassigned":
                raw_clades[clade_value] = None
        
        if "children" in node:
            for child in node["children"]:
                walk_tree(child)
    
    walk_tree(tree)
    
    # Second pass: parse clades, tracking existing hierarchies to avoid collisions
    aliases = {}
    for clade_str in raw_clades.keys():
        hierarchy = parse_clade_string(clade_str, existing_hierarchies)
        if hierarchy:
            aliases[clade_str] = hierarchy
            existing_hierarchies.add(tuple(hierarchy))
    
    return aliases


def build_aliases_json(tree, clade_key="clade_membership"):
    """
    Extract and return aliases, sorted by depth and alphabetically.
    """
    raw_aliases = extract_clade_hierarchy(tree, clade_key)
    
    # Sort: by tuple length (depth), then alphabetically
    sorted_aliases = {}
    for short_name in sorted(
        raw_aliases.keys(),
        key=lambda x: (len(parse_clade_string(x)), x)
    ):
        sorted_aliases[short_name] = raw_aliases[short_name]
    
    return sorted_aliases


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract starter aliases from Nextstrain tree"
    )
    parser.add_argument("--tree", type=str, required=True, help="Auspice JSON tree")
    parser.add_argument("--clade-key", type=str, default="clade_membership")
    parser.add_argument("--output", type=str, default="aliases.json")
    
    args = parser.parse_args()
    
    with open(args.tree) as fh:
        data = json.load(fh)
    
    aliases = build_aliases_json(data["tree"], args.clade_key)
    
    with open(args.output, "w") as fh:
        json.dump(aliases, fh, indent=2)
    
    print(f"✓ Extracted {len(aliases)} aliases → {args.output}")
    print(json.dumps(aliases, indent=2))