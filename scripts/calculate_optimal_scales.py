#!/usr/bin/env python3
"""
Script to calculate optimal normalization scales for EV-D68 clade suggestion algorithm.

This script analyzes a phylogenetic tree to determine empirically-based scale parameters
for the three main normalization factors used in clade designation:

1. bushiness_branch_scale: Based on average nucleotide mutations per branch
    We define s_ϕ to be the median value of ϕ across the tree.
2. branch_length_scale: Based on average weighted amino acid mutations per branch  
3. divergence_scale: Based on divergence patterns within and between clades

Author: Adapted for EV-D68 from influenza clade suggestion algorithm
Usage: python calculate_optimal_scales.py --tree auspice/original_D68.json --weights weights.json --output scales.json
"""

import json
import argparse
import numpy as np
from collections import defaultdict
import pandas as pd
       
def parse_tree_mutations(tree, weights=None): 
    """
    Parse the tree to extract mutation information from all branches.
    
    Args:
        tree: Tree node (dict) from Nextstrain JSON
        weights: Dictionary of mutation weights per position
    
    Returns:
        tuple: (nt_mutations_per_branch, aa_mutations_per_branch, all_branches)
    """
    if weights is None:
        weights = {}
    
    nt_mutations = []  # nucleotide mutations per branch
    aa_mutations = []  # weighted amino acid mutations per branch
    all_branches = []  # store all branch information
    
    def traverse_tree(node, parent_name="root"):
        """ 
        Recursively traverse tree and collect mutation data - starting at root of the tree 
        """
        
        # Get mutations for this branch
        mutations = node['branch_attrs']['mutations']
        
        # Count nucleotide mutations (excluding gaps and Ns)
        nt_count = len([m for m in mutations.get('nuc', [])
                       if m[0] not in ['N', '-'] and m[-1] not in ['N', '-']])

          
        # Count weighted amino acid mutations across all proteins
        aa_weight = 0
        for protein, protein_muts in mutations.items():
            if protein == 'nuc': continue  # Skip nucleotide mutations
            w= weights.get(protein, {})            
            for mut in protein_muts:
                if mut[-1] in ['-', 'X']: continue # Skip deletions and unknown AAs
                pos = str(mut[1:-1])  # Extract position number
                aa_weight += w[pos] if pos in w else w.get("default", 1)

        if node["name"] == 'NODE_0000000' and nt_count > 0:
            # ipdb.set_trace()  # Debugging breakpoint
            print(f"""\033[91m
    The root node has {nt_count} nucleotide mutations. Was the tree exported from Nextclade?
    Root mutations will be excluded here.
            \033[0m""")
            nt_count = 0  # Reset root node mutations to avoid double counting
            aa_weight = 0  # Reset root node amino acid mutations

        branch_info = {
            'node_name': node.get('name', f"internal_{len(all_branches)}"),
            'parent': parent_name,
            'nt_mutations': nt_count,
            'aa_weight': aa_weight,
            'total_proteins': len([p for p in mutations.keys() if p != 'nuc'])
        }
    
        all_branches.append(branch_info)
        nt_mutations.append(nt_count)
        aa_mutations.append(aa_weight)
        
        # Recursively process children
        for child in node.get('children', []):
            traverse_tree(child, node.get('name', 'internal'))
    
    traverse_tree(tree)
    return nt_mutations, aa_mutations, all_branches


def calculate_clade_divergences(tree, clade_key='clade_membership'): 
    """
    Calculate divergence patterns within and between existing clades.
    
    Args:
        tree: Tree root node
        clade_key: Key for existing clade annotations
    
    Returns:
        tuple: (within_clade_divs, between_clade_divs)
    """
    # First, collect all nodes with their clade assignments and divergence
    nodes_with_clades = []
    
    def collect_nodes(node, current_div=0):
        """Collect nodes with clade info and cumulative divergence"""
        # ipdb.set_trace()  # Debugging breakpoint
        # Get clade assignment
        clade = node['node_attrs'][clade_key].get('value', 'unassigned') 
        # Calculate divergence (mutations since root)
        mutations = node['branch_attrs']['mutations']
        # Count AA mutations across all proteins for divergence calculation
        branch_div = 0
        for protein, muts in mutations.items():
            if protein != 'nuc':  # Only count amino acid changes
                branch_div += len([m for m in muts if m[-1] not in ['-', 'X']])
        
        total_div = current_div + branch_div
        
        # Store node info
        node_info = {
            'name': node.get('name', 'internal'),
            'clade': clade,
            'divergence': total_div,
            'is_terminal': 'children' not in node or len(node['children']) == 0
        }
        nodes_with_clades.append(node_info)
        
        # Recursively process children
        for child in node.get('children', []):
            collect_nodes(child, total_div)
    
    collect_nodes(tree)
    
    # Group nodes by clade
    clade_groups = defaultdict(list)
    for node in nodes_with_clades:
        if node['clade'] != 'unassigned':
            clade_groups[node['clade']].append(node)
    
    # Calculate within-clade divergences
    within_clade_divs = []
    for clade, nodes in clade_groups.items():
        if len(nodes) > 1:
            divs = [n['divergence'] for n in nodes]
            # Calculate pairwise differences within clade
            for i in range(len(divs)):
                for j in range(i+1, len(divs)):
                    within_clade_divs.append(abs(divs[i] - divs[j]))
    
    # Calculate between-clade divergences
    between_clade_divs = []
    clade_names = list(clade_groups.keys())
    for i in range(len(clade_names)):
        for j in range(i+1, len(clade_names)):
            clade1_divs = [n['divergence'] for n in clade_groups[clade_names[i]]]
            clade2_divs = [n['divergence'] for n in clade_groups[clade_names[j]]]
            
            # Calculate pairwise differences between clades
            for div1 in clade1_divs:
                for div2 in clade2_divs:
                    between_clade_divs.append(abs(div1 - div2))

    return within_clade_divs, between_clade_divs


def calculate_optimal_scales(tree, weights=None, key='clade_membership'):
    """
    Calculate optimal scale parameters based on empirical tree data.
    
    The scales are chosen to normalize scores such that:
    - bushiness_branch_scale: represents typical nucleotide mutation load
    - branch_length_scale: represents typical weighted AA mutation load  
    - divergence_scale: distinguishes within vs between clade divergence
    
    Args:
        tree: Nextstrain tree JSON
        weights: Mutation weights dictionary
    
    Returns:
        dict: Optimal scale parameters with explanations
    """
    print("Analyzing tree structure for optimal scale calculation...")
    # Give number of tips (not nodes)
    count_tips = lambda tree: 1 if "children" not in tree or not tree["children"] else sum(count_tips(child) for child in tree["children"])
    num_tips = count_tips(tree)
    print(f"Total tips in tree: {num_tips}")
    # ipdb.set_trace()  # Debugging breakpoint
    
    # 1. Calculate bushiness_branch_scale (nucleotide mutations)
    print("\n1. Calculating bushiness_branch_scale...")
    nt_muts, aa_weights, branches = parse_tree_mutations(tree, weights)

    
    
    # Remove branches with zero mutations to avoid skewing averages
    nt_muts_nonzero = [x for x in nt_muts if x > 0]
    
    # Use median as it's more robust to outliers than mean
    bushiness_scale = np.median(nt_muts_nonzero) if nt_muts_nonzero else 1.0
    range_bushiness = np.percentile(nt_muts_nonzero,[40,60])

    print(f"   - Total number of tips: {num_tips}")
    print(f"   - Total branches analyzed: {len(nt_muts)}")
    print(f"   - Branches with mutations: {len(nt_muts_nonzero)}")
    print(f"   - Mean NT mutations per branch: {np.mean(nt_muts_nonzero):.2f}")
    print(f"   - Median NT mutations per branch: {np.median(nt_muts_nonzero):.2f}")
    print(f"   - Recommended bushiness_branch_scale: {bushiness_scale:.1f}")
    
    # 2. Calculate branch_length_scale (weighted AA mutations)
    print("\n2. Calculating branch_length_scale...")
    # ipdb.set_trace()  # Debugging breakpoint
    aa_weights_nonzero = [x for x in aa_weights if x > 0]

    branch_length_scale = np.percentile(aa_weights_nonzero,90)
    range_branch = np.percentile(aa_weights_nonzero,[85,95])
        
    print(f"   - Branches with AA mutations: {len(aa_weights_nonzero)}")
    print(f"   - Mean weighted AA mutations per branch: {np.mean(aa_weights_nonzero):.2f}")
    print(f"   - Median weighted AA mutations per branch: {np.median(aa_weights_nonzero):.2f}")
    print(f"   - Recommended branch_length_scale: {branch_length_scale:.2f}")
    
    # 3. Calculate divergence_scale (within vs between clade divergence)
    print("\n3. Calculating divergence_scale...")
    within_divs, between_divs = calculate_clade_divergences(tree, key)
    
    if within_divs and between_divs:
        # Use the difference between median between-clade and within-clade divergence
        # This represents the typical divergence that separates clades
        median_within = np.median(within_divs)
        median_between = np.median(between_divs)
        divergence_scale = median_between - median_within
        
        # Ensure positive scale
        divergence_scale = max(divergence_scale, 1.0)
        range_divergence = [min(median_within, 4.0),median_between]
        
        print(f"   - Within-clade divergence pairs: {len(within_divs)}")
        print(f"   - Between-clade divergence pairs: {len(between_divs)}")
        print(f"   - Median within-clade divergence: {median_within:.2f}")
        print(f"   - Median between-clade divergence: {median_between:.2f}")
        print(f"   - Recommended divergence_scale: {divergence_scale:.1f}")
    else:
        # Fallback if no clades are present
        divergence_scale = 4.0  # Default from original algorithm
        print(f"   - No existing clades found, using default: {divergence_scale}")
    
    # Compile results
    results = {
        "optimal_scales": {
            "bushiness_branch_scale": round(bushiness_scale, 1),
            "branch_length_scale": round(branch_length_scale, 1), 
            "divergence_scale": round(divergence_scale, 1)
        },
        "empirical_data": {
            "Tree size (no. tips)": num_tips,
            "Total no. of branches": len(nt_muts),
            "Branches with NT mutations": len(nt_muts_nonzero),
            "Branches with AA mutations": len(aa_weights_nonzero),
            "Mean NT muts per branch": round(np.mean(nt_muts_nonzero), 2) if nt_muts_nonzero else 0,
            "Median NT muts per branch": round(np.median(nt_muts_nonzero), 2) if nt_muts_nonzero else 0,
            "Mean weighted AA muts per branch": round(np.mean(aa_weights_nonzero), 2) if aa_weights_nonzero else 0,
            "Median weighted AA muts per branch": round(np.median(aa_weights_nonzero), 2) if aa_weights_nonzero else 0,
            "Within-clade divergence pairs": len(within_divs) if within_divs else 0,
            "Between-clade divergence pairs": len(between_divs) if between_divs else 0,
            "Median within-clade divergence": round(np.median(within_divs), 2) if within_divs else None,
            "Median between-clade divergence": round(np.median(between_divs), 2) if between_divs else None
        },
        "explanation": {
            "bushiness_branch_scale": "Based on median nucleotide mutations per branch. This normalizes the phylogenetic bushiness score.",
            "branch_length_scale": "Based on median weighted amino acid mutations per branch. This normalizes the branch mutation score.",
            "divergence_scale": "Based on difference between median inter- and intra-clade divergence. This helps distinguish meaningful clade boundaries."
        }
    }
    
    return results, [range_bushiness, range_branch, range_divergence]


def main():
    """Main function to run the scale calculation"""
    parser = argparse.ArgumentParser(
        description="Calculate optimal normalization scales for EV-D68 clade suggestion algorithm",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    
    parser.add_argument('--tree', required=True, 
                       help='Input Nextstrain tree JSON file')
    parser.add_argument('--weights', 
                       help='Optional weights JSON file for amino acid mutations')
    parser.add_argument('--output', required=True,
                       help='Output JSON file for calculated scales')
    parser.add_argument('--clade-key', default='clade_membership',
                       help='Key for existing clade annotations in tree (default: clade_membership)')
    
    args = parser.parse_args()
    
    # Load tree data
    print(f"Loading tree from: {args.tree}")
    with open(args.tree, 'r') as f:
        tree_data = json.load(f)
    
    tree = tree_data['tree']
    
    # Load weights if provided
    weights = {}
    if args.weights:
        print(f"Loading weights from: {args.weights}")
        with open(args.weights, 'r') as f:
            weights = json.load(f)
    else:
        print("No weights file provided, using default weights of 1.0")
    
    # Calculate optimal scales
    print(f"key for clade annotations: {args.clade_key}")
    results,ranges = calculate_optimal_scales(tree, weights, key=args.clade_key)
    
    # Save results
    print(f"\nSaving results to: {args.output}")
    with open(args.output, 'w') as f:
        json.dump(results, f, indent=2)
    
    # Print summary
    print("\n" + "="*80)
    print("OPTIMAL SCALE PARAMETERS")
    print("="*80)
    i=0
    for param, value in results['optimal_scales'].items():
        print(f"{param:25}: {value} [{round(ranges[i][0],1)}-{round(ranges[i][1],1)}]")
        i+=1
    
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    for stat, value in results['empirical_data'].items():
        if value is not None:
            print(f"{stat:35}: {value}")
    
    print(f"\nComplete results saved to: {args.output}")


if __name__ == "__main__":
    main()