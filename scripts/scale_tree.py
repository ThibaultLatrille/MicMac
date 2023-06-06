import argparse
from neutrality_index import open_tree, prune_tree
from pre_processed_traits_mammals import rename_tree


def main(primate_tree, mammal_neutral_tree, output_tree):
    """
    Because the primate tree is not based on neutral markers, we scaled it so that the total branch length of the
    primate tree and the mammal tree (based on neutral markers) are the same when both trees are pruned to the same
    intersection of taxa.
    :param primate_tree:
    :param mammal_neutral_tree:
    :param output_tree:
    :return:
    """
    t_primate = rename_tree(open_tree(primate_tree, format_ete3=1))
    discard_list = ["Galeopterus_variegatus", "Tupaia_belangeri", "Oryctolagus_cuniculus", "Mus_musculus"]
    keep_species = list(set(t_primate.get_leaf_names()).difference(set(discard_list)))
    t_primate = prune_tree(t_primate, keep_species)
    t_mammal = rename_tree(open_tree(mammal_neutral_tree, format_ete3=1))
    intersection = list(set(t_primate.get_leaf_names()).intersection(set(t_mammal.get_leaf_names())))
    print(f"Intersection between mammals and primates has {len(intersection)} leaves")
    pruned_primate = prune_tree(t_primate, intersection)
    t_mammal = prune_tree(t_mammal, intersection)
    assert set(pruned_primate.get_leaf_names()) == set(t_mammal.get_leaf_names())
    total_mammal_distance = sum([n.dist for n in t_mammal.traverse()])
    total_primate_distance = sum([n.dist for n in pruned_primate.traverse()])
    ratio = total_mammal_distance / total_primate_distance
    print(f"Total mammal distance: {total_mammal_distance}")
    print(f"Total primate distance: {total_primate_distance}")
    print(f"Ratio: {ratio}")
    # Scale the primate tree
    for n in t_primate.traverse():
        n.dist = n.dist * ratio
    # Scale the pruned primate tree
    for n in pruned_primate.traverse():
        n.dist = n.dist * ratio

    distance_scaled = sum([n.dist for n in pruned_primate.traverse()])
    assert abs(distance_scaled - total_mammal_distance) < 1e-6, f"{distance_scaled} != {total_mammal_distance}"
    t_primate.write(outfile=output_tree, format=1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--primate_tree", help="Input primate tree file", required=True)
    parser.add_argument("--mammal_neutral_tree", help="Input mammal neutral tree file", required=True)
    parser.add_argument("--output_tree", help="Output scaled primate tree file", required=True)
    args = parser.parse_args()
    main(args.primate_tree, args.mammal_neutral_tree, args.output_tree)
