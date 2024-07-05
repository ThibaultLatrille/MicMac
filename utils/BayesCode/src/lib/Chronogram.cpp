#include "Chronogram.hpp"
#include <cassert>
#include <fstream>
#include <sstream>

using namespace std;


NodeAges::NodeAges(const Tree &intree, const AnnotatedTree &annotated_tree)
    : SimpleNodeArray<double>(intree, 0.0) {
    (*this)[this->GetTree().root()] = 0.0;
    double max_age = 0.0;
    for (auto const &node : this->GetTree().root_to_leaves_iter()) {
        if (this->GetTree().is_root(node)) {
            continue;
        }
        auto parent = this->GetTree().parent(node);
        double age_parent = this->GetVal(parent);
        double branch_length = stod(annotated_tree.tag(node, "length"));
        (*this)[node] = age_parent + branch_length;
        max_age = max(max_age, this->GetVal(node));
    }
    for (auto const &node : this->GetTree().root_to_leaves_iter()) {
        (*this)[node] = max(max_age - (*this)[node], 0.0);
    }
}

NodeAges::NodeAges(const Tree &intree, const string &fossilsfile)
    : SimpleNodeArray<double>(intree, 0.0) {
    if (fossilsfile != "Null") {
        // for line in file
        ifstream input_stream{fossilsfile};
        if (!input_stream) {
            cerr << "Fossils file " << fossilsfile << " doesn't exist" << endl;
        } else {
            cout << "Opening Fossils file " << fossilsfile << "." << endl;
        }
        string line;
        // skip the header of the file
        getline(input_stream, line);
        char sep{'\t'};
        {
            istringstream line_stream(line);
            string word;
            // Skip the first column (taxon name)
            getline(line_stream, word, sep);
            assert(word == "NodeName");
            getline(line_stream, word, sep);
            assert(word == "Age");
            getline(line_stream, word, sep);
            assert(word == "LowerBound");
            getline(line_stream, word, sep);
            assert(word == "UpperBound");
        }
        while (getline(input_stream, line)) {
            istringstream line_stream(line);
            string word;
            // Skip the first column (taxon name)
            Tree::NodeIndex node_clamped{-1};
            getline(line_stream, word, sep);
            if (word == "Root") {
                node_clamped = this->GetTree().root();
            } else {
                for (auto const &node : this->GetTree().root_to_leaves_iter()) {
                    if (this->GetTree().node_name(node) == word) { node_clamped = node; }
                }
            }
            if (node_clamped == -1) {
                cerr << word << " node is unknown." << endl;
                continue;
            }
            node_clamped_set.insert(node_clamped);
            getline(line_stream, word, sep);
            clamped_ages[node_clamped] = stod(word);
            (*this)[node_clamped] = clamped_ages[node_clamped];
            getline(line_stream, word, sep);
            clamped_lower_bound[node_clamped] = max(stod(word), 0.0);
            assert(clamped_lower_bound[node_clamped] <= clamped_ages[node_clamped]);
            getline(line_stream, word, sep);
            clamped_upper_bound[node_clamped] = stod(word);
            assert(clamped_upper_bound[node_clamped] >= clamped_ages[node_clamped]);
            string name = node_clamped == this->GetTree().root()
                              ? "Root"
                              : this->GetTree().node_name(node_clamped);
            cout << "Node " << name << " clamped to " << clamped_ages[node_clamped]
                 << " with bounds [" << clamped_lower_bound[node_clamped] << ", "
                 << clamped_upper_bound[node_clamped] << "]." << endl;
        }
        if (node_clamped_set.count(this->GetTree().root()) == 0) {
            double max_root_age = 0.0;
            double root_ecc = EccentricityRecursive(this->GetTree().root());
            for (Tree::NodeIndex const &node_clamped : node_clamped_set) {
                double d =
                    this->GetVal(node_clamped) * root_ecc / EccentricityRecursive(node_clamped);
                if (d > max_root_age) { max_root_age = d; }
            }
            (*this)[this->GetTree().root()] = max_root_age;
        }
        clamped = true;
    } else {
        (*this)[this->GetTree().root()] = 1.0;
        node_clamped_set.insert(this->GetTree().root());
        clamped_ages[this->GetTree().root()] = 1.0;
        clamped_lower_bound[this->GetTree().root()] = 1.0;
        clamped_upper_bound[this->GetTree().root()] = 1.0;
    }
    set<Tree::NodeIndex> node_init_set = node_clamped_set;
    if (node_clamped_set.count(this->GetTree().root()) == 0) {
        node_init_set.insert(this->GetTree().root());
    }
    for (Tree::NodeIndex const &node_init : node_init_set) {
        vector<Tree::NodeIndex> internal_nodes{};
        vector<Tree::NodeIndex> frontier_nodes{};
        unordered_map<Tree::NodeIndex, double> distance{{node_init, 0.0}};

        queue<Tree::NodeIndex> bfs_queue{};
        bfs_queue.push(node_init);
        while (!bfs_queue.empty()) {
            Tree::NodeIndex node_bfs = bfs_queue.front();
            bfs_queue.pop();
            if ((node_clamped_set.count(node_bfs) == 1 and node_bfs != node_init) or
                this->GetTree().is_leaf(node_bfs)) {
                frontier_nodes.push_back(node_bfs);
            } else {
                if (node_bfs != node_init) { internal_nodes.push_back(node_bfs); }
                for (Tree::NodeIndex child : this->GetTree().children(node_bfs)) {
                    distance[child] = distance[node_bfs] + 1;
                    bfs_queue.push(child);
                }
            }
        }

        // Update the partial subtree (until the clamped frontier)
        double min_delta_t = numeric_limits<double>::infinity();
        for (Tree::NodeIndex const &frontier_node : frontier_nodes) {
            double frontier_age =
                node_clamped_set.count(frontier_node) == 1 ? clamped_ages[frontier_node] : 0.0;
            assert(GetVal(node_init) > 0.0);
            double delta_t = (this->GetVal(node_init) - frontier_age) / distance[frontier_node];
            assert(delta_t > 0);
            if (delta_t < min_delta_t) { min_delta_t = delta_t; }
        }

        for (Tree::NodeIndex const &internal_node : internal_nodes) {
            auto parent = this->GetTree().parent(internal_node);
            (*this)[internal_node] = this->GetVal(parent) - min_delta_t;
            assert(GetVal(parent) > 0.0);
            assert(GetVal(internal_node) > 0.0);
            assert(GetVal(parent) > GetVal(internal_node));
        }
    }
    Check();
}

double NodeAges::EccentricityRecursive(Tree::NodeIndex node) const {
    double max_eccent = 0.0;
    if (!this->GetTree().is_leaf(node)) {
        // proceed from leaves to root using recursive algorithm
        for (auto const &child : this->GetTree().children(node)) {
            double eccent = EccentricityRecursive(child) + 1.0;
            if (eccent > max_eccent) { max_eccent = eccent; }
        }
    }
    return max_eccent;
}

bool NodeAges::Check() const {
    for (auto const &node : this->GetTree().root_to_leaves_iter()) {
        auto name = this->GetTree().node_name(node);
        if (!this->GetTree().is_root(node)) {
            if (this->GetTree().is_leaf(node)) {
                if (GetVal(node) != 0.0) {
                    cerr << "The age of the leaf " << name << " is not 0." << endl;
                }
            } else {
                if (GetVal(node) == 0.0) {
                    cerr << "The age of the internal node " << name << " is 0." << endl;
                    exit(1);
                }
            }
            auto p = this->GetTree().parent(node);
            if (GetVal(p) < GetVal(node)) {
                cerr << "The node " << name << " is older (age=" << GetVal(node)
                     << ") than it's parent node " << this->GetTree().node_name(p)
                     << " (age=" << GetVal(p) << ")." << endl;
                exit(1);
            }
        } else {
            if (GetVal(node) == 0.0) {
                cerr << "The age of root " << name << " is 0." << endl;
                exit(1);
            }
        }
    }
    return true;
}

void NodeAges::SlidingMove(Tree::NodeIndex node, double scale) {
    assert(!this->GetTree().is_leaf(node));

    double lower_bound = node_clamped_set.count(node) == 1 ? clamped_lower_bound[node] : 0.0;
    for (auto const &child : this->GetTree().children(node)) {
        if (GetVal(child) > lower_bound) { lower_bound = GetVal(child); }
    }
    double upper_bound = node_clamped_set.count(node) == 1
                             ? clamped_upper_bound[node]
                             : std::numeric_limits<double>::infinity();
    if (!this->GetTree().is_root(node)) {
        double p = GetVal(this->GetTree().parent(node));
        if (p < upper_bound) { upper_bound = p; }
    }

    if (upper_bound < GetVal(node) or lower_bound > GetVal(node)) {
        cerr << "Error node : " << this->GetTree().node_name(node) << ". Age is " << node
             << ", lower bound is " << lower_bound << " and upper bound is " << upper_bound << '.'
             << endl;
    }
    if (upper_bound == lower_bound) { return; }

    double x = GetVal(node);
    x += scale * GetVal(this->GetTree().root());

    while ((x < lower_bound) || (x > upper_bound)) {
        if (x < lower_bound) { x = 2 * lower_bound - x; }
        if (x > upper_bound) { x = 2 * upper_bound - x; }
    }

    (*this)[node] = x;
}

Chronogram::Chronogram(const NodeAges &innodeages)
    : SimpleBranchArray<double>(innodeages.GetTree()), nodeages(innodeages) {
    Update();
}

void Chronogram::Update() {
    for (auto const &node : this->GetTree().root_to_leaves_iter()) {
        if (!this->GetTree().is_root(node)) { UpdateBranch(this->GetTree().parent(node), node); }
    }
}

void Chronogram::UpdateLocal(Tree::NodeIndex node) {
    // update all branch lengths around this node

    // for the branch attached to the node
    if (this->GetTree().is_root(node)) {
        Update();
    } else {
        UpdateBranch(this->GetTree().parent(node), node);
        // for all children
        for (auto const &child : this->GetTree().children(node)) { UpdateBranch(node, child); }
    }
}

void Chronogram::UpdateBranch(Tree::NodeIndex parent, Tree::NodeIndex node) {
    (*this)[this->GetTree().branch_index(node)] =
        (nodeages.GetVal(parent) - nodeages.GetVal(node)) / nodeages.GetVal(this->GetTree().root());
    assert((*this)[this->GetTree().branch_index(node)] > 0);
}

void Chronogram::Scale() {
    for (auto const &node : this->GetTree().root_to_leaves_iter()) {
        if (!this->GetTree().is_root(node)) {
            (*this)[this->GetTree().branch_index(node)] *= nodeages.GetVal(this->GetTree().root());
        }
    }
}
