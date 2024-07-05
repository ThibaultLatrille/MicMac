#pragma once

#include "BranchProduct.hpp"
#include "Chronogram.hpp"
#include "CodonSequenceAlignment.hpp"
#include "CodonSubMatrixBranchArray.hpp"
#include "CodonSuffStat.hpp"
#include "GTRSubMatrix.hpp"
#include "GammaSuffStat.hpp"
#include "IIDGamma.hpp"
#include "Move.hpp"
#include "MultivariateProcess.hpp"
#include "PhyloProcess.hpp"
#include "ScatterSuffStat.hpp"
#include "TaxonTraits.hpp"
#include "components/ChainComponent.hpp"
#include "components/Tracer.hpp"

/**
 * \brief A site-homogeneous and branch-heterogenous Muse and Gaut omega-codon model.
 * Mutation rate and omega are modeled as a browian process along the tree (like CoEvol)
 *
 * The model has the following structure:
 * - nucleotide relative exchangeabilities and stationaries are uniform
 * Dirichlet
 * - there is one omega=dN/dS for each branches but for all sites
 * - branch rates (mutation rates) and omega are brownian multivariate
 */


class DatedNodeOmegaModel : public ChainComponent {
    // tree and data
    std::string datafile, treefile, traitsfile{"Null"}, fossilsfile{"Null"};
    std::unique_ptr<const Tree> tree;
    FileSequenceAlignment *data;
    const TaxonSet *taxonset;
    CodonSequenceAlignment *codondata;
    TaxonTraits *taxon_traits{nullptr};

    int Nsite;
    int Ntaxa;
    int Nbranch;

    // Node ages
    NodeAges *nodeages;
    // Chronogram (diff between node ages)
    Chronogram *chronogram;

    // Precision matrix and prior
    int dimensions;
    int prior_cov_df;
    bool uniq_kappa;
    PriorCovariance *prior_matrix;
    PrecisionMatrix *precision_matrix;
    NodeMultivariateProcess *node_multivariate;

    // Branch rates (brownian process)
    NodeProcess *noderates;
    BranchProcess *branchrates;

    // Branch lengths (product of branch rates and chronogram)
    BranchwiseProduct *branchlength;

    // Branch omega (brownian process)
    NodeProcess *nodeomega;
    BranchProcess *branchomega;

    // Nucleotide rates
    std::vector<double> nucrelratehypercenter;
    double nucrelratehyperinvconc;
    std::vector<double> nucstathypercenter;
    double nucstathyperinvconc;

    std::vector<double> nucrelrate;
    std::vector<double> nucstat;
    GTRSubMatrix *nucmatrix;

    // codon matrices (parameterized by nucmatrix and omega)
    MGOmegaCodonSubMatrixBranchArray *codonmatrixbrancharray;
    MGOmegaCodonSubMatrix *rootcodonmatrix;

    // PhyloProcess
    PhyloProcess *phyloprocess;

    // suff stats for substitution paths
    // summed over all branches and over all sites
    PathSuffStatNodeArray *pathsuffstatarray;

    // path suff stat can be summarized in terms of 4x4 suff stats, as a function
    // of nucleotide rates
    NucPathSuffStat nucpathsuffstat;

    // or, alternatively, collected as a simple Poisson suff stat, as a function
    // of nodeomega
    OmegaPathSuffStatBranchArray *omegapathsuffstatarray;

    // Poisson suffstats for substitution histories, as a function of branch
    // lengths
    PoissonSuffStatBranchArray *lengthpathsuffstatarray;

    ScatterSuffStat *scattersuffstat;

  public:
    friend std::ostream &operator<<(std::ostream &os, DatedNodeOmegaModel &m);

    //-------------------
    // Construction and allocation
    // ------------------

    //! \brief constructor, parameterized by names of data and tree files
    //!
    //! Note: in itself, the constructor does not allocate the model;
    //! It only reads the data and tree file and register them together.
    DatedNodeOmegaModel(std::string indatafile, std::string intreefile, std::string intraitsfile,
        std::string infossilsfile, int inprior_cov_df, bool inuniq_kappa)
        : datafile(std::move(indatafile)),
          treefile(std::move(intreefile)),
          traitsfile(std::move(intraitsfile)),
          fossilsfile(std::move(infossilsfile)),
          prior_cov_df{inprior_cov_df},
          uniq_kappa{inuniq_kappa} {
        data = new FileSequenceAlignment(datafile);
        codondata = new CodonSequenceAlignment(data, true);

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        taxonset = codondata->GetTaxonSet();

        if (traitsfile != "Null") { taxon_traits = new TaxonTraits(traitsfile, *taxonset, false); }

        // get tree from file (newick format)
        std::ifstream tree_stream{treefile};
        NHXParser parser{tree_stream};
        tree = make_from_parser(parser);
        Nbranch = tree->nb_branches();

        // Node ages
        nodeages = new NodeAges(*tree, fossilsfile);
        // Chronogram (diff between node ages)
        chronogram = new Chronogram(*nodeages);

        dimensions = 2;
        if (taxon_traits != nullptr) { dimensions += taxon_traits->GetDim(); }
        prior_matrix = new PriorCovariance(dimensions, prior_cov_df, uniq_kappa);
        precision_matrix = new PrecisionMatrix(*prior_matrix);

        node_multivariate = new NodeMultivariateProcess(*chronogram, *precision_matrix, dimensions);

        // Branch omega (brownian process)
        nodeomega = new NodeProcess(*node_multivariate, dim_omega);
        branchomega = new BranchProcess(*nodeomega);

        // Branch rates (brownian process)
        noderates = new NodeProcess(*node_multivariate, dim_mut_rate);
        branchrates = new BranchProcess(*noderates);

        // Branch lengths (product of branch rates and chronogram)
        branchlength = new BranchwiseProduct(*chronogram, *branchrates);

        // Nucleotide rates
        nucrelratehypercenter.assign(Nrr, 1.0 / Nrr);
        nucrelratehyperinvconc = 1.0 / Nrr;
        nucrelrate.assign(Nrr, 0);
        Random::DirichletSample(nucrelrate, nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);

        nucstathypercenter.assign(Nnuc, 1.0 / Nnuc);
        nucstathyperinvconc = 1.0 / Nnuc;
        nucstat.assign(Nnuc, 0);
        Random::DirichletSample(nucstat, nucstathypercenter, 1.0 / nucstathyperinvconc);

        nucmatrix = new GTRSubMatrix(Nnuc, nucrelrate, nucstat, true);

        // Codon Matrices
        codonmatrixbrancharray =
            new MGOmegaCodonSubMatrixBranchArray(GetCodonStateSpace(), nucmatrix, branchomega);
        rootcodonmatrix = new MGOmegaCodonSubMatrix(
            GetCodonStateSpace(), nucmatrix, nodeomega->GetExpVal(tree->root()));

        // PhyloProcess
        phyloprocess = new PhyloProcess(
            tree.get(), codondata, branchlength, nullptr, codonmatrixbrancharray, rootcodonmatrix);
        phyloprocess->Unfold();

        if (taxon_traits != nullptr) {
            node_multivariate->ClampLeaves(*taxon_traits, phyloprocess->GetTaxonMap());
        }
        // Suff Stats
        lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);
        omegapathsuffstatarray = new OmegaPathSuffStatBranchArray(*tree);
        pathsuffstatarray = new PathSuffStatNodeArray(*tree);
        scattersuffstat = new ScatterSuffStat(*tree);
    }

    virtual ~DatedNodeOmegaModel() = default;

    void move(int it) override { Move(); }

    template <class Info>
    void declare_interface(Info info) {
        model_node(info, "nucstat", nucstat);
        model_node(info, "nucrelrate", nucrelrate);
        model_node(info, "nodeages", *nodeages);
        model_node(info, "node_multivariate", *node_multivariate);
        model_node(info, "prior_cov_matrix", *prior_matrix);
        model_node(info, "precision_matrix", *precision_matrix);

        model_stat(info, "logprior", [&]() { return DatedNodeOmegaModel::GetLogPrior(); });
        model_stat(info, "lnL", [&]() { return DatedNodeOmegaModel::GetLogLikelihood(); });
        model_stat(info, "blengthmean", [&]() { return branchlength->GetMean(); });
        model_stat(info, "blengthvar", [&]() { return branchlength->GetVar(); });
        for (int i = 0; i < dimensions; i++) {
            model_stat(info, "PriorCovariance_" + std::to_string(i), prior_matrix->coeffRef(i));
            for (int j = 0; j <= i; j++) {
                model_stat(info, "Precision_" + std::to_string(i) + "_" + std::to_string(j),
                    precision_matrix->coeffRef(i, j));
            }
        }
        // Descriptive statistics - for each branch of the tree
        for (Tree::BranchIndex branch = 0; branch < tree->nb_branches(); branch++) {
            string b_name = tree->node_name(tree->node_index(branch));
            model_stat(info, "*BranchTime_" + b_name, (*chronogram)[branch]);
            model_stat(info, "*BranchMutRate_" + b_name, (*branchrates)[branch]);
            model_stat(info, "*BranchLength_" + b_name, (*branchlength)[branch]);
            model_stat(info, "*BranchdNdS_" + b_name, (*branchomega)[branch]);
        }

        model_stat(info, "PredictedDNDS", [this]() { return branchomega->GetMean(); });
        model_stat(info, "statent", [&]() { return Random::GetEntropy(nucstat); });
        model_stat(info, "rrent", [&]() { return Random::GetEntropy(nucrelrate); });
    }

    //-------------------
    // Accessors
    // ------------------

    //! const access to codon state space
    CodonStateSpace *GetCodonStateSpace() const {
        return (CodonStateSpace *)codondata->GetStateSpace();
    }

    //! return number of aligned sites
    int GetNsite() const { return Nsite; }

    //! return tree
    const Tree &GetTree() const { return *tree; }

    //! return time of branch
    double GetBranchTime(Tree::NodeIndex node) const {
        assert(!tree->is_root(node));
        return chronogram->GetVal(tree->branch_index(node));
    };

    //! return branch length
    double GetBranchLength(Tree::NodeIndex node) const {
        assert(!tree->is_root(node));
        return branchlength->GetVal(tree->branch_index(node));
    };

    //! return precision matrix
    PrecisionMatrix GetPrecisionMatrix() const { return *precision_matrix; };


    //! return the value of the multivariate brownian process for a given node and a given
    //! dimensions of the process
    double GetExpBrownianEntry(Tree::NodeIndex node, int dim) const {
        return exp(node_multivariate->GetVal(node)(dim));
    }

    //! return number of dimensions of the multivariate brownian process
    int GetDimension() const { return dimensions; }

    std::string GetDimensionName(int dim) const {
        if (dim == dim_omega) {
            return "Omega";
        } else if (dim == dim_mut_rate) {
            return "MutationRatePerTime";
        } else {
            assert(taxon_traits != nullptr);
            return taxon_traits->GetHeader(dim);
        }
    }

    //! \brief tell the nucleotide matrix that its parameters have changed and
    //! that it should be updated
    //!
    //! The matrix is not directly updated at that step. Instead, corruption is
    //! notified, such that the matrix knows that it will have to recalculate
    //! whichever component is requested later on upon demand.
    void TouchNucMatrix() {
        nucmatrix->CopyStationary(nucstat);
        nucmatrix->CorruptMatrix();
    }

    //! \brief tell the codon matrix that its parameters have changed and that it
    //! should be updated
    //!
    //! The matrix is not directly updated at that step. Instead, corruption is
    //! notified, such that the matrix knows that it will have to recalculate
    //! whichever component is requested later on upon demand.
    void TouchCodonMatrices() {
        codonmatrixbrancharray->UpdateCodonMatrices();
        rootcodonmatrix->SetOmega(nodeomega->GetExpVal(tree->root()));
        rootcodonmatrix->CorruptMatrix();
    }

    //! \brief tell the nucleotide and the codon matrices that their parameters
    //! have changed and that they should be updated
    //!
    //! Just successive calls to TouchNucMatrix() and then TouchCodonMatrix();
    void TouchMatrices() {
        TouchNucMatrix();
        TouchCodonMatrices();
    }

    //! \brief dummy function that does not do anything.
    //!
    //! Used for the templates of ScalingMove, SlidingMove and ProfileMove
    //! (defined in Move), all of which require a void (*f)(void) function
    //! pointer to be called after changing the value of the focal parameter.
    void NoUpdate() {}

    //! \brief Update the BranchProcess (rates and omega) with the underlying NodeProcess, Update
    //! the Chronogram with the underlying NodeAges. And finally update the branch lengths with the
    //! Chronogram and the BranchProcess (rates).
    //!
    //! Used when the model is restarted or for the posterior predictif.
    void UpdateBranches(bool scale = false) {
        chronogram->Update();
        branchomega->Update();
        branchrates->Update();
        branchlength->Update();
        if (scale) {chronogram->Scale();}
    }

    //! \brief Update the chronogram (branch time) and branch lengths around the focal node.
    //!
    //! Update needed when the age (NodeAges) of the focal node is changed.
    void UpdateLocalChronogram(Tree::NodeIndex node) {
        chronogram->UpdateLocal(node);
        branchlength->UpdateLocal(node);
    }

    //! \brief Update the branch rates and lengths around the focal node.
    //!
    //! Update needed when the rate (NodeProcess) of the focal node is changed.
    void UpdateLocalBranchRates(Tree::NodeIndex node) {
        branchrates->UpdateLocal(node);
        branchlength->UpdateLocal(node);
    }

    //! \brief Update the branch omega for a focal node.
    //!
    //! Update needed when the omega (NodeProcess) of the focal node is changed.
    void UpdateBranchOmega(Tree::NodeIndex node) {
        if (tree->is_root(node)) {
            rootcodonmatrix->SetOmega(nodeomega->GetExpVal(tree->root()));
            rootcodonmatrix->CorruptMatrix();
        } else {
            Tree::BranchIndex branch = tree->branch_index(node);
            (*codonmatrixbrancharray)[branch].SetOmega(branchomega->GetVal(branch));
            (*codonmatrixbrancharray)[branch].CorruptMatrix();
        }
    }

    //! \brief Update the branch omega around the focal node.
    //!
    //! Update needed when the omega (NodeProcess) of the focal node is changed.
    void UpdateLocalBranchOmega(Tree::NodeIndex node) {
        branchomega->UpdateLocal(node);
        UpdateBranchOmega(node);
        for (Tree::NodeIndex const &child : tree->children(node)) { UpdateBranchOmega(child); }
    }

    //! \brief global update function (includes the stochastic mapping of
    //! character history)
    void Update() {
        UpdateBranches();
        TouchMatrices();
    }

    //-------------------
    // Posterior Predictive
    // ------------------

    //! \brief post pred function (does the update of all fields before doing the
    //! simulation)
    void PostPred(std::string name) {
        Update();
        phyloprocess->PostPredSample(name);
    }

    //-------------------
    // Priors and likelihood
    //-------------------

    //! \brief return total log prior
    //!
    //! Note: up to some multiplicative constant
    double GetLogPrior() const {
        double total = PrecisionMatrixLogProb();
        total += NodeMultivariateLogPrior();
        total += NucRatesLogPrior();
        return total;
    }

    //! return current value of likelihood (pruning-style, i.e. integrated over
    //! all substitution histories)
    double GetLogLikelihood() const { return phyloprocess->GetLogLikelihood(); }

    //! return joint log prob (log prior + log likelihood)
    double GetLogProb() const { return GetLogPrior() + GetLogLikelihood(); }

    //! log prob of precision matrix
    double PrecisionMatrixLogProb() const {
        return prior_matrix->GetLogProb() + precision_matrix->GetLogProb(*prior_matrix);
    }

    //! log prior over branch rates (brownian process)
    double NodeMultivariateLogPrior() const { return node_multivariate->GetLogProb(); }

    //! log prior of branch rate (brownian process) around of focal node
    double LocalNodeMultivariateLogPrior(Tree::NodeIndex node) const {
        return node_multivariate->GetLocalLogProb(node);
    }

    // Nucleotide rates prior
    //! log prior over nucleotide relative exchangeabilities (nucrelrate) and eq.
    //! freqs. (nucstat) -- uniform Dirichlet in both cases
    double NucRatesLogPrior() const {
        double total = 0;
        total += Random::logDirichletDensity(
            nucrelrate, nucrelratehypercenter, 1.0 / nucrelratehyperinvconc);
        total +=
            Random::logDirichletDensity(nucstat, nucstathypercenter, 1.0 / nucstathyperinvconc);
        return total;
    }

    //-------------------
    // Collect Suff Stat
    //-------------------

    // Paths
    //! collect generic sufficient statistics from substitution mappings
    void CollectPathSuffStat() {
        pathsuffstatarray->Clear();
        pathsuffstatarray->AddSuffStat(*phyloprocess);
    }

    // Branch lengths
    //! collect sufficient statistics for moving branch lengths
    void CollectLengthSuffStat() {
        lengthpathsuffstatarray->Clear();
        lengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
    }

    // Nucleotide rates
    //! collect sufficient statistics for moving nucleotide rates (based on
    //! generic sufficient statistics stored in pathsuffstat)
    void CollectNucPathSuffStat() {
        nucpathsuffstat.Clear();
        nucpathsuffstat.AddSuffStat(*codonmatrixbrancharray, *rootcodonmatrix, *pathsuffstatarray);
    }

    // Omega
    //! collect sufficient statistics for moving omega (based on
    //! generic sufficient statistics stored in pathsuffstat)
    void CollectOmegaPathSuffStat() {
        omegapathsuffstatarray->Clear();
        omegapathsuffstatarray->AddSuffStat(
            *codonmatrixbrancharray, *rootcodonmatrix, *pathsuffstatarray);
    }

    // Scatter (brownian process)
    void CollectScatterSuffStat() {
        scattersuffstat->Clear();
        scattersuffstat->AddSuffStat(*node_multivariate);
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    // Node ages and branch rates
    //! \brief log prob to be recomputed when moving age of focal node
    double LocalNodeAgeLogProb(Tree::NodeIndex node) const {
        if (GetTree().is_root(node)) {
            return NodeMultivariateLogPrior() + BranchLengthSuffStatLogProb();
        } else {
            return LocalNodeMultivariateLogPrior(node) + LocalBranchLengthSuffStatLogProb(node);
        }
    }

    //! \brief log prob to be recomputed when moving branch rates (brownian process) around of focal
    //! node
    double LocalNodeRatesLogProb(Tree::NodeIndex node) const {
        return LocalNodeMultivariateLogPrior(node) + LocalBranchLengthSuffStatLogProb(node);
    }

    //! \brief log prob factor (without prior) to be recomputed when moving age of focal node, or
    //! when moving branch rates (brownian process) around of focal node.
    double LocalBranchLengthSuffStatLogProb(Tree::NodeIndex node) const {
        double tot = 0;
        // for all children
        for (auto const &child : tree->children(node)) {
            tot += BranchLengthSuffStatLogProb(tree->branch_index(child));
        }
        if (!tree->is_root(node)) {
            // for the branch attached to the node
            tot += BranchLengthSuffStatLogProb(tree->branch_index(node));
        }
        assert(tot != 0);
        return tot;
    }

    //! \brief return log prob of current substitution mapping (on focal branch), as a function of
    //! the length of a given branch
    double BranchLengthSuffStatLogProb(Tree::BranchIndex branch) const {
        return lengthpathsuffstatarray->GetVal(branch).GetLogProb(branchlength->GetVal(branch));
    }

    //! \brief return log prob of current substitution mapping (on all branches), as a function of
    //! the length all branches
    double BranchLengthSuffStatLogProb() const {
        return lengthpathsuffstatarray->GetLogProb(*branchlength);
    }

    // Omega
    //! \brief log prob to be recomputed when moving omega (brownian process) around of focal node
    double LocalNodeOmegaLogProb(Tree::NodeIndex node) const {
        return LocalNodeMultivariateLogPrior(node) + LocalNodeOmegaSuffStatLogProb(node);
    }

    //! \brief log prob factor (without prior) to be recomputed when moving omega (brownian process)
    //! around of focal node
    double LocalNodeOmegaSuffStatLogProb(Tree::NodeIndex node) const {
        double tot = 0;
        // for all children
        for (auto const &child : tree->children(node)) {
            tot += NodeOmegaSuffStatLogProb(tree->branch_index(child));
        }
        if (!tree->is_root(node)) {
            // for the branch attached to the node
            tot += NodeOmegaSuffStatLogProb(tree->branch_index(node));
        }
        assert(tot != 0);
        return tot;
    }

    //! \brief return log prob of current substitution mapping (on focal branch), as a function of
    //! omega of a given branch
    double NodeOmegaSuffStatLogProb(Tree::BranchIndex branch) const {
        return omegapathsuffstatarray->GetVal(branch).GetLogProb(branchomega->GetVal(branch));
    }

    // Nucleotide rates
    //! \brief log prob factor to be recomputed when moving nucleotide mutation
    //! rate parameters (nucrelrate and nucstat)
    //! Calculated using nucpathsuffstat
    //! (which summarizes all information about how the probability of the
    //! substitution mapping depends on nucleotide mutation rates) and the
    //! nucmatrix. Both nucpathsuffstat and nucmatrix are assumed to be updated.
    double NucRatesLogProb() const {
        return NucRatesLogPrior() + nucpathsuffstat.GetLogProb(*nucmatrix, *GetCodonStateSpace());
    }

    //-------------------
    //  Moves
    //-------------------

    //! \brief complete MCMC move schedule
    double Move() {
        ResampleSub(1.0);
        MoveParameters(30);
        return 1.0;
    }

    //! Gibbs resample substitution mappings conditional on current parameter
    //! configuration
    void ResampleSub(double frac) {
        TouchMatrices();
        phyloprocess->Move(frac);
    }

    //! complete series of MCMC moves on all parameters (repeated nrep times)
    void MoveParameters(int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            CollectLengthSuffStat();
            MoveNodeAges(0.1, 3);
            MoveNodeAges(0.02, 3);

            MoveNodeRates(1.0, 3);
            MoveNodeRates(0.1, 3);

            CollectPathSuffStat();
            CollectNucPathSuffStat();
            MoveNucRates();
            TouchMatrices();

            CollectOmegaPathSuffStat();
            MoveNodeOmega(1.0, 3);
            MoveNodeOmega(0.1, 3);

            MoveNodeTraits(0.5, 3, true);
            MoveNodeTraits(0.5, 3, false);
            MoveNodeTraits(0.05, 3, true);
            MoveNodeTraits(0.05, 3, false);

            CollectScatterSuffStat();
            SamplePrecisionMatrix();
            MovePriorMatrix(0.1, 3);
            MovePriorMatrix(0.01, 3);
            TouchMatrices();
        }
    }

    void SamplePrecisionMatrix() {
        scattersuffstat->SamplePrecisionMatrix(*precision_matrix, *prior_matrix);
    };

    //! MH moves on the invert wishart matrix (prior of the covariance matrix)
    void MovePriorMatrix(double tuning, int nrep) {
        if (uniq_kappa) {
            for (int rep = 0; rep < nrep; rep++) {
                double deltalogprob = -PrecisionMatrixLogProb();
                double m = tuning * (Random::Uniform() - 0.5);
                double e = exp(m);
                (*prior_matrix) *= e;
                deltalogprob += PrecisionMatrixLogProb();
                deltalogprob += m;
                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (!accepted) { (*prior_matrix) /= e; }
            }
        } else {
            for (int i = 0; i < dimensions; ++i) {
                Move::Scaling((*prior_matrix)(i), tuning, nrep,
                    &DatedNodeOmegaModel::PrecisionMatrixLogProb, &DatedNodeOmegaModel::NoUpdate,
                    this);
            }
        }
    }

    //! MH moves on branch ages
    void MoveNodeAges(double tuning, int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            for (Tree::NodeIndex node : tree->root_to_leaves_iter()) {
                if (!tree->is_leaf(node) and !(nodeages->Unclamped() and tree->is_root(node))) {
                    MoveNodeAge(node, tuning);
                }
            }
        }
    }

    //! MH moves on branch ages for a focal node
    void MoveNodeAge(Tree::NodeIndex node, double tuning) {
        double logratio = -LocalNodeAgeLogProb(node);

        double bk = nodeages->GetVal(node);
        double sliding = tuning * (Random::Uniform() - 0.5);
        nodeages->SlidingMove(node, sliding);
        UpdateLocalChronogram(node);

        logratio += LocalNodeAgeLogProb(node);

        bool accept = (log(Random::Uniform()) < logratio);
        if (!accept) {
            (*nodeages)[node] = bk;
            UpdateLocalChronogram(node);
        }
    }

    //! MH moves on branch rates (brownian process)
    void MoveNodeRates(double tuning, int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            for (Tree::NodeIndex node : tree->root_to_leaves_iter()) {
                MoveNodeProcessRate(node, tuning);
            }
        }
    }

    //! MH moves on branch rates (brownian process) for a focal node
    void MoveNodeProcessRate(Tree::NodeIndex node, double tuning) {
        double logratio = -LocalNodeRatesLogProb(node);

        double m = tuning * noderates->GetSigma() * (Random::Uniform() - 0.5);
        noderates->SlidingMove(node, m);
        UpdateLocalBranchRates(node);

        logratio += LocalNodeRatesLogProb(node);

        bool accept = (log(Random::Uniform()) < logratio);
        if (!accept) {
            noderates->SlidingMove(node, -m);
            UpdateLocalBranchRates(node);
        }
    }

    //! MH moves on branch omega (brownian process)
    void MoveNodeOmega(double tuning, int nrep) {
        for (int rep = 0; rep < nrep; rep++) {
            for (Tree::NodeIndex node : tree->root_to_leaves_iter()) {
                MoveNodeOmega(node, tuning);
            }
        }
    }

    //! MH moves on branch omega (brownian process) for a focal node
    void MoveNodeOmega(Tree::NodeIndex node, double tuning) {
        double logratio = -LocalNodeOmegaLogProb(node);

        double m = tuning * nodeomega->GetSigma() * (Random::Uniform() - 0.5);
        nodeomega->SlidingMove(node, m);
        UpdateLocalBranchOmega(node);

        logratio += LocalNodeOmegaLogProb(node);

        bool accept = (log(Random::Uniform()) < logratio);
        if (!accept) {
            nodeomega->SlidingMove(node, -m);
            UpdateLocalBranchOmega(node);
        }
    }


    //! MH moves on traits (brownian process)
    void MoveNodeTraits(double tuning, int nrep, bool leaves_to_root) {
        if (taxon_traits == nullptr or taxon_traits->GetDim() == 0) { return; }
        for (int rep = 0; rep < nrep; rep++) {
            for (int trait_dim = 0; trait_dim < taxon_traits->GetDim(); trait_dim++) {
                for (Tree::NodeIndex node :
                    leaves_to_root ? tree->leaves_root_to_iter() : tree->root_to_leaves_iter()) {
                    if (!tree->is_leaf(node) or
                        !taxon_traits->DataPresence(
                            phyloprocess->GetTaxonMap().NodeToTaxon(node), trait_dim)) {
                        MoveNodeTrait(node, tuning, trait_dim);
                    }
                }
            }
        }
    }

    //! MH moves on traits (brownian process) for a focal node
    void MoveNodeTrait(Tree::NodeIndex node, double tuning, int trait_dim) {
        int multi_dim = taxon_traits->TraitDimToMultivariateDim(trait_dim);
        double logratio = -LocalNodeMultivariateLogPrior(node);

        double m = tuning * (Random::Uniform() - 0.5);
        (*node_multivariate)[node](multi_dim) += m;

        logratio += LocalNodeMultivariateLogPrior(node);

        bool accept = (log(Random::Uniform()) < logratio);
        if (!accept) { (*node_multivariate)[node](multi_dim) -= m; }
    }

    // Nucleotide rates
    //! MH moves on nucleotide rate parameters (nucrelrate and nucstat: using
    //! ProfileMove)
    void MoveNucRates() {
        Move::Profile(nucrelrate, 0.1, 1, 3, &DatedNodeOmegaModel::NucRatesLogProb,
            &DatedNodeOmegaModel::TouchNucMatrix, this);
        Move::Profile(nucrelrate, 0.03, 3, 3, &DatedNodeOmegaModel::NucRatesLogProb,
            &DatedNodeOmegaModel::TouchNucMatrix, this);
        Move::Profile(nucrelrate, 0.01, 3, 3, &DatedNodeOmegaModel::NucRatesLogProb,
            &DatedNodeOmegaModel::TouchNucMatrix, this);
        Move::Profile(nucstat, 0.1, 1, 3, &DatedNodeOmegaModel::NucRatesLogProb,
            &DatedNodeOmegaModel::TouchNucMatrix, this);
        Move::Profile(nucstat, 0.01, 1, 3, &DatedNodeOmegaModel::NucRatesLogProb,
            &DatedNodeOmegaModel::TouchNucMatrix, this);
    }

    void ToStream(std::ostream &os) { os << *this; }
};

std::istream &operator>>(std::istream &is, std::unique_ptr<DatedNodeOmegaModel> &m) {
    std::string model_name, datafile, treefile, traitsfile, fossilsfile;
    int prior_cov_df;
    bool uniq_kappa;

    is >> model_name;
    if (model_name != "DatedNodeOmega") {
        std::cerr << "Expected DatedNodeOmega for model name, got " << model_name << "\n";
        exit(1);
    }

    is >> datafile >> treefile >> traitsfile >> fossilsfile >> prior_cov_df >> uniq_kappa;
    m = std::make_unique<DatedNodeOmegaModel>(
        datafile, treefile, traitsfile, fossilsfile, prior_cov_df, uniq_kappa);
    Tracer tracer{*m};
    tracer.read_line(is);
    m->Update();
    return is;
}

std::ostream &operator<<(std::ostream &os, DatedNodeOmegaModel &m) {
    Tracer tracer{m};
    os << "DatedNodeOmega" << '\t';
    os << m.datafile << '\t';
    os << m.treefile << '\t';
    assert(!m.traitsfile.empty());
    os << m.traitsfile << '\t';
    assert(!m.fossilsfile.empty());
    os << m.fossilsfile << '\t';
    os << m.prior_cov_df << '\t';
    os << m.uniq_kappa << '\t';
    tracer.write_line(os);
    return os;
}
