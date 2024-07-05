#pragma once

#include <fstream>
#include <sstream>
#include "Random.hpp"
#include "TaxonMapping.hpp"

class TaxonTraits {
  public:
    //! constructor, based on a vector of taxon names
    explicit TaxonTraits(std::string const &filename, TaxonSet const &taxon_set,
        bool polymorphism_aware, int constrained_dim = 2)
        : constrained_dimensions{constrained_dim} {
        taxon_presence.resize(taxon_set.GetNtaxa(), false);
        int gentime_dim{-1};

        // for line in file
        std::ifstream input_stream{filename};
        if (!input_stream) {
            std::cerr << "Traits file " << filename << " doesn't exist" << std::endl;
        } else {
            std::cout << "Opening traits file " << filename << "." << std::endl;
        }

        std::string line;
        // skip the header of the file
        getline(input_stream, line);
        char sep{'\t'};
        {
            std::istringstream line_stream(line);
            std::string word;
            // Skip the first column (taxon name)
            getline(line_stream, word, sep);
            assert(word == "TaxonName" or word == "Name");
            int counter = 0;
            while (getline(line_stream, word, sep)) {
                if (word == "LogGenerationTime" and polymorphism_aware) {
                    gentime = true;
                    taxon_gentime_presence.resize(taxon_set.GetNtaxa(), false);
                    taxon_gentime.resize(taxon_set.GetNtaxa(), 0.0);
                    gentime_dim = counter;
                } else {
                    header.push_back(word);
                }
                counter++;
            }
        }
        dimensions = header.size();
        if (!header.empty()) {
            taxon_traits_presence.resize(taxon_set.GetNtaxa());
            taxon_traits.resize(taxon_set.GetNtaxa());
        }
        while (getline(input_stream, line)) {
            std::istringstream line_stream(line);
            std::string taxon{};
            getline(line_stream, taxon, sep);
            if (!taxon_set.TaxonPresent(taxon)) {
                std::cerr << "Taxon " << taxon
                          << " found in traits file (.tsv) does not have a match (skipping line)."
                          << std::endl;
                continue;
            }
            int id = taxon_set.GetTaxonIndex(taxon);
            taxon_presence.at(id) = true;
            if (!header.empty()) {
                taxon_traits_presence[id] = std::vector<bool>(dimensions, false);
                taxon_traits[id] = std::vector<double>(dimensions, 0.0);
            }
            std::string word;
            int counter = 0;
            int dim_counter = 0;
            while (getline(line_stream, word, sep)) {
                if (!word.empty() and convertible_to_float(word)) {
                    double val = std::stod(word);
                    if (!std::isnan(val)) {
                        if (counter == gentime_dim) {
                            taxon_gentime_presence.at(id) = true;
                            taxon_gentime.at(id) = val;
                        } else {
                            taxon_traits_presence.at(id).at(dim_counter) = true;
                            taxon_traits.at(id).at(dim_counter) = val;
                        }
                    }
                }
                if (counter != gentime_dim) { dim_counter++; }
                counter++;
            }
        }
        int nbr_presence{0};
        for (bool p : taxon_presence) { nbr_presence += p; }
        if (nbr_presence == 0) {
            std::cerr << "No taxa was found in trait file. Either it is empty or the names don't "
                         "match the names in the alignment."
                      << std::endl;
        } else {
            std::cerr << nbr_presence << " taxa found in trait file matched name the alignment."
                      << std::endl;
        }
        assert(polymorphism_aware == gentime);
    }

    bool static convertible_to_float(const string &in) {
        std::stringstream sstr(in);
        float f;
        return bool(sstr >> f);
    }

    int GetDim() const { return dimensions; }

    int TraitDimToMultivariateDim(int trait_dim) const {
        return constrained_dimensions + gentime + trait_dim;
    }

    int MultivariateDimToTraitDim(int multi_dim) const {
        return multi_dim - constrained_dimensions - gentime;
    }

    string GetHeader(int multi_dim) { return header.at(MultivariateDimToTraitDim(multi_dim)); }

    bool DataPresence(int taxon, int dim) const {
        if (taxon_presence[taxon]) {
            return taxon_traits_presence[taxon].at(dim);
        } else {
            return false;
        }
    }

    double Data(int taxon, int dim) const {
        assert(taxon_traits_presence[taxon][dim]);
        return taxon_traits[taxon].at(dim);
    }

    bool GenTimePresence() const { return gentime; }

    bool GenTimePresence(int taxon) const { return taxon_gentime_presence[taxon]; }

    double GenTime(int taxon) const {
        assert(taxon_gentime_presence[taxon]);
        return taxon_gentime[taxon];
    }

    //! default constructor
    ~TaxonTraits() = default;

  private:
    int dimensions{0};
    int constrained_dimensions;
    std::vector<std::string> header;
    std::vector<bool> taxon_presence;
    std::vector<std::vector<bool>> taxon_traits_presence;
    std::vector<std::vector<double>> taxon_traits;

    bool gentime{false};
    std::vector<bool> taxon_gentime_presence;
    std::vector<double> taxon_gentime;
};
